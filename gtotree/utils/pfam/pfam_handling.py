import os
import sys
from gtotree.utils.misc.general import (concat_files, decode_pyhmmer_text,
                                        read_single_column_file)
from gtotree.utils.misc.messaging import wprint, color_text, spinner
from gtotree.utils.pfam.get_pfam_data import HMM_FILENAME, INFO_FILENAME
import pyhmmer  # type: ignore


def _check_pfam_data_dir():
    try:
        return os.environ["Pfam_data_dir"]
    except KeyError:
        wprint(color_text("The environment variable 'Pfam_data_dir' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `gtt data locations check`.")
        sys.exit(1)


def _load_pfam_info(info_path):
    """
    Build a lookup from the local pfamA.txt info table.

    Returns a dict keyed by version-less accession (e.g., 'PF00001') mapping to
    the profile NAME (e.g., '7tm_1'). pfamA.txt is tab-delimited with accession
    in column 1 and name in column 2.
    """
    acc_to_name = {}
    with open(info_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            acc, name = parts[0].strip(), parts[1].strip()
            if acc:
                acc_to_name[acc] = name
    return acc_to_name


def get_additional_pfam_targets(run_data):
    """
    Extract the requested Pfam profiles from the pre-downloaded master
    Pfam-A.hmm. Matching is done on the version-less
    accession; the full accession (with version) as it appears in the master
    HMM is recorded as the 'pulled' id.
    """
    pfam_data_dir = _check_pfam_data_dir()
    master_hmm_path = os.path.join(pfam_data_dir, HMM_FILENAME)
    info_path = os.path.join(pfam_data_dir, INFO_FILENAME)

    acc_to_name = _load_pfam_info(info_path)

    # read requested targets, grouped by version-less accession: the same Pfam can be
    # listed under more than one spelling (PF00001 and PF00001.3), and each spelling
    # used to get its own entry while only the last was ever matched, so the others
    # were reported as failed even though the profile was pulled
    pfam_dict = {}            # key: input_id (as given), value: pulled full acc or None
    inputs_by_core = {}       # key: version-less acc, value: [input_ids, in listed order]
    for input_id in read_single_column_file(run_data.target_pfams_file):
        if input_id in pfam_dict:
            continue
        pfam_dict[input_id] = None
        inputs_by_core.setdefault(input_id.split(".")[0], []).append(input_id)

    _report_duplicate_pfam_targets(inputs_by_core)

    # for profiles with no accession set, matched by NAME via the info table; built
    # once rather than scanning every wanted target for each of ~20k profiles
    core_by_name = {}
    for core in inputs_by_core:
        name = acc_to_name.get(core)
        if name:
            core_by_name.setdefault(name, core)

    profiles_dir = f"{run_data.pfam_results_dir}/target-pfam-profiles"
    os.makedirs(profiles_dir, exist_ok=True)

    combined_pfam_hmm_path = f"{profiles_dir}/all-pfam-targets.hmm"
    run_data.all_pfam_targets_hmm_path = combined_pfam_hmm_path

    found_pfam_targets = []
    found_pfam_paths = []
    pulled = set()

    with spinner("Collecting Pfam targets...", "", clear_on_done=True):
        # single streaming pass over the master HMM, pulling only wanted profiles
        with pyhmmer.plan7.HMMFile(master_hmm_path) as hmm_file:
            for hmm in hmm_file:
                acc = decode_pyhmmer_text(hmm.accession) or ""
                core = acc.split(".")[0]

                # fall back to matching on NAME if accession isn't set on the profile
                if core not in inputs_by_core:
                    name = decode_pyhmmer_text(hmm.name) or ""
                    core = core_by_name.get(name, core)

                if core not in inputs_by_core:
                    continue

                if core in pulled:
                    continue  # already pulled (guards against dup accessions)

                out_path = f"{profiles_dir}/{core}.hmm"
                with open(out_path, "wb") as out_fh:
                    hmm.write(out_fh)

                for input_id in inputs_by_core[core]:
                    pfam_dict[input_id] = acc if acc else core
                pulled.add(core)
                found_pfam_targets.append(core)
                found_pfam_paths.append(out_path)

        failed_pfam_targets = [iid for iid, pulled in pfam_dict.items() if pulled is None]

        run_data.found_pfam_targets = found_pfam_targets
        run_data.failed_pfam_targets = failed_pfam_targets
        run_data.pfam_dict = pfam_dict

        write_out_failed_pfams(run_data)
        write_requested_and_pulled_pfams(run_data, pfam_dict)

        if found_pfam_paths:
            concat_files(found_pfam_paths, combined_pfam_hmm_path)

    return run_data


def write_out_failed_pfams(run_data):
    if len(run_data.failed_pfam_targets) > 0:
        with open(run_data.pfam_results_dir + "/failed-pfam-targets.txt", "w") as failed_pfams_file:
            for entry in run_data.failed_pfam_targets:
                failed_pfams_file.write(entry + "\n")


def _report_duplicate_pfam_targets(inputs_by_core):
    """
    Say once, for each Pfam listed more than once under different spellings, that
    it's only searched once
    """
    for core, input_ids in inputs_by_core.items():
        if len(input_ids) > 1:
            wprint(color_text(
                f"    Note: {', '.join(input_ids)} are all the same Pfam ({core}), so "
                "it will be searched for just once.", "yellow"))


def write_requested_and_pulled_pfams(run_data, pfam_dict):
    with open(run_data.pfam_results_dir + "/info/requested-and-pulled-pfams.tsv", "w") as pfam_file:
        pfam_file.write("specified_pfam\tpulled_pfam\n")
        for input_id, pulled_id in pfam_dict.items():
            if pulled_id is None:
                pulled_id = "NA"
            pfam_file.write(f"{input_id}\t{pulled_id}\n")
