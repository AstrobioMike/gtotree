import os
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pyarrow.compute as pc # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, accession_core
from gtotree.utils.misc.general import newick_safe_label_part
from gtotree.utils.ncbi.get_ncbi_assembly_data import ncbi_data_table_path

_WANTED_COLUMNS = ["assembly_accession", "infraspecific_name"] + list(RANKS)

# extracts the numeric core (the run of digits) from a full assembly accession like
# "GCA_000739065.1" -> "000739065"; used Arrow-side to build the join key column
    # i'm 99% sure i'm not worried about any possible diff between a GCA's tax and it's GCF's tax
    # but if that comes up, here's where you want to fix it, Mike...
_ACCESSION_CORE_REGEX = r'^[A-Za-z]+_([0-9]+)\.[0-9]+$'


def update_mapping_dict_with_ncbi_tax_info(args, run_data):

    tax_rows = subset_ncbi_info(run_data)

    run_data.tax_info_dict = tax_rows

    for input_acc, ranks in tax_rows.items():
        if input_acc in run_data.mapping_dict:
            continue

        parts = [newick_safe_label_part(ranks.get(rank.lower(), "NA"))
                 for rank in args.lineage.split(",")]
        run_data.mapping_dict[input_acc] = "_".join([input_acc] + parts)

    return run_data


def subset_ncbi_info(run_data):
    """
    Pull NCBI lineage info for the remaining NCBI input accessions straight from the
    hosted Parquet asset

    Returns a dict: input_acc -> {rank: value for the 7 ranks, plus 'strain'}, with the
    species value de-duplicated of its leading genus and the strain de-duplicated of its
    leading species (mirroring the old taxonkit reformat, and the GTDB annotation path).

    The join is done on the version-stripped numeric accession core, which is shared
    between an assembly's GCA_ and GCF_ forms -- so an input GCF_ accession matches the
    GCA_ row (and vice-versa) even though NCBI lists both.
    """
    input_accs = run_data.remaining_ncbi_accs()

    # map each wanted accession to its numeric core for the join
    # a core can be shared by more than one input (the GCA_ and GCF_ of one assembly);
    # each of them gets the lineage, rather than only the last one seen
    core_to_inputs = {}
    for acc in input_accs:
        core = accession_core(acc)
        if core:
            core_to_inputs.setdefault(core, []).append(acc)

    if not core_to_inputs:
        return {}

    ncbi_path = ncbi_data_table_path(os.environ['NCBI_ASSEMBLY_DATA_DIR'])

    # Arrow-native filter: read the join column, extract each row's numeric core, keep
    # only rows whose core is one we want, then materialize just those matched rows.
    acc_col = pq.read_table(ncbi_path, columns=["assembly_accession"]).column(0)
    cores = pc.replace_substring_regex(acc_col, _ACCESSION_CORE_REGEX, r'\1')
    mask = pc.is_in(cores, value_set=pa.array(list(core_to_inputs.keys())))
    keep_idx = pc.indices_nonzero(mask)

    subset = pq.read_table(ncbi_path, columns=_WANTED_COLUMNS).take(keep_idx)

    acc_matched = subset.column("assembly_accession").to_pylist()
    infra_col = subset.column("infraspecific_name").to_pylist()
    rank_cols = {rank: subset.column(rank).to_pylist() for rank in RANKS}

    result = {}
    for i, acc in enumerate(acc_matched):
        for input_acc in core_to_inputs.get(accession_core(acc), ()):
            # a GCA and GCF row can share a core; keep the first match, but prefer a
            # row that actually has lineage (non-NA domain) over an NA one
            if input_acc in result and result[input_acc].get("domain", "NA") != "NA":
                continue
            result[input_acc] = _ranks_for_row(i, rank_cols, infra_col)

    return result


def _ranks_for_row(i, rank_cols, infra_col):
    """
    One matched NCBI row as {rank: value, ..., 'strain': value} for labelling
    """
    ranks = {rank: rank_cols[rank][i] for rank in RANKS}

    # capture the FULL "Genus species" before genus-stripping, for strain de-dup
    full_species = ranks.get("species")

    # species stores "Genus species"; strip the leading genus so labels don't
    # carry the genus twice
    genus = ranks.get("genus")
    species = ranks.get("species")
    if isinstance(species, str) and isinstance(genus, str):
        ranks["species"] = species.removeprefix(f"{genus} ")

    ranks["strain"] = _strain_from_infraspecific(infra_col[i], full_species)

    return ranks


def _strain_from_infraspecific(infraspecific_name, full_species):
    """
    Turn an NCBI infraspecific_name ("strain=<value>", "isolate=<value>", or "na") into
    a strain label, stripping a leading species prefix and collapsing spaces to hyphens.
    Returns "NA" when there's no usable strain.
    """
    if not isinstance(infraspecific_name, str):
        return "NA"
    value = infraspecific_name.strip()
    if not value or value.lower() == "na":
        return "NA"

    # infraspecific_name is "key=value" (strain=); taking the value
    if "=" in value:
        value = value.split("=", 1)[1].strip()
    if not value:
        return "NA"

    if isinstance(full_species, str) and full_species not in ("", "NA"):
        value = value.removeprefix(f"{full_species} ")

    return value.replace(" ", "-")
