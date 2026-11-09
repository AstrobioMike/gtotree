"""
Unit tests for gtotree/main_stages/processing_genomes.py.

The stage downloads/preprocesses and searches each genome in one worker call, then drops
that genome's sequence files. The tests below cover the observable results of that:
which genomes get processed, and what is left on disk afterwards.
"""

import datetime
import os
import shutil

import pytest  # type: ignore

from gtotree.tests.paths import DATA_DIR
from gtotree.utils.misc.general import GenomeData
from gtotree.main_stages.processing_genomes import (SearchPlan,
                                                    genomes_needing_processing)
from gtotree.utils.misc.stages import GenomeRemovalStage

# pyrodigal is what the nucleotide tests themselves need (gene calling is in-process).
# The rest are what `preflight_checks()` insists on before any stage runs -- the stage
# under test never calls them, but _process() goes through preflight, so a machine
# without them gets failures that look like real breakage instead of skips.
PREFLIGHT_REQUIRED = ["muscle", "trimal", "VeryFastTree"]

try:
    import pyrodigal  # type: ignore  # noqa: F401
    _missing = []
except ImportError:
    _missing = ["pyrodigal"]
_missing_preflight = [b for b in PREFLIGHT_REQUIRED if shutil.which(b) is None]

needs_preflight_deps = pytest.mark.skipif(
    bool(_missing_preflight),
    reason=f"missing binaries preflight_checks requires: {_missing_preflight}")


def _plan(do_pfam=False, do_ko=False, keep=False):
    return SearchPlan(do_pfam=do_pfam, do_ko=do_ko, keep_genome_files=keep)


def _genome(gid, **flags):
    gd = GenomeData.from_acc(gid)
    for k, v in flags.items():
        setattr(gd, k, v)
    return gd


class TestGenomesNeedingProcessing:
    """
    A genome's sequence files are deleted once its worker returns, so only fully
    processed genomes can be skipped, a processed-but-unsearched genome has nothing
    left on disk to search and is redone from the top.
    """

    def test_selects_unprocessed_and_partially_processed_genomes(self):
        done = _genome("G1", processing_done=True, hmm_search_done=True)
        partial = _genome("G2", processing_done=True, hmm_search_done=False)
        fresh = _genome("G3")
        removed = _genome("G4")
        removed.mark_removed("HMM search failed", GenomeRemovalStage.HMM_SEARCH)

        assert genomes_needing_processing([done, partial, fresh, removed], _plan()) == \
            [partial, fresh]

    def test_extra_searches_are_required_only_when_planned(self):
        gd = _genome("G1", processing_done=True, hmm_search_done=True,
                     pfam_search_done=False)
        assert genomes_needing_processing([gd], _plan(do_pfam=True)) == [gd]
        assert genomes_needing_processing([gd], _plan(do_pfam=False)) == []


# ---------------------------------------------------------------------------
# end-to-end behaviour of the stage
# ---------------------------------------------------------------------------

def _process(tmp_path, source_flag, listing_paths, hmm, extra=()):
    from gtotree.cli.parser import parser
    from gtotree.utils.misc.preflight_checks import preflight_checks
    from gtotree.main_stages.processing_genomes import process_genomes
    from gtotree.utils.misc.context import log_file_var

    listing = tmp_path / "inputs.txt"
    listing.write_text("\n".join(str(p) for p in listing_paths) + "\n")

    argv = [source_flag, str(listing), "-H", str(hmm),
            "-o", str(tmp_path / "out"), "-j", "2", *extra]

    token = log_file_var.set(log_file_var.get())
    try:
        args, run_data = preflight_checks(parser().parse_args(argv))
        run_data.start_time = datetime.datetime.now()
        return args, process_genomes(args, run_data)
    finally:
        log_file_var.reset(token)


AA_GENOMES = [DATA_DIR / f"mock-{i}.faa" for i in range(1, 5)]
MOCK_HMM = DATA_DIR / "mock.hmm"


@needs_preflight_deps
def test_processing_searches_every_genome_and_drops_its_sequence_files(tmp_path):
    _, run_data = _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM)

    assert len(run_data.all_input_genomes) == 4
    assert all(gd.hmm_search_done for gd in run_data.all_input_genomes)
    assert os.listdir(run_data.ready_genome_files_dir) == []


@needs_preflight_deps
def test_keep_working_dir_keeps_the_sequence_files(tmp_path):
    _, run_data = _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM, extra=["--keep-working-dir"])
    assert sorted(os.listdir(run_data.ready_genome_files_dir)) == \
        ["mock-1.faa", "mock-2.faa", "mock-3.faa", "mock-4.faa"]


@needs_preflight_deps
def test_running_the_stage_again_skips_everything_and_keeps_outputs_complete(tmp_path):
    from gtotree.main_stages.processing_genomes import process_genomes

    args, run_data = _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM)
    table = os.path.join(run_data.output_dir, "SCG-hit-counts.tsv")
    first = open(table).read()

    process_genomes(args, run_data)

    assert open(table).read() == first
    assert len([ln for ln in first.splitlines() if ln.strip()]) == 5  # header + 4


@needs_preflight_deps
def test_a_resume_after_losing_per_genome_results_reprocesses_those_genomes(tmp_path):
    # with the per-genome files gone (scratch --tmp-dir wiped), the rebuild skipped
    # the genome but its flags still said done, so its hits silently vanished
    import shutil
    from gtotree.main_stages.processing_genomes import process_genomes
    from gtotree.utils.misc.stages import PipelineStage

    args, run_data = _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM)
    table = os.path.join(run_data.output_dir, "SCG-hit-counts.tsv")
    first = open(table).read()
    counts = {s.id: s.num_genomes_with_hit for s in run_data.SCG_targets}

    shutil.rmtree(os.path.join(run_data.hmm_results_dir, "mock-2"))
    run_data = process_genomes(args, run_data)

    assert open(table).read() == first
    assert {s.id: s.num_genomes_with_hit for s in run_data.SCG_targets} == counts
    assert run_data.stage_is_complete(PipelineStage.PROCESS_GENOMES)


@needs_preflight_deps
def test_a_resume_that_regenerates_results_brings_back_genomes_removed_later(tmp_path):
    # the `-G` removal was judged against the run as it stood; with per-genome results
    # to regenerate, every later stage reruns, so the genome has to be back in the run
    # (and have its own missing results regenerated) rather than staying dropped
    import shutil
    from gtotree.main_stages.processing_genomes import process_genomes
    from gtotree.utils.misc.stages import PipelineStage

    args, run_data = _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM)
    table = os.path.join(run_data.output_dir, "SCG-hit-counts.tsv")
    first = open(table).read()
    for stage in (PipelineStage.FILTER_GENES, PipelineStage.FILTER_GENOMES):
        run_data.mark_stage_complete(stage)

    dropped = next(gd for gd in run_data.all_input_genomes if gd.id == "mock-2")
    dropped.mark_removed("too few unique SCG hits", GenomeRemovalStage.SCG_HIT_FILTER)
    shutil.rmtree(os.path.join(run_data.hmm_results_dir, "mock-2"))

    run_data = process_genomes(args, run_data)

    assert not dropped.removed
    assert dropped.hmm_search_done
    assert open(table).read() == first
    assert not run_data.stage_is_complete(PipelineStage.FILTER_GENOMES)


@needs_preflight_deps
def test_a_resume_does_not_rebuild_counts_without_genomes_removed_later(tmp_path):
    # re-running the combine on a resume dropped genomes that `-G` removed AFTER this
    # stage, so the SCG counts and SCG-hit-counts.tsv drifted from an uninterrupted run
    from gtotree.main_stages.processing_genomes import process_genomes
    from gtotree.utils.misc.stages import GenomeRemovalStage

    args, run_data = _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM)
    table = os.path.join(run_data.output_dir, "SCG-hit-counts.tsv")
    first = open(table).read()
    counts = {s.id: s.num_genomes_with_hit for s in run_data.SCG_targets}

    run_data.all_input_genomes[0].mark_removed("too few unique SCG hits",
                                               GenomeRemovalStage.SCG_HIT_FILTER)
    run_data = process_genomes(args, run_data)

    assert open(table).read() == first
    assert {s.id: s.num_genomes_with_hit for s in run_data.SCG_targets} == counts


@needs_preflight_deps
def test_a_stop_partway_through_a_pool_resumes_from_the_genomes_already_done(
        tmp_path, monkeypatch):
    # run-data used to be saved only once a whole source's pool finished, so a run
    # stopped partway through (e.g. hours into 20k genomes) redid every genome
    import gtotree.main_stages.processing_genomes as pg
    from gtotree.utils.misc.general import read_run_data

    # the search results are the last thing applied for a genome, so stopping right
    # after the second one leaves exactly two genomes fully done
    real_apply = pg._apply_searches
    applied = []

    def stop_after_two(gd, searches, run_data, plan):
        real_apply(gd, searches, run_data, plan)
        applied.append(gd.id)
        if len(applied) == 2:
            raise KeyboardInterrupt

    monkeypatch.setattr(pg, "_apply_searches", stop_after_two)
    with pytest.raises(KeyboardInterrupt):
        _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM, extra=["-j", "1"])

    saved = read_run_data(str(tmp_path / "out" / "run-files" / "run-data.json"))
    done = sorted(gd.id for gd in saved.amino_acid_files if gd.hmm_search_done)
    assert done == sorted(applied)

    monkeypatch.setattr(pg, "_apply_searches", real_apply)
    redone = []
    real_worker = pg._process_one_amino_acid_file

    def counting_worker(gd, run_data):
        redone.append(gd.id)
        return real_worker(gd, run_data)

    monkeypatch.setattr(pg, "_process_one_amino_acid_file", counting_worker)
    _, run_data = _process(tmp_path, "-A", AA_GENOMES, MOCK_HMM, extra=["-R"])
    assert sorted(redone) == sorted(set(f"mock-{i}" for i in range(1, 5)) - set(applied))
    assert all(gd.hmm_search_done for gd in run_data.all_input_genomes)


@needs_preflight_deps
def test_a_resume_with_nothing_saved_starts_fresh_and_honors_tmp_dir(tmp_path):
    from gtotree.cli.parser import parser
    from gtotree.utils.misc.preflight_checks import preflight_checks
    from gtotree.utils.misc.context import log_file_var

    listing = tmp_path / "inputs.txt"
    listing.write_text("\n".join(str(p) for p in AA_GENOMES) + "\n")
    scratch = tmp_path / "scratch"
    argv = ["-A", str(listing), "-H", str(MOCK_HMM), "-o", str(tmp_path / "out"),
            "-R", "--tmp-dir", str(scratch)]
    token = log_file_var.set(log_file_var.get())
    try:
        args, run_data = preflight_checks(parser().parse_args(argv))
    finally:
        log_file_var.reset(token)
    assert os.path.dirname(run_data.tmp_dir) == str(scratch)
    assert args.tmp_dir == run_data.tmp_dir


@needs_preflight_deps
def test_gzipped_genomes_are_read_in_place_whatever_they_are_named(tmp_path):
    # gzipped inputs used to be decompressed to a scratch copy first, chosen by a
    # `.gz` suffix, so a gzipped file without one was read as text and dropped
    import gzip
    inputs = []
    for i, name in enumerate(["mock-1.faa.gz", "mock-2.faa"], start=1):
        dest = tmp_path / name
        with open(DATA_DIR / f"mock-{i}.faa", "rb") as src, gzip.open(dest, "wb") as out:
            out.write(src.read())
        inputs.append(dest)
    inputs += AA_GENOMES[2:]

    _, run_data = _process(tmp_path, "-A", inputs, MOCK_HMM, extra=["--keep-working-dir"])

    assert not [gd.id for gd in run_data.all_input_genomes if gd.removed]
    assert all(gd.hmm_search_done for gd in run_data.all_input_genomes)
    gz = {gd.id: gd.was_gzipped for gd in run_data.all_input_genomes}
    assert gz == {"mock-1": True, "mock-2": True, "mock-3": False, "mock-4": False}
    assert os.listdir(run_data.AA_processing_dir) == []


# ---------------------------------------------------------------------------
# nucleotide mode
# ---------------------------------------------------------------------------

@needs_preflight_deps
@pytest.mark.skipif(bool(_missing), reason=f"missing binaries: {_missing}")
def test_gene_calling_intermediates_are_not_left_behind(tmp_path, nt_genomes, nt_hmm):
    # prodigal's per-genome protein/CDS files (and raw NCBI downloads) were never
    # cleaned up, building up ~5 MB a genome in the tmp dir for the whole run
    _, run_data = _process(tmp_path, "-f", nt_genomes, nt_hmm)
    assert os.listdir(run_data.fasta_processing_dir) == []


@needs_preflight_deps
@pytest.mark.skipif(bool(_missing), reason=f"missing binaries: {_missing}")
def test_keep_working_dir_keeps_the_gene_calling_intermediates(tmp_path, nt_genomes,
                                                               nt_hmm):
    _, run_data = _process(tmp_path, "-f", nt_genomes, nt_hmm,
                           extra=["--keep-working-dir"])
    assert os.listdir(run_data.fasta_processing_dir)


@needs_preflight_deps
@pytest.mark.skipif(bool(_missing), reason=f"missing binaries: {_missing}")
def test_nucleotide_input_has_genes_called_and_targets_found(tmp_path, nt_genomes, nt_hmm):
    _, run_data = _process(tmp_path, "-f", nt_genomes, nt_hmm)

    assert all(gd.prodigal_used for gd in run_data.all_input_genomes)
    assert all(gd.num_genes > 50 for gd in run_data.all_input_genomes)
    assert all(gd.num_unique_SCG_hits == 1 for gd in run_data.all_input_genomes)
    # prodigal output must have its stop characters stripped
    for name in os.listdir(run_data.hmm_results_dir):
        hits = os.path.join(run_data.hmm_results_dir, name, "SCG-hits.faa")
        if os.path.isfile(hits):
            assert "*" not in open(hits).read()


@needs_preflight_deps
@pytest.mark.skipif(bool(_missing), reason=f"missing binaries: {_missing}")
def test_nucleotide_mode_yields_nucleotide_hits_and_drops_both_files(
        tmp_path, nt_genomes, nt_hmm):
    """In `-z` a genome has both a .faa and a .fasta, and both are dropped."""
    _, run_data = _process(tmp_path, "-f", nt_genomes, nt_hmm, extra=["-z"])

    assert run_data.general_ext == ".fasta"
    hits = os.path.join(run_data.found_SCG_seqs_dir, "NT_MOCK.fasta")
    seqs, cur = {}, None
    for line in open(hits):
        if line.startswith(">"):
            cur = line[1:].strip(); seqs[cur] = ""
        else:
            seqs[cur] += line.strip()

    assert len(seqs) == 4
    for name, seq in seqs.items():
        assert set(seq) <= set("ACGTN"), f"{name} is not nucleotide"
        assert len(seq) % 3 == 0

    assert os.listdir(run_data.ready_genome_files_dir) == []


################################################################################
# do_scg gating
################################################################################

"""
`SearchPlan.do_scg` exists so `gtt search-annotations` can reuse this
module's fused preprocess-then-search worker without an SCG set or a tree. It defaults
to True, so the main pipeline never opts in -- which is exactly why it needs testing
from this side too: a regression that ignored the flag would leave the subcommands
trying to press and search `run_data.hmm_path`, which they never set.
"""


def test_search_plan_defaults_to_doing_the_scg_search():
    """The main driver must never have to opt in."""
    assert _plan().do_scg is True


def test_do_scg_false_is_carried_on_the_plan():
    plan = SearchPlan(do_pfam=True, do_ko=False, keep_genome_files=False, do_scg=False)
    assert plan.do_scg is False
    assert plan.do_pfam is True


def test_genome_is_fully_processed_ignores_the_hmm_search_when_do_scg_is_off():
    """
    The gate that decides whether a genome still needs work. With do_scg off, a
    processed-and-pfam-searched genome is done; if the flag were ignored here, every
    genome would be re-processed on every run because hmm_search_done never becomes
    True.
    """
    from gtotree.main_stages.processing_genomes import genome_is_fully_processed

    gd = _genome("g1", processing_done=True, pfam_search_done=True)

    scg_plan = SearchPlan(do_pfam=True, do_ko=False, keep_genome_files=False)
    no_scg_plan = SearchPlan(do_pfam=True, do_ko=False, keep_genome_files=False,
                             do_scg=False)

    assert genome_is_fully_processed(gd, scg_plan) is False
    assert genome_is_fully_processed(gd, no_scg_plan) is True


def test_genomes_needing_processing_respects_do_scg():
    no_scg_plan = SearchPlan(do_pfam=True, do_ko=False, keep_genome_files=False,
                             do_scg=False)

    done = _genome("done", processing_done=True, pfam_search_done=True)
    needs_search = _genome("needs", processing_done=True, pfam_search_done=False)

    remaining = genomes_needing_processing([done, needs_search], no_scg_plan)
    assert [gd.id for gd in remaining] == ["needs"]


def test_run_searches_skips_the_hmm_worker_when_do_scg_is_off(monkeypatch):
    """
    The load-bearing assertion: with do_scg off the SCG worker must not be called at
    all. Calling it would fail on the unset hmm_path rather than no-op, so this is the
    difference between the subcommands working and crashing.
    """
    from gtotree.main_stages import processing_genomes as pg

    called = []

    monkeypatch.setattr(pg, "_hmm_search_worker",
                        lambda *a, **kw: called.append("hmm") or {})
    monkeypatch.setattr(pg, "_pfam_search_worker",
                        lambda *a, **kw: called.append("pfam") or {})

    plan = SearchPlan(do_pfam=True, do_ko=False, keep_genome_files=False, do_scg=False)
    status = {"final_AA_path": "/tmp/x.faa", "final_nt_path": None}

    results = pg._run_searches(_genome("g1"), object(), plan, status)

    assert called == ["pfam"]
    assert "hmm" not in results


def test_run_searches_still_runs_the_hmm_worker_by_default(monkeypatch):
    """The other direction, so the gate can't be left permanently off."""
    from gtotree.main_stages import processing_genomes as pg

    called = []
    monkeypatch.setattr(pg, "_hmm_search_worker",
                        lambda *a, **kw: called.append("hmm") or {})

    plan = SearchPlan(do_pfam=False, do_ko=False, keep_genome_files=False)
    status = {"final_AA_path": "/tmp/x.faa", "final_nt_path": None}

    pg._run_searches(_genome("g1"), object(), plan, status)
    assert called == ["hmm"]


def test_apply_searches_skips_the_hmm_result_when_do_scg_is_off(monkeypatch):
    from gtotree.main_stages import processing_genomes as pg

    called = []
    monkeypatch.setattr(pg, "_apply_hmm_search_result",
                        lambda *a, **kw: called.append("hmm"))

    plan = SearchPlan(do_pfam=False, do_ko=False, keep_genome_files=False, do_scg=False)
    # a stale "hmm" key must be ignored rather than applied
    pg._apply_searches(_genome("g1"), {"hmm": {}}, object(), plan)

    assert called == []
