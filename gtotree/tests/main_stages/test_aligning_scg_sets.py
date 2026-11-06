"""
Alignment scheduling vs output order: SCG-sets start aligning biggest-first (so a big
set doesn't run alone at the end), while the concatenated alignment always follows
the HMM-file order of the targets, whatever order the alignments finished in.
"""

import os
from datetime import datetime

from gtotree.main_stages import aligning_and_preparing_SCG_sets as stage
from gtotree.main_stages.aligning_and_preparing_SCG_sets import largest_first
from gtotree.main_stages.concatenating_SCG_sets import concatenate_SCG_sets
from gtotree.utils.misc.general import GenomeData, RunData, SCGset

GENOME_IDS = ("G0", "G1")


def _run_data(tmp_path, scg_ids):
    scg_dir = tmp_path / "out" / "run-files" / "SCG-seqs"
    scg_dir.mkdir(parents=True)
    rd = RunData()
    rd.start_time = datetime.now()
    rd.output_dir = str(tmp_path / "out")
    rd.run_files_dir = str(tmp_path / "out" / "run-files")
    rd.found_SCG_seqs_dir = str(scg_dir)
    rd.run_data_path = rd.run_files_dir + "/run-data.json"
    rd.general_ext = ".faa"
    rd.ncbi_accs = [GenomeData.from_acc(g) for g in GENOME_IDS]
    rd.update_all_input_genomes()
    rd.SCG_targets = [SCGset.from_id(s) for s in scg_ids]
    return rd


def _write_input(rd, scg_id, seq_len):
    path = f"{rd.found_SCG_seqs_dir}/{scg_id}-genome-filtered.faa"
    with open(path, "w") as fh:
        for g in GENOME_IDS:
            fh.write(f">{g}\n{'M' * seq_len}\n")


def test_largest_inputs_are_scheduled_first_and_ties_keep_hmm_order(tmp_path):
    rd = _run_data(tmp_path, ["small", "big", "mid_a", "mid_b", "missing"])
    for scg_id, n in [("small", 10), ("big", 500), ("mid_a", 100), ("mid_b", 100)]:
        _write_input(rd, scg_id, n)
    order = [s.id for s in largest_first(rd.SCG_targets, rd)]
    assert order == ["big", "mid_a", "mid_b", "small", "missing"]


def test_the_stage_hands_the_pool_the_largest_first(tmp_path, monkeypatch):
    rd = _run_data(tmp_path, ["small", "big"])
    _write_input(rd, "small", 10)
    _write_input(rd, "big", 500)
    seen = []

    def fake_pool(items, worker, apply_result, args, run_data, **kw):
        seen.extend(s.id for s in items)
        return run_data

    monkeypatch.setattr(stage, "run_pooled_stage", fake_pool)
    monkeypatch.setattr(stage, "write_SCG_info_table", lambda rd: None)
    monkeypatch.setattr(stage, "report_SCG_alignment_update", lambda rd: None)
    monkeypatch.setattr(stage, "report_processing_stage", lambda *a: None)
    args = type("A", (), {"keep_gene_alignments": False, "num_jobs": 1})()
    stage.align_and_prepare_SCG_sets(args, rd)
    assert seen == ["big", "small"]


def test_concatenation_follows_hmm_order_not_finishing_order(tmp_path):
    rd = _run_data(tmp_path, ["A", "B", "C"])
    # finished (and were marked ready) in the reverse of HMM order
    for scg in reversed(rd.SCG_targets):
        with open(f"{rd.found_SCG_seqs_dir}/{scg.id}-final.faa", "w") as fh:
            fh.write("".join(f">{g}\n{scg.id * 3}\n" for g in GENOME_IDS))
        scg.aligned = scg.trimmed = scg.ready_for_cat = True

    rd = concatenate_SCG_sets(rd)

    rows = open(rd.concatenated_alignment_path).read().split()
    assert rows == [">G0", "AAAXXXXXBBBXXXXXCCC", ">G1", "AAAXXXXXBBBXXXXXCCC"]
    parts = open(os.path.join(rd.run_files_dir, "partitions.txt")).read()
    assert [ln.split()[1] for ln in parts.splitlines()] == ["A", "B", "C"]
