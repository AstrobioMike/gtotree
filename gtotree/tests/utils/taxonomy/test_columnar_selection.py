"""
The columnar pieces of selection: the vectorized accession core, the cached NCBI
liveness set, and dereplication streaming its pool in batches.

These replaced per-row Python over multi-million-row pools, so the property that
matters is that they give exactly the answers the scalar versions did.
"""

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore

from gtotree.utils.taxonomy import tax_derep, tax_select
from gtotree.utils.taxonomy.tax_ranks import (RANKS, accession_core,
                                              accession_core_array)
from gtotree.utils.taxonomy.exclusion_list import filter_table_by_exclusion


TRICKY = ["RS_GCF_000005845.2", "GCA_000005845.3", " gcf_123.1 ", "na", "", None,
          "GB_GCA_ABC_99.4", "GCF_1", "RS_X", "GCF_000001.1.2", "foo_GCA_1",
          "GCA_", "GCA_.1", "GCF_000000001.1"]


class TestAccessionCoreArray:

    def test_matches_the_scalar_version_on_awkward_input(self):
        got = accession_core_array(pa.array(TRICKY)).to_pylist()
        assert got == [accession_core(a if a is not None else "") for a in TRICKY]

    def test_the_canonical_fast_path_matches_the_scalar_version(self):
        canonical = ["GCF_000005845.2", "GCA_000001405.29", "GCA_900000005.1"]
        got = accession_core_array(pa.array(canonical)).to_pylist()
        assert got == [accession_core(a) for a in canonical]

    def test_chunked_input_is_handled(self):
        arr = pa.chunked_array([["GCF_000005845.2"], ["RS_GCF_000000001.1"]])
        assert accession_core_array(arr).to_pylist() == ["000005845", "000000001"]


def _ncbi_table(path, accs):
    pq.write_table(pa.table({"assembly_accession": pa.array(accs, type=pa.string())}),
                   str(path))


class TestLiveness:

    def test_live_cores_match_the_scalar_set(self, tmp_path):
        path = tmp_path / "ncbi.parquet"
        accs = ["GCF_000000001.1", "GCA_000000002.3", "", None, "GCF_000000001.2"]
        _ncbi_table(path, accs)
        got = set(tax_select.live_accession_cores(str(path)).to_pylist())
        assert got == {accession_core(a) for a in accs if a}

    def test_a_rewritten_table_is_not_served_from_the_cache(self, tmp_path):
        path = tmp_path / "ncbi.parquet"
        _ncbi_table(path, ["GCF_000000001.1"])
        assert tax_select.live_accession_cores(str(path)).to_pylist() == ["000000001"]
        _ncbi_table(path, ["GCF_000000001.1", "GCF_000000002.1"])
        assert sorted(tax_select.live_accession_cores(str(path)).to_pylist()) == \
            ["000000001", "000000002"]

    def test_filter_drops_dead_and_empty_accessions(self, tmp_path):
        path = tmp_path / "ncbi.parquet"
        _ncbi_table(path, ["GCF_000000001.1", "GCF_000000003.1"])
        # GTDB-style prefixes and a GCA/GCF pairing still match by core
        tab = pa.table({"acc": ["RS_GCF_000000001.1", "GB_GCA_000000002.1",
                                "GB_GCA_000000003.1", ""]})
        kept, dropped = tax_select.filter_table_by_liveness(tab, "acc", str(path))
        assert kept.column("acc").to_pylist() == ["RS_GCF_000000001.1",
                                                  "GB_GCA_000000003.1"]
        assert dropped == 2


class TestTableExclusion:

    def test_matches_the_row_version(self):
        tab = pa.table({"acc": TRICKY})
        exclude = {"000005845", "99"}
        kept, n = filter_table_by_exclusion(tab, "acc", exclude)
        expected = [a for a in TRICKY
                    if accession_core(a or "") not in exclude]
        assert kept.column("acc").to_pylist() == expected
        assert n == len(TRICKY) - len(expected)


def _gtdb_rows(n_groups, per_group):
    rows = []
    for g in range(n_groups):
        for i in range(per_group):
            acc = f"GCA_{g * per_group + i:09d}.1"
            lineage = ("Bacteria", "Testophyla", f"Class{g}", f"Ord{g}", f"Fam{g}",
                       f"Gen{g}", f"Gen{g} sp{i}")
            # completeness varies so each group has a unique, non-first best
            comp = str(50 + ((i * 7 + g) % 50))
            row = {"ncbi_genbank_assembly_accession": acc, "gtdb_representative": "t",
                   "ncbi_refseq_category": "", "checkm2_completeness": comp,
                   "checkm2_contamination": "1.0", "genome_size": "4000000",
                   "contig_count": "10"}
            row.update(dict(zip(RANKS, lineage)))
            rows.append(row)
    return rows


def _scalar_best(rows, rank, spec):
    best = {}
    for row in rows:
        group = row[rank]
        k = tax_derep.quality_key(row, spec)
        if group not in best or k < tax_derep.quality_key(best[group], spec):
            best[group] = row
    return [best[g][spec.acc_col] for g in sorted(best)]


def _loop_quality(row, spec):
    return tax_derep.quality_key(row, spec)


def test_batched_derep_picks_what_the_all_at_once_loop_did(tmp_path, monkeypatch):
    """Batches far smaller than the pool, so groups straddle batch boundaries"""
    rows = _gtdb_rows(n_groups=7, per_group=13)
    path = tmp_path / "gtdb.parquet"
    pq.write_table(pa.Table.from_pylist(rows), str(path))
    monkeypatch.setattr(tax_derep, "DEREP_BATCH_ROWS", 5)

    # a callable picker takes the row-batch loop (the built-in one is columnar)
    accs, n_groups, _w, attrition = tax_derep.derep(
        str(path), "gtdb", "phylum", "Testophyla", "class", pick=_loop_quality)

    assert accs == _scalar_best(rows, "class", tax_select.SOURCES["gtdb"])
    assert n_groups == 7
    assert attrition.n_after_floor == len(rows)


def test_batched_derep_counts_the_quality_floor_across_batches(tmp_path, monkeypatch):
    rows = _gtdb_rows(n_groups=3, per_group=10)
    path = tmp_path / "gtdb.parquet"
    pq.write_table(pa.Table.from_pylist(rows), str(path))
    monkeypatch.setattr(tax_derep, "DEREP_BATCH_ROWS", 4)

    _accs, _n, _w, attrition = tax_derep.derep(
        str(path), "gtdb", "phylum", "Testophyla", "class", min_completeness=80,
        pick=_loop_quality)

    n_below = sum(1 for r in rows if float(r["checkm2_completeness"]) < 80)
    assert attrition.n_below_floor == n_below
    assert attrition.n_after_floor == len(rows) - n_below


################################################################################
# the columnar quality pick
################################################################################

def _awkward_ncbi_rows(n, seed):
    """NCBI-shaped rows hitting every branch of quality_key and the floor"""
    import random
    rng = random.Random(seed)
    pick = rng.choice
    rows = []
    for i in range(n):
        g = rng.randrange(12)
        row = {
            # duplicate accessions too, so whole-key ties fall to row order
            "assembly_accession": pick([f"GCF_{rng.randrange(40):09d}.1", "", None]),
            "refseq_category": pick(["reference genome", "na", "Reference Genome ",
                                     None, ""]),
            "checkm_completeness": pick(["NA", "na", "", None, "84.9", "85", "99.5",
                                         "100", "50", " 97.1 ", "1e2"]),
            "checkm_contamination": pick(["NA", None, "0", "10", "10.01", "2.5",
                                          "0.0"]),
            "assembly_level": pick(["Complete Genome", "Chromosome", "Scaffold",
                                    "Contig", "contig", None, "weird"]),
            "contig_count": pick(["0", "1", "250", "NA", None, "-3"]),
            "genome_size_ungapped": pick(["5000000", "", None, "0", "NA"]),
            "genome_size": pick(["5100000", None, "4000000"]),
        }
        lineage = ["Bacteria", "Testophyla", f"Class{g % 4}", f"Ord{g}",
                   pick([f"Fam{g}", "NA", "", None]), f"Gen{g}", f"Gen{g} sp"]
        row.update(dict(zip(RANKS, lineage)))
        rows.append(row)
    return rows


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("floor", [{}, {"min_completeness": 90},
                                   {"max_contamination": 5},
                                   {"min_completeness": 85, "max_contamination": 10}])
def test_the_columnar_pick_matches_the_python_loop(tmp_path, seed, floor):
    rows = _awkward_ncbi_rows(400, seed)
    path = tmp_path / "ncbi.parquet"
    pq.write_table(pa.Table.from_pylist(rows), str(path))

    fast = tax_derep.derep(str(path), "ncbi", "phylum", "Testophyla", "family",
                           **floor)
    slow = tax_derep.derep(str(path), "ncbi", "phylum", "Testophyla", "family",
                           pick=_loop_quality, **floor)

    assert fast[0] == slow[0]
    assert fast[1] == slow[1]
    assert fast[2] == slow[2]
    assert vars(fast[3]) == vars(slow[3])


def test_the_columnar_pick_prefers_a_sound_reference_genome(tmp_path):
    lineage = dict(zip(RANKS, ["Bacteria", "Testophyla", "C", "O", "F", "G", "G sp"]))

    def row(acc, ref, comp, cont):
        return {"assembly_accession": acc, "refseq_category": ref,
                "checkm_completeness": comp, "checkm_contamination": cont,
                "assembly_level": "Contig", "contig_count": "5",
                "genome_size_ungapped": "1000", "genome_size": "1000", **lineage}

    rows = [row("GCF_000000001.1", "na", "100", "0"),
            row("GCF_000000002.1", "reference genome", "90", "3"),
            row("GCF_000000003.1", "reference genome", "70", "1")]  # too incomplete
    path = tmp_path / "ncbi.parquet"
    pq.write_table(pa.Table.from_pylist(rows), str(path))
    accs, *_ = tax_derep.derep(str(path), "ncbi", "phylum", "Testophyla", "family")
    assert accs == ["GCF_000000002.1"]


def test_taxon_lookups_reuse_the_rank_names_until_the_table_changes(tmp_path):
    rows = [dict(zip(RANKS, ["Bacteria", "Testophyla", "C", "O", "F", "G", "G sp"]))]
    path = tmp_path / "t.parquet"
    pq.write_table(pa.Table.from_pylist(rows), str(path))
    assert tax_select.find_ranks_for_taxon(str(path), "testophyla")[:2] == \
        ("Testophyla", ["phylum"])

    rows[0]["phylum"] = "Renamedota"
    pq.write_table(pa.Table.from_pylist(rows), str(path))
    import os
    st = os.stat(path)
    os.utime(path, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
    assert tax_select.find_ranks_for_taxon(str(path), "renamedota")[1] == ["phylum"]
    with pytest.raises(tax_select.TaxonNotFound):
        tax_select.find_ranks_for_taxon(str(path), "testophyla")


def test_the_default_pick_never_turns_the_pool_into_python_rows(tmp_path, monkeypatch):
    # the whole point: a 4M-row NCBI pool ranked as columns, not dicts
    rows = _awkward_ncbi_rows(200, seed=1)
    path = tmp_path / "ncbi.parquet"
    pq.write_table(pa.Table.from_pylist(rows), str(path))

    def no_rows(*a, **k):
        raise AssertionError("the pool was materialized as Python rows")

    monkeypatch.setattr(tax_derep, "apply_quality_floor", no_rows)
    accs, *_ = tax_derep.derep(str(path), "ncbi", "phylum", "Testophyla", "family")
    assert accs


def test_a_second_taxon_lookup_does_not_reread_the_rank_columns(tmp_path, monkeypatch):
    rows = [dict(zip(RANKS, ["Bacteria", "Testophyla", "C", "O", "F", "G", "G sp"]))]
    path = tmp_path / "t.parquet"
    pq.write_table(pa.Table.from_pylist(rows), str(path))
    tax_select.find_ranks_for_taxon(str(path), "Testophyla")

    reads = []
    real_read = pq.read_table

    def counting_read(p, columns=None, **kw):
        if not kw.get("filters"):
            reads.append(columns)
        return real_read(p, columns=columns, **kw)

    monkeypatch.setattr(pq, "read_table", counting_read)
    assert tax_select.find_ranks_for_taxon(str(path), "G")[1] == ["genus"]
    assert reads == []
