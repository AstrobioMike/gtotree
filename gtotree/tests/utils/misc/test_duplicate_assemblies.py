"""
One assembly reaching a run under two accessions (its GCA_ and GCF_, or two versions)
is the same genome twice: a zero-length pair in the tree, and only one of the two
used to get a taxonomy label.
"""

import pytest  # type: ignore

from gtotree.utils.misc.general import GenomeData, RunData


def _run_data(*accs):
    rd = RunData()
    rd.ncbi_accs = [GenomeData.from_acc(a) for a in accs]
    rd.update_all_input_genomes()
    return rd


class TestWantedRefTaxMerge:

    def test_a_w_accession_for_an_assembly_already_given_is_not_added(self):
        # GTDB hands back GCA_; the user's -a list often has the GCF_ or another version
        rd = _run_data("GCF_000005845.2")
        added = rd.merge_wanted_ref_tax_accessions(
            ["GCA_000005845.2", "GCA_000005845.1", "GCA_000009065.1"])
        assert added == 1
        assert [gd.id for gd in rd.ncbi_accs] == ["GCF_000005845.2", "GCA_000009065.1"]

    def test_two_w_selections_sharing_an_assembly_add_it_once(self):
        rd = _run_data()
        rd.merge_wanted_ref_tax_accessions(["GCA_000009065.1"])
        assert rd.merge_wanted_ref_tax_accessions(["GCF_000009065.1"]) == 0


class TestDuplicateAssembliesInTheInput:

    def test_they_are_collected_by_accession_core(self):
        from gtotree.utils.misc.preflight_checks import collect_duplicate_assemblies
        rd = _run_data("GCA_000005845.2", "GCF_000005845.2", "GCF_000009065.1")
        dups = collect_duplicate_assemblies(rd)
        assert list(dups) == ["000005845"]
        assert [gd.id for gd in dups["000005845"]] == ["GCA_000005845.2",
                                                       "GCF_000005845.2"]

    def test_the_main_run_refuses_them(self, capsys):
        from gtotree.utils.misc.preflight_checks import check_for_duplicate_assemblies
        rd = _run_data("GCA_000005845.2", "GCF_000005845.2")
        with pytest.raises(SystemExit):
            check_for_duplicate_assemblies(rd)
        assert "GCA_000005845.2  /  GCF_000005845.2" in capsys.readouterr().out

    def test_the_subcommands_shared_resolver_refuses_them(self, tmp_path):
        import argparse
        from gtotree.utils.misc.general import resolve_input_genomes

        class _Err(Exception):
            pass

        args = argparse.Namespace(ncbi_accessions=None, wanted_ref_tax=None,
                                  genbank_files=None, fasta_files=None,
                                  amino_acid_files=None, source="gtdb",
                                  target_rank=None, derep_rank="off",
                                  min_completeness=None, max_contamination=None)
        rd = _run_data("GCA_000005845.2", "GCF_000005845.1")
        with pytest.raises(_Err, match="GCA_000005845.2  /  GCF_000005845.1"):
            resolve_input_genomes(args, rd, _Err)


def test_both_accessions_of_one_assembly_get_gtdb_taxonomy(tmp_path):
    # the join keyed each core to a single input, so the other got no lineage
    import pyarrow as pa  # type: ignore
    import pyarrow.parquet as pq  # type: ignore
    from gtotree.utils.taxonomy.tax_ranks import RANKS
    from gtotree.utils.gtdb import handle_gtdb_tax_info as mod

    row = {"ncbi_genbank_assembly_accession": "GCA_000005845.2"}
    row.update(dict(zip(RANKS, ["Bacteria", "Pseudomonadota", "Gamma", "Entero",
                                "Enterobacteriaceae", "Escherichia",
                                "Escherichia coli"])))
    cols = [c for c in mod._WANTED_COLUMNS]
    table = pa.table({c: pa.array([row.get(c, "NA")], type=pa.string()) for c in cols})
    path = tmp_path / "gtdb.parquet"
    pq.write_table(table, str(path))

    result = mod.gtdb_lineages_for_accessions(["GCF_000005845.2", "GCA_000005845.2"],
                                              gtdb_path=str(path))
    assert set(result) == {"GCF_000005845.2", "GCA_000005845.2"}
