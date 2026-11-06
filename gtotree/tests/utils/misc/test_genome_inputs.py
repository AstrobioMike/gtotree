"""
The four genome-input flags each take one or more values: genome files themselves,
quoted patterns GToTree expands, files listing them, and for `-a` accessions too.
"""

import gzip
import os
import shutil
import types

import pytest  # type: ignore

from gtotree.tests.paths import DATA_DIR
from gtotree.utils.misc import genome_inputs as gi
from gtotree.utils.misc.genome_inputs import (GenomeInputError, ResolvedInputs,
                                              classify_file, input_entries,
                                              resolve_accessions, resolve_genome_files,
                                              resolve_genome_input_args)


@pytest.fixture
def genomes(tmp_path, monkeypatch):
    """Four small fasta genomes in tmp_path/fastas, with cwd at tmp_path"""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "fastas").mkdir()
    paths = []
    for i in range(1, 5):
        p = tmp_path / "fastas" / f"g{i}.fasta"
        p.write_text(f">c{i}\nACGTACGT\n")
        paths.append(f"fastas/g{i}.fasta")
    return paths


################################################################################
# telling genome files from list files
################################################################################

class TestClassify:

    def test_a_genome_extension_is_a_genome_without_opening_it(self, tmp_path):
        p = tmp_path / "big.fna.gz"
        p.write_bytes(b"x" * (gi.TINY_FILE_BYTES + 10))   # not even real gzip
        assert classify_file(str(p)) == "genome"

    def test_an_unknown_extension_is_judged_by_its_contents(self, tmp_path):
        genome = tmp_path / "genome.seqs"
        genome.write_text(">c\nACGT\n")
        listing = tmp_path / "inputs.txt"
        listing.write_text(f"{genome}\n")
        assert classify_file(str(genome)) == "genome"
        assert classify_file(str(listing)) == "listing"

    def test_gzipped_contents_are_looked_through(self, tmp_path):
        gb = tmp_path / "genome.dat.gz"
        with gzip.open(gb, "wt") as fh:
            fh.write("LOCUS       X\n")
        assert classify_file(str(gb)) == "genome"

    def test_a_list_file_saved_under_a_genome_extension_is_still_a_list(self, genomes):
        listing = "list.fasta"
        with open(listing, "w") as fh:
            fh.write("\n".join(genomes) + "\n")
        assert classify_file(listing) == "listing"


################################################################################
# resolving file flags
################################################################################

class TestResolveGenomeFiles:

    def test_a_list_file_is_read_as_before(self, genomes):
        with open("list.txt", "w") as fh:
            fh.write("\n".join(genomes) + "\n")
        r = resolve_genome_files(["list.txt"], "-f")
        assert r.entries == genomes
        assert r.listing_files == ["list.txt"] and r.num_direct == 0

    def test_genome_files_can_be_given_directly(self, genomes):
        r = resolve_genome_files(genomes[:2], "-f")
        assert r.entries == genomes[:2] and r.num_direct == 2

    def test_a_quoted_pattern_is_expanded_in_sorted_order(self, genomes):
        r = resolve_genome_files(["fastas/*.fasta"], "-f")
        assert r.entries == sorted(genomes)

    def test_forms_mix_and_repeats_are_used_once(self, genomes):
        with open("list.txt", "w") as fh:
            fh.write("\n".join(genomes[:2]) + "\n")
        r = resolve_genome_files(["list.txt", genomes[1], genomes[3]], "-f")
        assert r.entries == [genomes[0], genomes[1], genomes[3]]
        assert r.duplicates == [genomes[1]]
        assert r.describe() == "listed in list.txt and given directly"

    def test_list_file_comments_and_blank_lines_are_skipped(self, genomes):
        with open("list.txt", "w") as fh:
            fh.write(f"# my genomes\n\n{genomes[0]}\r\n")
        assert resolve_genome_files(["list.txt"], "-f").entries == [genomes[0]]

    def test_a_pattern_matching_nothing_says_so(self, genomes):
        with pytest.raises(GenomeInputError, match="didn't match any files"):
            resolve_genome_files(["fastas/*.gbk"], "-g")

    def test_a_missing_value_says_so(self, genomes):
        with pytest.raises(GenomeInputError, match="doesn't exist"):
            resolve_genome_files(["nope.fasta"], "-f")

    def test_a_directory_suggests_a_pattern(self, genomes):
        with pytest.raises(GenomeInputError, match=r'fastas/\*\.fasta'):
            resolve_genome_files(["fastas"], "-f")

    def test_whitespace_inside_a_list_entry_is_refused(self, genomes):
        with open("list.txt", "w") as fh:
            fh.write("fastas/g1.fasta extra\n")
        with pytest.raises(GenomeInputError, match="whitespace"):
            resolve_genome_files(["list.txt"], "-f")

    def test_an_empty_list_file_is_refused(self, genomes):
        open("list.txt", "w").close()
        with pytest.raises(GenomeInputError, match="No entries"):
            resolve_genome_files(["list.txt"], "-f")


################################################################################
# resolving -a
################################################################################

class TestResolveAccessions:

    def test_accessions_can_be_given_directly(self, tmp_path):
        r = resolve_accessions(["GCF_000005845.2", "GCA_000009065"])
        assert r.entries == ["GCF_000005845.2", "GCA_000009065"]

    def test_files_and_accessions_mix(self, tmp_path):
        accs = tmp_path / "accs.txt"
        accs.write_text("GCF_000005845.2\nGCA_000009065.1\n")
        r = resolve_accessions([str(accs), "GCF_000006945.2", "GCF_000005845.2"])
        assert r.entries == ["GCF_000005845.2", "GCA_000009065.1", "GCF_000006945.2"]
        assert r.duplicates == ["GCF_000005845.2"]

    def test_something_that_is_neither_says_so(self):
        with pytest.raises(GenomeInputError, match="neither an existing file"):
            resolve_accessions(["accs.txt"])


################################################################################
# the args-level entry points
################################################################################

def test_resolving_args_is_idempotent_and_skips_absent_flags(genomes):
    args = types.SimpleNamespace(ncbi_accessions=None, genbank_files=None,
                                 fasta_files=["fastas/*.fasta"], amino_acid_files=[])
    resolve_genome_input_args(args)
    first = args.fasta_files
    resolve_genome_input_args(args)
    assert args.fasta_files is first
    assert isinstance(first, ResolvedInputs) and len(first) == 4
    assert args.ncbi_accessions is None and args.amino_acid_files == []


def test_input_entries_accepts_every_form(genomes):
    with open("list.txt", "w") as fh:
        fh.write("\n".join(genomes) + "\n")
    assert input_entries("list.txt") == genomes
    assert input_entries(resolve_genome_files(["list.txt"], "-f")) == genomes
    assert input_entries(None) == []


################################################################################
# end to end through the main program's preflight
################################################################################

needs_preflight = pytest.mark.skipif(
    any(shutil.which(b) is None for b in ("muscle", "trimal", "VeryFastTree")),
    reason="preflight needs muscle, trimal and VeryFastTree")


def _preflight(tmp_path, *argv):
    from gtotree.cli.parser import parser
    from gtotree.utils.misc.preflight_checks import preflight_checks
    from gtotree.utils.misc.context import log_file_var
    token = log_file_var.set(log_file_var.get())
    try:
        return preflight_checks(parser().parse_args(
            [*argv, "-H", str(DATA_DIR / "mock.hmm"), "-o", str(tmp_path / "out")]))
    finally:
        log_file_var.reset(token)


@pytest.fixture
def aa_genomes(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "aa").mkdir()
    for i in range(1, 5):
        shutil.copy(DATA_DIR / f"mock-{i}.faa", tmp_path / "aa" / f"mock-{i}.faa")
    return [f"aa/mock-{i}.faa" for i in range(1, 5)]


@needs_preflight
@pytest.mark.parametrize("form", ["list", "direct", "pattern", "mixed"])
def test_every_input_form_yields_the_same_genomes(tmp_path, aa_genomes, form):
    with open("aa.txt", "w") as fh:
        fh.write("\n".join(aa_genomes) + "\n")
    values = {"list": ["aa.txt"],
              "direct": aa_genomes,
              "pattern": ["aa/*.faa"],
              "mixed": ["aa.txt", aa_genomes[0]]}[form]
    _args, run_data = _preflight(tmp_path, "-A", *values)
    assert sorted(gd.id for gd in run_data.all_input_genomes) == \
        ["mock-1", "mock-2", "mock-3", "mock-4"]


@needs_preflight
def test_the_resume_fingerprint_follows_the_genomes_not_the_form(tmp_path, aa_genomes):
    from gtotree.utils.misc.preflight_checks import RESUME, build_fingerprint
    from gtotree.cli.parser import parser
    with open("aa.txt", "w") as fh:
        fh.write("\n".join(sorted(aa_genomes)) + "\n")

    def fp(*values):
        args = parser().parse_args(["-A", *values, "-H", "x"])
        resolve_genome_input_args(args)
        return build_fingerprint(args)

    assert not RESUME.compare(fp("aa.txt"), fp("aa/*.faa"))
    assert RESUME.compare(fp("aa/*.faa"), fp(*sorted(aa_genomes)[:3]))


@needs_preflight
def test_a_bad_input_exits_with_a_message(tmp_path, aa_genomes, capsys):
    with pytest.raises(SystemExit):
        _preflight(tmp_path, "-A", "aa/*.gbk")
    assert "didn't match any files" in capsys.readouterr().out


################################################################################
# the subcommands go through the same resolver
################################################################################

def test_search_annotations_takes_files_and_patterns(tmp_path, aa_genomes):
    from gtotree.utils.misc.general import RunData
    from gtotree.utils.target_search.target_search_cli import build_parser
    from gtotree.utils.target_search.target_search_setup import (
        validate_genome_input_files, _populate_input_genomes)
    args = build_parser().parse_args(["-A", "aa/*.faa", aa_genomes[0], "-p", "PF00001"])
    args = validate_genome_input_files(args)
    rd = RunData()
    _populate_input_genomes(args, rd)
    assert sorted(gd.id for gd in rd.amino_acid_files) == \
        ["mock-1", "mock-2", "mock-3", "mock-4"]


def test_gen_scg_hmms_takes_accessions_and_patterns(tmp_path, aa_genomes):
    from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_cli as cli
    from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_genomes import build_run_data
    args = cli.build_parser().parse_args(["-a", "GCF_000005845.2", "-A", "aa/*.faa"])
    cli.check_args(args)
    cli.resolve_genome_inputs(args)
    rd = build_run_data(args, str(tmp_path / "out"), str(tmp_path / "work"))
    assert [gd.id for gd in rd.ncbi_accs] == ["GCF_000005845.2"]
    assert len(rd.amino_acid_files) == 4


def test_gen_scg_hmms_reports_an_unresolvable_input(tmp_path, aa_genomes):
    from gtotree.utils.hmms.gen_scg_hmms import gen_scg_hmms_cli as cli
    from gtotree.utils.hmms.gen_scg_hmms.gen_scg_hmms_module import GenSCGHMMsError
    args = cli.build_parser().parse_args(["-A", "aa/*.gbk"])
    with pytest.raises(GenSCGHMMsError, match="didn't match any files"):
        cli.resolve_genome_inputs(args)


def test_dl_ncbi_assemblies_takes_accessions_directly(tmp_path):
    from gtotree.utils.ncbi import dl_ncbi_assemblies as dl
    from unittest.mock import patch
    args = dl.build_parser().parse_args(["-a", "GCF_000005845.2", "GCA_000009065.1"])
    with patch("gtotree.utils.ncbi.dl_ncbi_assemblies.get_ncbi_assembly_data"):
        dl.preflight_checks(args)
    assert args.ncbi_accessions.entries == ["GCF_000005845.2", "GCA_000009065.1"]
