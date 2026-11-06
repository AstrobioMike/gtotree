"""
In-process gene calling (pyrodigal) in place of the prodigal binary.

The contract is "the same genes prodigal calls" for genome-sized input, so where the
binary happens to be installed the outputs are compared against it directly.
"""

import random
import shutil
import subprocess

import pytest  # type: ignore

pytest.importorskip("pyrodigal")

from gtotree.utils.misc.gene_calling import (GenomeTooShort, call_genes,  # noqa: E402
                                             read_contigs)

needs_prodigal = pytest.mark.skipif(shutil.which("prodigal") is None,
                                    reason="the prodigal binary isn't installed")


def _seqs(path):
    out, cur = [], []
    for line in open(path):
        if line.startswith(">"):
            if cur:
                out.append("".join(cur).rstrip("*"))
            cur = []
        else:
            cur.append(line.strip())
    if cur:
        out.append("".join(cur).rstrip("*"))
    return out


def _random_genome(path, n_contigs, contig_len, wrap=60, seed=0):
    rng = random.Random(seed)
    with open(path, "w") as fh:
        for c in range(n_contigs):
            seq = "".join(rng.choice("ACGTACGTGC") for _ in range(contig_len))
            fh.write(f">contig_{c} some description\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")


def test_contig_ids_are_the_first_word_of_the_header(tmp_path):
    path = tmp_path / "g.fna"
    path.write_text(">c1 desc here\nACGT\n>c2\nGG\nTT\n")
    assert read_contigs(str(path)) == [("c1", b"ACGT"), ("c2", b"GGTT")]


def test_too_little_sequence_raises_and_writes_nothing(tmp_path):
    src = tmp_path / "g.fna"
    src.write_text(">a\n" + "ACGT" * 100 + "\n")
    out = tmp_path / "out.faa"
    with pytest.raises(GenomeTooShort):
        call_genes(str(src), str(out))
    assert not out.exists()
    assert not (tmp_path / "out.faa.part").exists()


def test_an_empty_file_is_too_short(tmp_path):
    src = tmp_path / "g.fna"
    src.write_text("")
    with pytest.raises(GenomeTooShort):
        call_genes(str(src), str(tmp_path / "out.faa"))


def test_proteins_carry_no_stop_characters(tmp_path, nt_genomes):
    out = tmp_path / "out.faa"
    assert call_genes(str(nt_genomes[0]), str(out)) > 0
    assert "*" not in out.read_text()


@needs_prodigal
def test_same_proteins_and_cds_as_the_prodigal_binary(tmp_path, nt_genomes):
    for genome in nt_genomes:
        subprocess.run(["prodigal", "-c", "-q", "-i", str(genome),
                        "-a", str(tmp_path / "p.faa"), "-d", str(tmp_path / "p.fna"),
                        "-o", "/dev/null"], check=True)
        call_genes(str(genome), str(tmp_path / "y.faa"), str(tmp_path / "y.fna"))
        assert _seqs(tmp_path / "y.faa") == _seqs(tmp_path / "p.faa")
        assert _seqs(tmp_path / "y.fna") == _seqs(tmp_path / "p.fna")


@needs_prodigal
def test_same_calls_as_prodigal_on_a_multi_contig_genome(tmp_path):
    # several contigs (so the training linkers matter) and long unwrapped lines
    path = tmp_path / "g.fna"
    _random_genome(path, n_contigs=3, contig_len=40_000, wrap=25_000, seed=3)
    subprocess.run(["prodigal", "-c", "-q", "-i", str(path), "-a",
                    str(tmp_path / "p.faa"), "-o", "/dev/null"], check=True)
    call_genes(str(path), str(tmp_path / "y.faa"))
    assert _seqs(tmp_path / "y.faa") == _seqs(tmp_path / "p.faa")


def test_a_gzipped_genome_gives_the_same_genes(tmp_path, nt_genomes):
    import gzip
    gz = tmp_path / "g.fna.gz"
    with open(nt_genomes[0], "rb") as src, gzip.open(gz, "wb") as out:
        out.write(src.read())
    assert read_contigs(str(gz)) == read_contigs(str(nt_genomes[0]))
    call_genes(str(nt_genomes[0]), str(tmp_path / "plain.faa"))
    call_genes(str(gz), str(tmp_path / "gz.faa"))
    assert (tmp_path / "gz.faa").read_text() == (tmp_path / "plain.faa").read_text()


def test_a_gzipped_genbank_is_read_in_place(tmp_path):
    import gzip
    import types
    from Bio.Seq import Seq  # type: ignore
    from Bio.SeqFeature import FeatureLocation, SeqFeature  # type: ignore
    from Bio.SeqRecord import SeqRecord  # type: ignore
    from Bio import SeqIO  # type: ignore
    from gtotree.utils.misc.seqs import (extract_fasta_from_gb,
                                         extract_filter_and_rename_cds_amino_acids_from_gb)
    rec = SeqRecord(Seq("ATGAAATAA" * 3), id="c1", annotations={"molecule_type": "DNA"})
    rec.features.append(SeqFeature(FeatureLocation(0, 9, strand=1), type="CDS",
                                   qualifiers={"translation": ["MK"]}))
    gz = tmp_path / "g.gbk.gz"
    with gzip.open(gz, "wt") as fh:
        SeqIO.write(rec, fh, "genbank")
    rd = types.SimpleNamespace(ready_genome_files_dir=str(tmp_path),
                               genbank_processing_dir=str(tmp_path))
    done, aa_path, num = extract_filter_and_rename_cds_amino_acids_from_gb("g", str(gz), rd)
    assert done and num == 1 and open(aa_path).read() == ">g_1\nMK\n"
    extract_fasta_from_gb("g", str(gz), rd)
    assert (tmp_path / "g.fasta").read_text() == f">g_1\n{rec.seq}\n"
