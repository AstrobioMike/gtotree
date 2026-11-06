"""
In-process gene calling with pyrodigal, in place of the `prodigal` binary

pyrodigal is a Cython port of Prodigal 2.6.3 (same author as pyhmmer). It calls the
same genes, roughly 2-5x faster, without a subprocess or the extra pass rewriting
prodigal's output to strip its `*` stop characters.

One deliberate difference: prodigal's single-genome mode trains on at most its first
32,000,000 bases, while here training always uses the whole input. For any real
prokaryotic genome (well under 32 Mb) the two are identical; only metagenome-sized
inputs differ, and for those training on all of it is the better model anyway.
"""

import os
import warnings

# pyrodigal warns below 100 kb of training sequence, as prodigal does; the prodigal
# calls this replaced ran with -q. Set once at import (catch_warnings() per call
# isn't thread-safe, and this runs in worker threads)
warnings.filterwarnings("ignore", message=r"sequence should be at least",
                        category=UserWarning, module=r"gtotree\.utils\.misc\.gene_calling")


class GenomeTooShort(Exception):
    """Less sequence than gene calling will train on (prodigal exits 10 for this)"""


def read_contigs(path):
    """
    [(id, sequence bytes)] from a FASTA; the id is the header up to the first
    whitespace, and whitespace inside the sequence is dropped. Reads gzipped input
    directly
    """
    from gtotree.utils.misc.general import open_maybe_gzipped
    contigs = []
    header = None
    chunks = []
    with open_maybe_gzipped(path, "rb") as fh:
        for line in fh:
            if line[:1] == b">":
                if header is not None:
                    contigs.append((header, b"".join(chunks)))
                fields = line[1:].split()
                header = fields[0].decode(errors="replace") if fields else ""
                chunks = []
            elif header is not None:
                chunks.append(b"".join(line.split()))
    if header is not None:
        contigs.append((header, b"".join(chunks)))
    return contigs


def call_genes(nt_path, aa_out_path, nt_out_path=None):
    """
    Call genes on `nt_path` (prodigal's `-c` closed-ends, single mode, trained on the
    whole input) and write the proteins (no `*` stop characters) and optionally the
    CDS nucleotides, each atomically via `.part` + os.replace(). Returns the number
    of genes.

    Raises GenomeTooShort when there's too little sequence to train on.
    """
    import pyrodigal  # type: ignore

    contigs = read_contigs(nt_path)
    if not contigs:
        raise GenomeTooShort("no sequences to train gene calling on")

    finder = pyrodigal.GeneFinder(closed=True)
    try:
        # contigs are joined with prodigal's stop-codon linker for training, as
        # prodigal does, so the minimum is judged on the same length prodigal uses
        finder.train(*(seq for _, seq in contigs))
    except ValueError as err:
        if "at least" in str(err):
            raise GenomeTooShort(str(err)) from err
        raise

    outputs = [aa_out_path] + ([nt_out_path] if nt_out_path else [])
    parts = [p + ".part" for p in outputs]
    num_genes = 0
    try:
        aa_fh = open(parts[0], "w")
        nt_fh = open(parts[1], "w") if nt_out_path else None
        try:
            for contig_id, seq in contigs:
                genes = finder.find_genes(seq)
                num_genes += len(genes)
                genes.write_translations(aa_fh, sequence_id=contig_id,
                                         include_stop=False)
                if nt_fh is not None:
                    genes.write_genes(nt_fh, sequence_id=contig_id)
        finally:
            aa_fh.close()
            if nt_fh is not None:
                nt_fh.close()
        for part, final in zip(parts, outputs):
            os.replace(part, final)
    except BaseException:
        for part in parts:
            try:
                os.remove(part)
            except OSError:
                pass
        raise

    return num_genes
