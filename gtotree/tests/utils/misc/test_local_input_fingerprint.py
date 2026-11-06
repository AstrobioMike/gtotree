"""
The main run's resume fingerprint used to hash only the -g/-f/-A listing files and the
-H argument string, so editing or replacing a listed genome file, or a local HMM file,
in place still resumed against the old results.
"""

import os

from gtotree.utils.misc.preflight_checks import RESUME, build_fingerprint


class _Args:
    def __init__(self, **kw):
        for name in ("wanted_ref_tax", "ncbi_accessions", "genbank_files", "fasta_files",
                     "amino_acid_files", "mapping_file", "target_pfams_file",
                     "target_kos_file", "exclusion_list", "target_rank", "target_domain",
                     "gtdb_section"):
            setattr(self, name, None)
        self.hmm = "Universal"
        self.derep_rank = "auto"
        self.source = "gtdb"
        self.ncbi_section = "both"
        self.add_gtdb_tax = self.add_ncbi_tax = False
        self.lineage = ""
        self.__dict__.update(kw)


def _listing(tmp_path, *names):
    paths = []
    for name in names:
        p = tmp_path / name
        if not p.exists():
            p.write_text(">c\nMKV\n")
        paths.append(str(p))
    listing = tmp_path / "aa.txt"
    listing.write_text("\n".join(paths) + "\n")
    return str(listing)


def test_an_edited_genome_file_invalidates_a_resume(tmp_path):
    listing = _listing(tmp_path, "g1.faa", "g2.faa")
    before = build_fingerprint(_Args(amino_acid_files=listing))

    target = tmp_path / "g1.faa"
    target.write_text(">c\nMKVLLLLLLLL\n")
    os.utime(target, (1, 1))
    after = build_fingerprint(_Args(amino_acid_files=listing))

    assert RESUME.compare(before, after)


def test_untouched_genome_files_still_resume(tmp_path):
    listing = _listing(tmp_path, "g1.faa", "g2.faa")
    assert not RESUME.compare(build_fingerprint(_Args(amino_acid_files=listing)),
                              build_fingerprint(_Args(amino_acid_files=listing)))


def test_an_edited_local_hmm_file_invalidates_a_resume(tmp_path):
    hmm = tmp_path / "mine.hmm"
    hmm.write_text("HMMER3/f\nNAME a\n//\n")
    before = build_fingerprint(_Args(hmm=str(hmm)))
    hmm.write_text("HMMER3/f\nNAME b\n//\n")
    after = build_fingerprint(_Args(hmm=str(hmm)))

    assert RESUME.compare(before, after)
