"""
Resolving the four genome-input flags (`-a`, `-g`, `-f`, `-A`) into their entries

Each flag takes one or more values, in any mix of:

  - genome files themselves (`-f fastas/*.fasta`, expanded by the shell)
  - quoted patterns GToTree expands itself (`-f "fastas/*.fasta"`), which sidesteps the
    shell's limit on argument-list length for very large sets
  - files listing genome files, one path per line (the original form)
  - for `-a`: accessions themselves (`-a GCF_000005845.2`), or files listing them

What a value is gets decided by its extension when that settles it (a `.fasta` is a
genome, no file opened), and otherwise by its first bytes (`>`, `LOCUS`, or gzip that
decompresses to either is a genome; anything else is a list file). So an ordinary glob
over genome files costs no I/O here, and only oddly named values get looked inside.

The guard against the one case that rule gets wrong -- a list of paths saved under a
genome extension -- is a cheap size check: a genome-extension file under
`TINY_FILE_BYTES` that doesn't start like a genome, and whose lines are all existing
paths, is treated as the list file it is.

Every program resolves its flags once, right after parsing (`resolve_genome_input_args`),
and stores a `ResolvedInputs` back on args. Downstream code reads `.entries` and never
needs to know which form was used.

Library code: raises GenomeInputError; each program's CLI layer translates it.
"""

import glob
import gzip
import os
import re

GENOME_EXTENSIONS = (".fa", ".fasta", ".fna", ".faa", ".gb", ".gbk", ".gbff")
TINY_FILE_BYTES = 4096
ACCESSION_PATTERN = re.compile(r"^GC[AF]_\d{9}(\.\d+)?$")
_GLOB_CHARS = re.compile(r"[*?\[]")

# (args dest, flag, kind)
GENOME_INPUT_FLAGS = (
    ("ncbi_accessions", "-a", "accessions"),
    ("genbank_files", "-g", "files"),
    ("fasta_files", "-f", "files"),
    ("amino_acid_files", "-A", "files"),
)


class GenomeInputError(Exception):
    """A genome-input value that can't be resolved"""


class ResolvedInputs:
    """
    One flag's resolved entries, plus where they came from (for reporting)

    Truthy when there are entries, so the existing `if args.fasta_files:` checks keep
    working unchanged.
    """

    def __init__(self, entries, listing_files=(), num_direct=0, duplicates=()):
        self.entries = list(entries)
        self.listing_files = list(listing_files)
        self.num_direct = num_direct
        self.duplicates = list(duplicates)

    def __bool__(self):
        return bool(self.entries)

    def __len__(self):
        return len(self.entries)

    def __iter__(self):
        return iter(self.entries)

    def __eq__(self, other):
        if isinstance(other, ResolvedInputs):
            return self.entries == other.entries
        return NotImplemented

    def __repr__(self):
        return (f"ResolvedInputs({len(self.entries)} entries, "
                f"listing_files={self.listing_files}, num_direct={self.num_direct})")

    def describe(self):
        """
        Where the entries came from, for a run-info line: 'listed in X', 'given
        directly', or both
        """
        parts = []
        if self.listing_files:
            parts.append("listed in " + ", ".join(self.listing_files))
        if self.num_direct:
            parts.append("given directly")
        return " and ".join(parts) if parts else "given"


def input_entries(value):
    """
    The entries behind a genome-input argument, whatever form it's in: a
    ResolvedInputs, a list of raw values, or a single path to a list file (the form
    library callers and older code pass)
    """
    if value is None:
        return []
    if isinstance(value, ResolvedInputs):
        return value.entries
    if isinstance(value, str):
        return read_listing(value)[0]
    return list(value)


################################################################################
# reading
################################################################################

def read_listing(path, flag=None):
    """
    (entries, duplicates) from a list file: one entry per line, blank lines and `#`
    comment lines skipped, whitespace inside an entry refused, duplicates dropped
    with first-seen order kept
    """
    entries = []
    duplicates = []
    seen = set()
    with open(path, errors="replace") as fh:
        for lineno, line in enumerate(fh, start=1):
            entry = line.strip()
            if not entry or entry.startswith("#"):
                continue
            if any(ch.isspace() for ch in entry):
                raise GenomeInputError(
                    f'The list file "{path}"{_flag_note(flag)} has whitespace inside an '
                    f"entry on line {lineno}, which isn't expected (one entry per line):"
                    f"\n\n    {entry}")
            if entry in seen:
                if entry not in duplicates:
                    duplicates.append(entry)
                continue
            seen.add(entry)
            entries.append(entry)
    if not entries:
        raise GenomeInputError(f'No entries were found in "{path}"{_flag_note(flag)}.')
    return entries, duplicates


def _flag_note(flag):
    return f" (passed to `{flag}`)" if flag else ""


def _head_bytes(path, n=64):
    """The first `n` bytes, decompressed if the file is gzipped"""
    with open(path, "rb") as fh:
        head = fh.read(n)
    if head[:2] == b"\x1f\x8b":
        try:
            with gzip.open(path, "rb") as fh:
                head = fh.read(n)
        except (OSError, EOFError):
            return b""
    return head


def starts_like_genome(path):
    """True if the file (gzipped or not) starts the way a FASTA or GenBank file does"""
    head = _head_bytes(path).lstrip()
    return head[:1] == b">" or head[:5] == b"LOCUS"


def has_genome_extension(path):
    name = path.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return name.endswith(GENOME_EXTENSIONS)


def _is_small_list_of_paths(path):
    """A tiny file whose non-blank, non-comment lines are all existing file paths"""
    try:
        with open(path, errors="replace") as fh:
            lines = [ln.strip() for ln in fh
                     if ln.strip() and not ln.strip().startswith("#")]
    except OSError:
        return False
    return bool(lines) and all(os.path.isfile(ln) for ln in lines)


def classify_file(path):
    """
    'genome' or 'listing' for an existing file: by extension when that settles it,
    otherwise by its first bytes
    """
    if has_genome_extension(path):
        if (os.path.getsize(path) < TINY_FILE_BYTES and not starts_like_genome(path)
                and _is_small_list_of_paths(path)):
            return "listing"
        return "genome"
    return "genome" if starts_like_genome(path) else "listing"


################################################################################
# resolving
################################################################################

def _expand(value, flag):
    """
    The existing paths one command-line value stands for: itself, or, for a pattern
    that isn't itself a file, its (sorted) matches
    """
    if os.path.exists(value):
        if os.path.isdir(value):
            raise GenomeInputError(
                f'"{value}" (passed to `{flag}`) is a directory. To use the genomes in '
                f'it, give them with a pattern, e.g. `{flag} "{value.rstrip("/")}/*.fasta"`.')
        return [value]
    if _GLOB_CHARS.search(value):
        matches = sorted(m for m in glob.glob(value) if os.path.isfile(m))
        if not matches:
            raise GenomeInputError(
                f'The pattern "{value}" (passed to `{flag}`) didn\'t match any files.')
        return matches
    raise GenomeInputError(f'"{value}" (passed to `{flag}`) doesn\'t exist.')


def resolve_genome_files(values, flag):
    """
    Resolve one file flag's values (`-g`/`-f`/`-A`) into genome-file paths, as given
    or as expanded. Paths listed inside list files are not checked here; each program
    keeps its own policy for a listed file that's missing.
    """
    entries, listing_files, duplicates = [], [], []
    seen = set()
    num_direct = 0

    def add(entry):
        if entry in seen:
            if entry not in duplicates:
                duplicates.append(entry)
            return False
        seen.add(entry)
        entries.append(entry)
        return True

    for value in values:
        for path in _expand(value, flag):
            if classify_file(path) == "genome":
                if add(path):
                    num_direct += 1
            else:
                listed, dups = read_listing(path, flag)
                listing_files.append(path)
                duplicates.extend(d for d in dups if d not in duplicates)
                for entry in listed:
                    add(entry)

    return ResolvedInputs(entries, listing_files, num_direct, duplicates)


def resolve_accessions(values, flag="-a"):
    """
    Resolve `-a` values: accessions themselves, or files listing them (which may also
    be given as patterns)
    """
    entries, listing_files, duplicates = [], [], []
    seen = set()
    num_direct = 0

    def add(entry):
        if entry in seen:
            if entry not in duplicates:
                duplicates.append(entry)
            return False
        seen.add(entry)
        entries.append(entry)
        return True

    for value in values:
        if not os.path.exists(value) and ACCESSION_PATTERN.match(value):
            if add(value):
                num_direct += 1
            continue
        try:
            paths = _expand(value, flag)
        except GenomeInputError:
            raise GenomeInputError(
                f'"{value}" (passed to `{flag}`) is neither an existing file of '
                "accessions nor an NCBI assembly accession (like GCF_000005845.2).")
        for path in paths:
            listed, dups = read_listing(path, flag)
            listing_files.append(path)
            duplicates.extend(d for d in dups if d not in duplicates)
            for entry in listed:
                add(entry)

    return ResolvedInputs(entries, listing_files, num_direct, duplicates)


def resolve_genome_input_args(args, flags=GENOME_INPUT_FLAGS):
    """
    Replace each present genome-input arg with its ResolvedInputs. Idempotent, and an
    absent flag (None or empty) is left as it is.
    """
    for dest, flag, kind in flags:
        value = getattr(args, dest, None)
        if not value or isinstance(value, ResolvedInputs):
            continue
        values = [value] if isinstance(value, str) else list(value)
        resolver = resolve_accessions if kind == "accessions" else resolve_genome_files
        setattr(args, dest, resolver(values, flag))
    return args


def hash_input_entries(value):
    """
    A content hash of a genome-input arg's entries, in order, for resume fingerprints
    (None when the flag wasn't given)
    """
    import hashlib
    entries = input_entries(value)
    if not entries:
        return None
    return hashlib.sha256("\n".join(entries).encode()).hexdigest()
