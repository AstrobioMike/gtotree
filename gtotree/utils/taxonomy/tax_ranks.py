RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# fill value for a rank that is absent from a lineage
NA = "NA"

REFERENCE_VALUE = "reference genome"

def accession_core(acc):
    """
    Canonical key for an assembly, shared by everything that has to match accessions
    ACROSS sources: GTDB's RS_/GB_ prefixes, the GCA_/GCF_ pairing, and version
    suffixes all normalise away.

        'RS_GCF_000005845.2' -> '000005845'
        'GCA_000005845.3'    -> '000005845'
        'na' / '' / junk     -> ''
    """
    s = str(acc).strip()
    if s.startswith(("RS_", "GB_")):
        s = s[3:]
    if not s.upper().startswith(("GCA_", "GCF_")):
        return ""
    return s.split(".")[0].split("_")[-1]


def accession_core_array(accessions):
    """
    `accession_core` over a whole Arrow string column/array at once, same rules

    For the multi-million-row paths (liveness screening against the NCBI table,
    exclusion over a bulk pool) where calling the scalar version per row through
    `.to_pylist()` costs seconds and a large Python-object footprint. Nulls come back
    as '' just like junk does.
    """
    import pyarrow as pa  # type: ignore
    import pyarrow.compute as pc  # type: ignore

    accessions = pc.cast(accessions, pa.string())

    # fast path: when every value is a canonical NCBI accession (GCF_000005845.2), the
    # core is just characters 4-13, and one regex check + a slice is far cheaper than
    # the general rewrite below. Same answer either way for such values
    canonical = pc.match_substring_regex(accessions,
                                         pattern=r"^GC[AF]_[0-9]{9}\.[0-9]+$")
    if len(accessions) and accessions.null_count == 0 and pc.all(canonical).as_py():
        return pc.utf8_slice_codeunits(accessions, 4, 13)

    s = pc.utf8_trim_whitespace(accessions)
    s = pc.replace_substring_regex(s, pattern=r"^(?:RS_|GB_)", replacement="")
    is_assembly = pc.match_substring_regex(s, pattern=r"^[Gg][Cc][AaFf]_")
    core = pc.replace_substring_regex(s, pattern=r"\..*$", replacement="")   # split(".")[0]
    core = pc.replace_substring_regex(core, pattern=r"^.*_", replacement="")  # split("_")[-1]
    core = pc.if_else(pc.fill_null(is_assembly, False), core, "")
    return pc.fill_null(core, "")


# ---------------------------------------------------------------------------
# taxon-name aliases
# ---------------------------------------------------------------------------

# Names users reach for that aren't the string in the asset
TAXON_ALIASES = {
    "eukarya": "Eukaryota",
    "eukaryote": "Eukaryota",
    "eukaryotes": "Eukaryota",
}


def normalize_taxon_name(taxon):
    """
    Map a user-supplied taxon name onto the name to actually look up
    """
    s = str(taxon).strip()
    return TAXON_ALIASES.get(s.lower(), s)


def rank_index(rank):
    """0 for domain ... 6 for species. Raises ValueError on an unknown rank."""
    r = str(rank).strip().lower()
    if r not in RANKS:
        raise ValueError(
            f"'{rank}' is not one of the supported ranks: {', '.join(RANKS)}")
    return RANKS.index(r)


def validate_derep_rank(wanted_rank, derep_rank):
    """
    The --derep-rank rank must be the SAME as, or FINER than, the rank
    of --wanted-ref-tax

    Returns None if OK, else an explanatory string.
    """
    w = rank_index(wanted_rank)
    d = rank_index(derep_rank)
    if d < w:
        return (f"--derep-rank '{derep_rank}' is a broader rank than "
                f"the target taxon's rank '{wanted_rank}'. The --derep-rank rank "
                f"must be the same or finer (e.g., target 'domain' + "
                f"--derep-rank 'class').")
    return None
