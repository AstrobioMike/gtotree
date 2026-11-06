"""
Taxon-based genome/accession selection from the GTDB / NCBI Parquet assets.

The filter/select piece of the taxonomy toolkit: resolve a user-supplied taxon to
a (canonical, rank) pair and pull the rows under it from the hosted Parquet assets
(bit's ncbi-data.parquet / gtdb-data.parquet). Column names are the only thing that
differs between sources; that difference is isolated in SourceSpec / SOURCES.
"""

import os
import pyarrow.compute as pc # type: ignore
import pyarrow.parquet as pq # type: ignore

import pyarrow as pa # type: ignore

from gtotree.utils.taxonomy.tax_ranks import (RANKS, NA, REFERENCE_VALUE, accession_core,
                                              accession_core_array,
                                              normalize_taxon_name, rank_index)


class SourceSpec:
    """Per-source column names. The only thing that differs between GTDB and NCBI."""

    def __init__(self, name, acc_col, rep_filter, quality_cols,
                 ref_col, default_reps_only, level_col=None, contig_col=None,
                 size_col=None, size_fallback_col=None):
        self.name = name
        self.acc_col = acc_col
        self.rep_filter = rep_filter
        self.quality_cols = quality_cols
        self.ref_col = ref_col
        self.default_reps_only = default_reps_only
        self.level_col = level_col
        self.contig_col = contig_col
        self.size_col = size_col
        self.size_fallback_col = size_fallback_col


SOURCES = {
    "gtdb": SourceSpec(
        name="gtdb",
        acc_col="ncbi_genbank_assembly_accession",
        rep_filter=("gtdb_representative", "t"),
        quality_cols=("checkm2_completeness", "checkm2_contamination"),
        ref_col="ncbi_refseq_category",
        default_reps_only=True,
        size_col="genome_size",
        contig_col="contig_count",
    ),
    "ncbi": SourceSpec(
        name="ncbi",
        acc_col="assembly_accession",
        rep_filter=("refseq_category", REFERENCE_VALUE),
        quality_cols=("checkm_completeness", "checkm_contamination"),
        ref_col="refseq_category",
        default_reps_only=False,
        level_col="assembly_level",
        contig_col="contig_count",
        size_col="genome_size_ungapped",
        size_fallback_col="genome_size",
    ),
}


class TaxonNotFound(Exception):
    pass


class AmbiguousTaxon(Exception):
    """The taxon name exists at more than one rank; caller must disambiguate."""

    def __init__(self, taxon, ranks_found):
        self.taxon = taxon
        self.ranks_found = ranks_found
        super().__init__(
            f"'{taxon}' occurs at more than one rank ({', '.join(ranks_found)}); "
            f"specify which rank is wanted")


class CrossDomainTaxon(Exception):
    """
    The taxon name occurs in more than one domain; caller must disambiguate.

    e.g., `Bacillus` is both a bacterial and eukaryotic genus
    """

    def __init__(self, taxon, domains_found):
        self.taxon = taxon
        self.domains_found = list(domains_found)
        example = self.domains_found[0] if self.domains_found else "<domain>"
        super().__init__(
            f"'{taxon}' occurs in more than one domain "
            f"({', '.join(self.domains_found)}). Specify which domain is wanted with "
            f"`--target-domain` (e.g., `--target-domain {example}`).")


# ---------------------------------------------------------------------------
# resolving a user-supplied taxon
# ---------------------------------------------------------------------------

def find_ranks_for_taxon(path, taxon):
    """
    Which of the 7 ranks contain `taxon`, and which domains it spans. Case-insensitive.

    Returns (canonical_name, [ranks], [domains]). Reads one column at a time, so this
    stays cheap even on the 4M-row NCBI table. `domains` is the sorted list of distinct
    assigned domains the name appears in (across whatever rank(s) it occurs at); it is
    what lets resolution catch a name shared across domains.
    """
    target = normalize_taxon_name(taxon).lower()
    canonical = None
    found = []
    for rank in RANKS:
        name = _rank_names(path, rank).get(target)
        if name is not None:
            canonical = name
            found.append(rank)
    if canonical is None:
        raise TaxonNotFound(f"'{taxon}' doesn't exist at any rank in this source")

    domains = _domains_for_taxon(path, canonical, found)
    return canonical, found, domains


_RANK_NAMES_CACHE = {}


def _rank_names(path, rank):
    """
    {lowercased name: name} for one rank's distinct assigned names, keeping the
    first-seen spelling when two differ only in case. Cached per (path, size,
    mtime): reading and lowercasing a rank column is the slow part of resolving a
    taxon, and a run resolves one per `-w` taxon (and again for each size check)
    """
    st = os.stat(path)
    key = (os.path.abspath(path), st.st_size, st.st_mtime_ns, rank)
    names = _RANK_NAMES_CACHE.get(key)
    if names is None:
        names = {}
        uniq = pc.unique(pq.read_table(path, columns=[rank]).column(rank))
        for name in uniq.to_pylist():
            if name is not None and name != NA:
                names.setdefault(name.lower(), name)
        _RANK_NAMES_CACHE[key] = names
    return names


def _domains_for_taxon(path, canonical, ranks):
    """
    The sorted distinct assigned domains a resolved taxon appears in.

    One filtered scan of the `domain` column per rank the name occupies (usually one).
    A domain-rank taxon is its own domain, handled without a scan.
    """
    domain_rank = RANKS[0]
    domains = set()
    for rank in ranks:
        if rank == domain_rank:
            domains.add(canonical)
            continue
        tab = pq.read_table(path, columns=[domain_rank],
                            filters=[(rank, "=", canonical)])
        for d in pc.unique(tab.column(domain_rank)).to_pylist():
            if d and d != NA:
                domains.add(d)
    return sorted(domains)


def resolve_taxon(path, taxon, rank=None, domain=None):
    """
    Resolve a user-supplied taxon (+ optional explicit rank and/or domain) to
    (canonical, rank, domain).

    Raises AmbiguousTaxon if the name lives at multiple ranks and no `rank` was given
    (e.g. a name used as both an order and a family).

    Raises CrossDomainTaxon if the name occurs in more than one domain and no `domain`
    was given (e.g. `Bacillus`, both a bacterial and a eukaryotic genus). Selecting on
    the bare name would mix domains into one tree.

    `rank` and `domain` are independent disambiguators and may be given together. A
    `domain` is normalized through the same alias table as `-w` (so 'eukarya' ->
    'Eukaryota'), and must be one the taxon actually occurs in.

    The returned `domain` is the taxon's sole domain (or the one selected), or None
    when the taxon has no assigned domain at all (viral/metagenome names).
    """
    canonical, found, domains = find_ranks_for_taxon(path, taxon)

    if rank:
        r = str(rank).strip().lower()
        rank_index(r)
        if r not in found:
            raise TaxonNotFound(
                f"'{canonical}' exists at rank(s) {', '.join(found)}, not '{r}'")
        resolved_rank = r
    elif len(found) > 1:
        raise AmbiguousTaxon(canonical, found)
    else:
        resolved_rank = found[0]

    # `domains` spans every rank the name occurs at. Once `-r` has picked one of
    # several, only that rank's domains matter: e.g. GTDB's JANQNQ01 is an archaeal
    # family and a bacterial genus, and `-r family` already says which is meant
    if len(found) > 1:
        domains = _domains_for_taxon(path, canonical, [resolved_rank])

    resolved_domain = _resolve_domain_choice(canonical, domains, domain)

    return canonical, resolved_rank, resolved_domain


def _resolve_domain_choice(canonical, domains, domain):
    """
    Reconcile the domains a taxon spans with an optional user --target-domain.

    - no --target-domain: sole domain if there's exactly one, None if there are none
      (domain-less viral/metagenome names), CrossDomainTaxon if there's more than one.
    - --target-domain given: it must be one the taxon actually occurs in; returns it
      (as the asset spells it). Aliases are honored via normalize_taxon_name.
    """
    if domain is None:
        if len(domains) > 1:
            raise CrossDomainTaxon(canonical, domains)
        return domains[0] if domains else None

    wanted = normalize_taxon_name(domain).strip().lower()
    for d in domains:
        if d.lower() == wanted:
            return d
    if not domains:
        raise TaxonNotFound(
            f"'{canonical}' has no assigned domain, so it can't be scoped with "
            f"--target-domain.")
    raise TaxonNotFound(
        f"'{canonical}' doesn't occur in domain '{domain}'. It occurs in: "
        f"{', '.join(domains)}.")


# ---------------------------------------------------------------------------
# selection
# ---------------------------------------------------------------------------

UNASSIGNED = (NA, "")


def prefix_mask(acc_col, prefixes):
    mask = None
    for p in prefixes:
        m = pc.starts_with(acc_col, p)
        mask = m if mask is None else pc.or_(mask, m)
    return mask


def assigned_mask(col, assigned=True):
    """
    Mask for rows whose value in `col` is (or isn't) an actually-assigned taxon.

    Null, 'NA' and '' all mean "nothing assigned at this rank"
    """
    unassigned = pa.array(list(UNASSIGNED), type=col.type)
    is_assigned = pc.and_(pc.is_valid(col),
                          pc.invert(pc.fill_null(pc.is_in(col, value_set=unassigned),
                                                 False)))
    return is_assigned if assigned else pc.invert(is_assigned)


def rep_predicate(spec, reps_only):
    """
    The (column, value) representatives predicate `reps_only` asks for, or None

    True -> the source's default predicate; a 2-tuple -> used as given; falsy -> None
    """
    if not reps_only:
        return None
    if isinstance(reps_only, tuple):
        return reps_only
    return spec.rep_filter


def select(path, source, rank, taxon, reps_only=False, columns=None,
           accession_prefixes=None, assembly_levels=None, domain=None):
    """
    All rows under `taxon` at `rank`. Returns a pyarrow Table

    `domain`, when given, additionally scopes to that domain -- needed when the taxon
    name is shared across domains (e.g. the genus `Bacillus`) and the caller has
    resolved which domain is wanted. Ignored when the target rank IS domain (the rank
    filter already pins it).

    `reps_only` is normally a bool (True -> the source's own representatives
    predicate), but can also be an explicit `(column, value)` predicate, for sources
    with more than one notion of "representative" (GTDB: its own species reps vs
    RefSeq reference genomes, see `rep_predicate`)
    """
    spec = SOURCES[source]
    filters = [(rank, "=", taxon)]
    if domain and rank != RANKS[0]:
        filters.append((RANKS[0], "=", domain))
    rep = rep_predicate(spec, reps_only)
    if rep:
        filters.append((rep[0], "=", rep[1]))
    if assembly_levels and spec.level_col:
        filters.append((spec.level_col, "in", set(assembly_levels)))

    cols = columns or [spec.acc_col]
    # always need the accession back
    if spec.acc_col not in cols:
        cols = [spec.acc_col] + list(cols)

    tab = pq.read_table(path, columns=cols, filters=filters)

    if accession_prefixes:
        tab = tab.filter(prefix_mask(tab.column(spec.acc_col), accession_prefixes))

    return tab


# ---------------------------------------------------------------------------
# diagnosing an empty candidate pool
# ---------------------------------------------------------------------------
#
# The pool filters below are Parquet pushdown predicates (or, for prefixes, a mask
# applied right after the read), so nothing downstream ever learns what any one of
# them cost. Counting in stages on the HOT path would mean giving that pushdown up,
# which isn't worth it to serve an error message.
#
# So instead: when a selection comes back empty, peel the optional filters back one
# at a time and see which one was load-bearing. Cost is irrelevant here (we are
# already on our way to exiting), and it's a handful of single-column reads against
# an already-narrow taxon slice.

POOL_FILTERS = ("assembly_levels", "accession_prefixes", "reps_only")


def count_pool(path, source, rank, taxon, domain=None, reps_only=False,
               accession_prefixes=None, assembly_levels=None):
    """How many rows the pool `select()` would build holds. One narrow read."""
    tab = select(path, source, rank, taxon, reps_only=reps_only,
                 accession_prefixes=accession_prefixes,
                 assembly_levels=assembly_levels, domain=domain)
    return tab.num_rows


def present_assembly_levels(path, source, rank, taxon, domain=None, reps_only=False,
                            accession_prefixes=None):
    """
    The assembly_level values a taxon actually HAS, sorted best-to-worst.

    For the very common "you asked for a level this taxon doesn't have" case, so the
    message can say what is there instead of only what isn't. Empty list for a source
    with no assembly-level column (GTDB).
    """
    spec = SOURCES[source]
    if not spec.level_col:
        return []

    tab = select(path, source, rank, taxon, reps_only=reps_only,
                 columns=[spec.level_col],
                 accession_prefixes=accession_prefixes, domain=domain)
    col = tab.column(spec.level_col)
    if len(col) == 0:
        return []

    present = pc.unique(col.filter(assigned_mask(col))).to_pylist()
    return sorted(present, key=lambda v: LEVEL_REPORT_ORDER.get(str(v).strip().lower(),
                                                               len(LEVEL_REPORT_ORDER)))


#: best-to-worst, so a "present for it: ..." list reads the way a person would say it
LEVEL_REPORT_ORDER = {
    "complete genome": 0,
    "chromosome": 1,
    "scaffold": 2,
    "contig": 3,
}


def diagnose_empty_pool(path, source, rank, taxon, domain=None, reps_only=False,
                        accession_prefixes=None, assembly_levels=None):
    """
    Work out why the candidate pool for `taxon` came back empty.

    Returns (n_unfiltered, culprits, present_levels):

        n_unfiltered  -- rows under the taxon with NONE of the optional pool filters
                         applied. Zero means the taxon is genuinely empty here and no
                         filter is to blame.
        culprits      -- the POOL_FILTERS names that, dropped on their own, bring the
                         pool back. Usually exactly one. EMPTY when no single filter
                         explains it, which is the signal that a COMBINATION did it
                         and the caller should list them all rather than accuse one.
        present_levels-- the assembly levels the taxon does have, but only when
                         assembly_levels is implicated; [] otherwise.

    `domain` is deliberately not treated as a droppable filter: it is a disambiguator
    resolved FROM the taxon, not something the user narrowed with, so it can't be the
    thing that emptied the pool.
    """
    active = {
        "assembly_levels": assembly_levels or None,
        "accession_prefixes": accession_prefixes or None,
        "reps_only": reps_only or None,
    }

    def count(**overrides):
        kwargs = {
            "reps_only": ((active["reps_only"] or False)
                          if overrides.get("reps_only", True) else False),
            "accession_prefixes": (active["accession_prefixes"]
                                   if overrides.get("accession_prefixes", True) else None),
            "assembly_levels": (active["assembly_levels"]
                                if overrides.get("assembly_levels", True) else None),
        }
        return count_pool(path, source, rank, taxon, domain=domain, **kwargs)

    n_unfiltered = count(reps_only=False, accession_prefixes=False,
                         assembly_levels=False)

    culprits = []
    if n_unfiltered:
        for name in POOL_FILTERS:
            if active[name] and count(**{name: False}):
                culprits.append(name)

    present_levels = []
    if culprits == ["assembly_levels"]:
        present_levels = present_assembly_levels(
            path, source, rank, taxon, domain=domain,
            reps_only=active["reps_only"] or False,
            accession_prefixes=active["accession_prefixes"])

    return n_unfiltered, culprits, present_levels


# ---------------------------------------------------------------------------
# liveness screening (suppressed / removed / version-drifted assemblies)
# ---------------------------------------------------------------------------

_LIVE_CORES_CACHE = {}


def live_accession_cores(ncbi_table_path):
    """
    The assembly core accs currently present in the NCBI assembly summary, as a unique
    Arrow array (for `pc.is_in`)

    NCBI drops suppressed/removed assemblies from the assembly summary files, so
    ABSENCE == suppressed / removed. GTDB is a snapshot and can still have them

    Computed in Arrow and cached per (path, size, mtime): this used to build a Python
    set of all ~4.4M cores (~1 GB, ~9 s) on every call, and `select_all_domains` and
    multi-taxon `-w` runs call it once per domain/taxon
    """
    st = os.stat(ncbi_table_path)
    key = (os.path.abspath(ncbi_table_path), st.st_size, st.st_mtime_ns)
    cached = _LIVE_CORES_CACHE.get(key)
    if cached is None:
        col = pq.read_table(ncbi_table_path,
                            columns=["assembly_accession"]).column("assembly_accession")
        # the scalar version skipped empty/null accessions entirely
        col = pc.filter(col, pc.fill_null(pc.not_equal(col, ""), False))
        cached = pc.unique(accession_core_array(col))
        _LIVE_CORES_CACHE.clear()
        _LIVE_CORES_CACHE[key] = cached
    return cached


def live_mask(accessions, live_cores):
    """
    Boolean mask: which `accessions` are still live at NCBI (non-empty, and their
    core is among `live_cores` from `live_accession_cores`)
    """
    present = pc.fill_null(pc.not_equal(accessions, ""), False)
    return pc.and_(present, pc.is_in(accession_core_array(accessions),
                                     value_set=live_cores))


def filter_table_by_liveness(tab, acc_col, ncbi_table_path):
    """
    Drop rows whose assembly is no longer live at NCBI. Returns (kept, num_dropped)
    """
    before = tab.num_rows
    tab = tab.filter(live_mask(tab.column(acc_col),
                               live_accession_cores(ncbi_table_path)))
    return tab, before - tab.num_rows
