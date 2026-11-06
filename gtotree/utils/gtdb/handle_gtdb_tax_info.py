import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pyarrow.compute as pc # type: ignore

from gtotree.utils.taxonomy.tax_ranks import RANKS, accession_core
from gtotree.utils.misc.general import newick_safe_label_part
from gtotree.utils.gtdb.get_gtdb_data import gtdb_data_table_path

# the 7 rank columns plus the join key, read from the Parquet asset
_WANTED_COLUMNS = ["ncbi_genbank_assembly_accession"] + list(RANKS)

# extracts the numeric core (the run of digits) from a full assembly accession like
    # "GCA_000739065.1" -> "000739065"
_ACCESSION_CORE_REGEX = r'^[A-Za-z]+_([0-9]+)\.[0-9]+$'


def update_mapping_dict_with_gtdb_tax_info(args, run_data):

    tax_rows = subset_gtdb_info(run_data)

    run_data.tax_info_dict = tax_rows

    for input_acc, ranks in tax_rows.items():
        if input_acc in run_data.mapping_dict:
            continue

        parts = [newick_safe_label_part(ranks[rank.lower()])
                 for rank in args.lineage.split(",") if rank.lower() != "strain"]
        run_data.mapping_dict[input_acc] = "_".join([input_acc] + parts)

    return run_data


def subset_gtdb_info(run_data):
    """
    Pull GTDB lineage info for the remaining NCBI input accessions from the hosted
    Parquet asset.

    Returns a dict: input_acc -> {rank: value for the 7 ranks}, with the species value
    de-duplicated of its leading genus
    """
    lineages = gtdb_lineages_for_accessions(run_data.remaining_ncbi_accs())

    for ranks in lineages.values():
        # removing 'genus' from 'species' otherwise if the user wants genus and species,
        # it would list the genus twice
        genus = ranks["genus"]
        species = ranks["species"]
        if isinstance(species, str):
            ranks["species"] = species.removeprefix(f"{genus} ")

    return lineages


def gtdb_lineages_for_accessions(accessions, gtdb_path=None):
    """
    Map assembly accessions to their GTDB lineages, verbatim from the Parquet asset.

    Accessions with no GTDB row are simply absent from the result

    Returns a dict: input_acc -> {rank: value for the 7 ranks}.

    Shared by the two callers that need this join:
    header re-labelling (`-D`) and the `-H` auto-pick, which links an NCBI-sourced `-w`
    selection back to the GTDB taxonomy the pre-built SCG-sets are built from
    """
    # map each wanted accession to its numeric core for the join
    # a core can be shared by more than one input (the GCA_ and GCF_ of one assembly);
    # each of them gets the lineage, rather than only the last one seen
    core_to_inputs = {}
    for acc in accessions:
        core = accession_core(acc)
        if core:
            core_to_inputs.setdefault(core, []).append(acc)

    if not core_to_inputs:
        return {}

    if gtdb_path is None:
        gtdb_path = gtdb_data_table_path()

    wanted_cores = pa.array(list(core_to_inputs.keys()))

    result = {}
    parquet_file = pq.ParquetFile(gtdb_path)

    for group_index in range(parquet_file.metadata.num_row_groups):
        group = parquet_file.read_row_group(group_index, columns=_WANTED_COLUMNS)

        cores = pc.replace_substring_regex(
            group.column("ncbi_genbank_assembly_accession"),
            _ACCESSION_CORE_REGEX, r'\1')
        matched = group.filter(pc.is_in(cores, value_set=wanted_cores))

        if matched.num_rows == 0:
            continue

        gb_matched = matched.column("ncbi_genbank_assembly_accession").to_pylist()
        rank_cols = {rank: matched.column(rank).to_pylist() for rank in RANKS}

        for i, gb_acc in enumerate(gb_matched):
            for input_acc in core_to_inputs.get(accession_core(gb_acc), ()):
                result[input_acc] = {rank: rank_cols[rank][i] for rank in RANKS}

    return result
