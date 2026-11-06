from gtotree.utils.misc.general import (write_run_data,
                                   required_count,
                                   run_pooled_stage)
from gtotree.utils.misc.messaging import (report_message,
                                     report_processing_stage,
                                     report_SCG_set_filtering_update)
from gtotree.utils.misc.seqs import check_target_SCGs_have_seqs, filter_seqs_by_length
from gtotree.utils.misc.stages import PipelineStage, SCGRemovalStage
from gtotree.utils.misc.summary_info import write_SCG_info_table

def count_genome_hits_in_remaining_SCGs(run_data):
    """
    Set each genome's `num_SCG_hits_after_filtering` to the number of still-remaining
    SCG-sets whose length-filtered file holds a hit for it
    """
    genome_dict = {gd.id: gd for gd in run_data.get_all_input_genomes_for_filtering()}
    counts = dict.fromkeys(genome_dict, 0)

    for scg in run_data.get_all_SCG_targets_remaining():
        path = run_data.found_SCG_seqs_dir + f"/{scg.id}-gene-filtered{run_data.general_ext}"
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    genome_id = line[1:].strip()
                    if genome_id in counts:
                        counts[genome_id] += 1

    for genome_id, count in counts.items():
        genome_dict[genome_id].num_SCG_hits_after_filtering = count


def filter_genes(args, run_data):

    report_processing_stage("filter-genes", run_data)
    cutoff = f"{run_data.seq_length_cutoff * 100:.0f}"
    in_genomes_cutoff_for_report = f"{args.gene_representation_cutoff * 100:.0f}"

    message = (f"Keeping genes with lengths within {cutoff}% of the median for each gene set (controlled by `-c`), "
               f"and keeping gene sets with hits in at least {in_genomes_cutoff_for_report}% of "
               f"the currently retained genomes (controlled by `-r`).")
    report_message(message, ii="    ", si="    ", width=80)

    if run_data.stage_is_complete(PipelineStage.FILTER_GENES):
        write_SCG_info_table(run_data)
        report_SCG_set_filtering_update(run_data)
        return run_data

    scgs_to_filter = run_data.get_all_SCG_targets_remaining_but_not_filtered()

    if len(scgs_to_filter) > 0:

        def worker(scg, run_data):
            path = run_data.found_SCG_seqs_dir + f"/{scg.id}{run_data.general_ext}"
            # returns the ids of genomes whose hit survived length-filtering; the
            # length-filtered FASTA is written alongside as a side effect (per-SCG
            # file only -- no shared state touched here)
            return filter_seqs_by_length(path, run_data.seq_length_cutoff)

        def apply_result(scg, genomes_with_hits, run_data):
            scg.gene_length_filtered = True
            scg.num_genomes_after_length_filtering = len(genomes_with_hits)

        run_data = run_pooled_stage(scgs_to_filter, worker, apply_result, args, run_data)

        write_run_data(run_data)

    run_data = check_target_SCGs_have_seqs(run_data,
                                          f"-gene-filtered{run_data.general_ext}",
                                          SCGRemovalStage.GENE_FILTER)

    total_genomes_remaining = len(run_data.get_all_input_genomes_for_filtering())
    min_genomes_required = required_count(total_genomes_remaining,
                                          args.gene_representation_cutoff)

    removed_any = False
    for scg in run_data.get_all_SCG_targets_remaining():
        count = getattr(scg, 'num_genomes_after_length_filtering', 0)
        if count < min_genomes_required:
            scg.mark_removed(
                f"too few genomes with hits ({count} < {min_genomes_required} required)",
                SCGRemovalStage.GENE_FILTER)
            removed_any = True

    # counted only now, over the SCG-sets that survived `-r` too: `-G` is "hits to the
    # target SCGs that end up in the final tree", and counting before `-r` (or only
    # the SCGs filtered in this invocation, on a resume) inflated it
    count_genome_hits_in_remaining_SCGs(run_data)

    write_SCG_info_table(run_data)

    run_data.mark_stage_complete(PipelineStage.FILTER_GENES)
    write_run_data(run_data)

    report_SCG_set_filtering_update(run_data)

    return run_data
