import os
import easygrid
import re
import numpy as np
import sys

# Stage names
DETECT_BARCODES = 'detect_barcodes'
PREPROCESS_CDS = 'preprocess'
SWEEP_CAPTURE_RATE_MOI = 'capture_rate_moi_sweep'
CAPTURE_RATE_MOI_FIGURES = 'capture_rate_moi_figures'
INITIAL_SCREENS_MATCHED_DEG = 'initial_screens_matched_deg'
INITIAL_SCREENS_ANALYSIS = 'initial_screens_analysis'
FACS_SWAP_RATE = 'facs_swap'
FACS_SWAP_RATE_FIGURE = 'facs_swap_figures'
SIMULATE_SWAP_RATES = 'swap_sim'
SWAP_RATE_SIM_FIGURE = 'swap_sim_figure'
GENERATE_DEG_JOBS = 'pairwise_deg_setup'
PAIRWISE_DEG = 'pairwise_deg'
BARCODE_ENRICHMENT = 'barcode_enrichment'
SUMMARIZE_BARCODE_ENRICHMENT = 'barcode_enrichment_figures'
MAIN_TEXT_TSNE = 'main_text_tsne'
REVIEWER_FIGURES = 'reviewer_figures'

# Prep directories
easygrid.mkdir('temp_data')
easygrid.mkdir('data')
easygrid.mkdir('figures')
easygrid.mkdir('tables')
easygrid.mkdir('supplemental_figures')
easygrid.mkdir('reviewer_figures')
easygrid.mkdir('diagnostic_plots')
easygrid.mkdir('temp_data/ko_barcodes')
easygrid.mkdir('temp_data/ko_barcodes/cropseq')
easygrid.mkdir('temp_data/ko_barcodes/initial_arrayed')
easygrid.mkdir('temp_data/ko_barcodes/initial_pooled')
easygrid.mkdir('temp_data/pairwise_deg/')
easygrid.mkdir('temp_data/pairwise_deg/raw_data')
easygrid.mkdir('temp_data/pairwise_deg/informative/')
easygrid.mkdir('temp_data/pairwise_deg/raw_data/singles')
easygrid.mkdir('temp_data/pairwise_deg/informative/singles/')
easygrid.mkdir('temp_data/external_single_cell_experiments_cds/')

####################################
# KO barcode detection
####################################
pipeline = easygrid.JobManager('.easygrid_hacks')

barcode_detection_commands = [line.strip() for line in open('barcode_detection.jobs')]
outputs = [re.search('[-]o (.+[.]txt) --whitelist', command).group(1).strip() for command in barcode_detection_commands]

for command, output in zip(barcode_detection_commands, outputs):
    pipeline.add(command, name=DETECT_BARCODES, memory='10G', outputs=[output])

####################################
# Merging RNA-seq and KO barcodes (CROP-seq)
####################################
metadata_file = 'data/cropseq.sample_metadata.txt'
cds = 'temp_data/aggregated_cds.rds'
pdata = 'temp_data/aggregated_cds_pdata.txt'
guide_gene_associations = 'data/guide_gene_associations.txt'
barcode_qc_plot = 'supplemental_figures/barcode_enrichment_qc.png'

command = "Rscript preprocess_cfg.R %s %s %s --guide_metadata %s --barcode_enrichment_qc_plot %s" % (metadata_file, cds, pdata, guide_gene_associations, barcode_qc_plot)

pipeline.add(command, name=PREPROCESS_CDS, memory='50G', dependencies=[DETECT_BARCODES], outputs=[cds, pdata, barcode_qc_plot])

####################################
# Merging RNA-seq and KO barcodes (CROP-seq with no barcode enrichment)
####################################
metadata_file_no_enrichment = 'data/cropseq.sample_metadata.no_enrichment.txt'
cds_no_enrichment = 'temp_data/aggregated_cds.no_enrichment.rds'
pdata_no_enrichment = 'temp_data/aggregated_cds_pdata.no_enrichment.txt'
guide_gene_associations = 'data/guide_gene_associations.txt'
barcode_qc_plot_no_enrichment = 'diagnostic_plots/barcode_enrichment_qc.no_enrichment.png'

command = "Rscript preprocess_cfg.R %s %s %s --guide_metadata %s --barcode_enrichment_qc_plot %s --ko_assignment_reads_threshold 3" % (metadata_file_no_enrichment, cds_no_enrichment, pdata_no_enrichment, guide_gene_associations, barcode_qc_plot_no_enrichment)

pipeline.add(command, name=PREPROCESS_CDS, memory='50G', dependencies=[DETECT_BARCODES], outputs=[cds_no_enrichment, pdata_no_enrichment, barcode_qc_plot_no_enrichment])

####################################
# Merging RNA-seq and KO barcodes (initial screens)
####################################
metadata_file_initial_screens = 'data/initial_screens.sample_metadata.txt'
cds_initial_screens = 'temp_data/aggregated_cds.initial_screens.rds'
pdata_initial_screens = 'temp_data/aggregated_cds_pdata.initial_screens.txt'
guide_gene_associations_initial_screens = 'data/barcode_gene_associations.initial_10_target.txt'
barcode_qc_plot_initial_screens = 'supplemental_figures/barcode_enrichment_qc.initial_10_target.png'

# Note running with low size factor cluster removal turned off
command = "Rscript preprocess_cfg.R %s %s %s --guide_metadata %s --barcode_enrichment_qc_plot %s --no_size_factor_filter --aggregated" % (metadata_file_initial_screens, cds_initial_screens, pdata_initial_screens, guide_gene_associations_initial_screens, barcode_qc_plot_initial_screens)

pipeline.add(command, name=PREPROCESS_CDS, memory='50G', dependencies=[DETECT_BARCODES], outputs=[cds_initial_screens, pdata_initial_screens, barcode_qc_plot_initial_screens])

##################################
# Munge external datasets
##################################
command_template = 'Rscript munge_adamson_to_cds.R --matrix %s --genes %s --barcodes %s --assignments %s --output_file %s --treatment %s'

adamson_upr_cds = 'temp_data/external_single_cell_experiments_cds/adamson.upr.rds'
matrix = 'data/external_single_cell_experiments/adamson/upr/GSM2406681_10X010_matrix.mtx.txt.gz'
genes = 'data/external_single_cell_experiments/adamson/upr/GSM2406681_10X010_genes.tsv.gz'
barcodes = 'data/external_single_cell_experiments/adamson/upr/GSM2406681_10X010_barcodes.tsv.gz'
assignments = 'data/external_single_cell_experiments/adamson/upr/GSM2406681_10X010_cell_identities.csv.gz'
command = command_template % (matrix, genes, barcodes, assignments, adamson_upr_cds, 'adamson.upr')
pipeline.add(command, name=PREPROCESS_CDS, memory='50G', dependencies=[DETECT_BARCODES], outputs=[adamson_upr_cds])

adamson_tf_pilot_cds = 'temp_data/external_single_cell_experiments_cds/adamson.tf_pilot.rds'
matrix = 'data/external_single_cell_experiments/adamson/tf_pilot/GSM2406675_10X001_matrix.mtx.txt.gz'
genes = 'data/external_single_cell_experiments/adamson/tf_pilot/GSM2406675_10X001_genes.tsv.gz'
barcodes = 'data/external_single_cell_experiments/adamson/tf_pilot/GSM2406675_10X001_barcodes.tsv.gz'
assignments = 'data/external_single_cell_experiments/adamson/tf_pilot/GSM2406675_10X001_cell_identities.csv.gz'
command = command_template % (matrix, genes, barcodes, assignments, adamson_tf_pilot_cds, 'adamson.tf_pilot')
pipeline.add(command, name=PREPROCESS_CDS, memory='50G', dependencies=[DETECT_BARCODES], outputs=[adamson_tf_pilot_cds])

##################################
# Capture rate estimates for CROP-seq
##################################
easygrid.mkdir('temp_data/moi_capture_rate_parameter_sweep')

log_likelihoods = 'temp_data/moi_capture_rate_parameter_sweep/crop_seq_log_likelihoods.txt'
command = 'python moi_capture_rate_generative_model.py %s %s --sample_column condition' % (pdata, log_likelihoods)
pipeline.add(command, name=SWEEP_CAPTURE_RATE_MOI, dependencies=[PREPROCESS_CDS], outputs=[log_likelihoods], memory='10G')

log_likelihood_plot = 'supplemental_figures/crop_seq_capture_rate_moi_estimate.png'
command = 'Rscript moi_capture_rate_figures.R'
pipeline.add(command, name=CAPTURE_RATE_MOI_FIGURES, dependencies=[SWEEP_CAPTURE_RATE_MOI], outputs=[log_likelihood_plot], memory='10G')

##################################
# Capture rate estimates for CROP-seq (no enrichment)
##################################
log_likelihoods = 'temp_data/moi_capture_rate_parameter_sweep/crop_seq_log_likelihoods.no_enrichments.txt'
command = 'python moi_capture_rate_generative_model.py %s %s --sample_column condition' % (pdata_no_enrichment, log_likelihoods)
pipeline.add(command, name=SWEEP_CAPTURE_RATE_MOI, dependencies=[PREPROCESS_CDS], outputs=[log_likelihoods], memory='10G')
# note that the capture rate here was just determined by examining this output file

##################################
# Matched DEG tests for arrayed vs. pooled screens
##################################
easygrid.mkdir('temp_data/initial_screen_analysis_deg')

for seed in range(0, 10):
    output_file = 'temp_data/initial_screen_analysis_deg/matched_deg.seed%s.rds' % seed

    command = 'Rscript do_arrayed_pooled_matched_deg.R %s %s --seed %s' % (cds_initial_screens, output_file, seed)
    pipeline.add(command, name=INITIAL_SCREENS_MATCHED_DEG, dependencies=[PREPROCESS_CDS], outputs=[output_file], memory='50G')

###################################
# Summary plots of initial screen data
###################################
command = 'Rscript initial_screen_analysis.R'

outputs = ["supplemental_figures/initial_screens.arrayed.markers.pdf", "supplemental_figures/initial_screens.pooled.markers.pdf", "supplemental_figures/initial_screens.matched_deg.pdf"]

pipeline.add(command, name=INITIAL_SCREENS_ANALYSIS, dependencies=[INITIAL_SCREENS_MATCHED_DEG], outputs=outputs, memory='40G')

####################################
# Swap rate figures for FACS
####################################
# Calculate proportions from sequencing data
command = 'python generate_facs_proportions.py --fastq_list data/facs_swap_experiment/*.fastq.gz --output_file facs_gfp_bfp_proportions.txt'
output_file = 'tables/facs_gfp_bfp_proportions.txt'
pipeline.add(command, name=FACS_SWAP_RATE, outputs=[output_file])

facs_calculations = "Rscript facs_experiment_swap_rate_figures.R"
sweep_matrix_plot = 'supplemental_figures/optimal_swap_rate_parameter_sweep_matrix.pdf'
sweep_curve_plot = 'figures/optimal_swap_rate_sweep_matrix.fixed_facs_green_proportion.pdf'

command = 'Rscript facs_experiment_swap_rate_figures.R'

pipeline.add(command, name=FACS_SWAP_RATE_FIGURE, dependencies=[FACS_SWAP_RATE], memory='2G', outputs=[sweep_matrix_plot, sweep_curve_plot])

####################################
# Swap rate simulations on crop-seq
####################################
swap_rate_simulation_dir = "temp_data/swap_rate_simulations"
swap_rate_simulation_no_enrichment_dir = "temp_data/swap_rate_simulations.no_enrichment" # this is an extra sim done in review
easygrid.mkdir(swap_rate_simulation_dir)
easygrid.mkdir(swap_rate_simulation_no_enrichment_dir)

for swap_rate in np.arange(0.0, 1.05, 0.05):
    for seed in range(1, 11): 
        for cds_file, condition in zip([cds, cds_no_enrichment, cds, cds_no_enrichment, adamson_tf_pilot_cds, adamson_upr_cds], ['dox_100nm', 'dox_100nm', 'mock', 'mock', 'adamson.tf_pilot', 'adamson.upr']):
            if condition == 'adamson.upr':
                output_file = os.path.join(swap_rate_simulation_dir, 'swap%s_seed%s.%s.5000cells.txt' % (swap_rate, seed, condition))
            elif cds_file == cds_no_enrichment:
                output_file = os.path.join(swap_rate_simulation_no_enrichment_dir, 'swap%s_seed%s.%s.txt' % (swap_rate, seed, condition))
            else:
                output_file = os.path.join(swap_rate_simulation_dir, 'swap%s_seed%s.%s.txt' % (swap_rate, seed, condition))

            command = 'Rscript swapped_assignment_deg.R %s %s --swap_rate %s --treatment %s --seed %s' % (cds_file, output_file, swap_rate, condition, seed)

            # these are large so I already estimated dispersions and size factors so don't repeat
            if condition == 'adamson.upr' or condition == 'adamson.tf_pilot':
                command = command + ' --use_as_is'

            pipeline.add(command, name=SIMULATE_SWAP_RATES, dependencies=[PREPROCESS_CDS], outputs=[output_file], memory='50G')


# Summarize those swap rate simulations with a figure
swap_rate_simulation_figure = 'figures/swap_rate_simulation.pdf'
pipeline.add('Rscript swapped_assignment_figure.R', name=SWAP_RATE_SIM_FIGURE, dependencies=[SIMULATE_SWAP_RATES], outputs=[swap_rate_simulation_figure], memory='5G')


####################################
# TSNE highlighting TP53 and associated markers
####################################
command = 'Rscript main_text_tsne.R'
outputs = ['figures/tsne_by_assignment.png', 'figures/tsne_by_assignment_tp53_markers.png']
pipeline.add(command, name=MAIN_TEXT_TSNE, dependencies=[PREPROCESS_CDS], outputs=outputs, memory='50G')

####################################
# Barcode enrichment
####################################
easygrid.mkdir('temp_data/barcode_enrichment')

# Regular barcode enrichment
command_template = 'Rscript jose_run_barcode_enrichment.R --cds %s --target_level_chisq %s --guide_level_chisq %s --processed_cds %s --cds_mock %s --cds_dox %s --processed_cds_pdata %s'
mock_cds_enrichment = 'temp_data/barcode_enrichment/mock_cds.rds'
dox_cds_enrichment = 'temp_data/barcode_enrichment/dox_100nm_cds.rds'
cds_enrichment = 'temp_data/barcode_enrichment/aggregated_cds.barcode_enrichment.rds'
cds_enrichment_pdata = 'temp_data/barcode_enrichment/aggregated_cds.barcode_enrichment.pdata.txt'
guide_level_chisq = 'temp_data/barcode_enrichment/initial.guide.level.chisq.qval.rds'
target_level_chisq = 'temp_data/barcode_enrichment/initial.target.level.chisq.qval.rds'

command = command_template % (cds, target_level_chisq, guide_level_chisq, cds_enrichment, mock_cds_enrichment, dox_cds_enrichment, cds_enrichment_pdata)
pipeline.add(command, name=BARCODE_ENRICHMENT, dependencies=[PREPROCESS_CDS], outputs=[mock_cds_enrichment, dox_cds_enrichment, cds_enrichment, cds_enrichment_pdata, guide_level_chisq, target_level_chisq], memory='50G')

# Barcode enrichment with a simulated swap rate
mock_cds_enrichment_swapped = 'temp_data/barcode_enrichment/mock_cds.swapped_50.rds'
dox_cds_enrichment_swapped = 'temp_data/barcode_enrichment/dox_100nm_cds.swapped_50.rds'
cds_enrichment_swapped = 'temp_data/barcode_enrichment/aggregated_cds.barcode_enrichment.swapped_50.rds'
cds_enrichment_pdata_swapped = 'temp_data/barcode_enrichment/aggregated_cds.barcode_enrichment.pdata.swapped_50.txt'
guide_level_chisq_swapped = 'temp_data/barcode_enrichment/initial.guide.level.chisq.qval.swapped_50.rds'
target_level_chisq_swapped = 'temp_data/barcode_enrichment/initial.target.level.chisq.qval.swapped_50.rds'

command = command_template % (cds, target_level_chisq_swapped, guide_level_chisq_swapped, cds_enrichment_swapped, mock_cds_enrichment_swapped, dox_cds_enrichment_swapped, cds_enrichment_pdata_swapped)
command = command + ' --swap_rate 0.50'
pipeline.add(command, name=BARCODE_ENRICHMENT, dependencies=[PREPROCESS_CDS], outputs=[mock_cds_enrichment_swapped, dox_cds_enrichment_swapped, cds_enrichment_swapped, cds_enrichment_pdata_swapped, guide_level_chisq_swapped, target_level_chisq_swapped], memory='50G')

####################################
# Barcode enrichment summary figures
####################################
command = 'Rscript jose_target_level_peturbation_figures.R'
pipeline.add(command, name=SUMMARIZE_BARCODE_ENRICHMENT, dependencies=[BARCODE_ENRICHMENT], memory='10G')

command = 'Rscript jose_run_barcode_enrichment.R' 
pipeline.add(command, name=SUMMARIZE_BARCODE_ENRICHMENT, dependencies=[BARCODE_ENRICHMENT], memory='10G')

####################################
# Pairwise DEG
####################################
command_file = 'pairwise_deg_cropseq.jobs'
output_directory = 'temp_data/pairwise_deg/raw_data/singles/cropseq'
easygrid.mkdir(output_directory)

# Raw data
command = 'python all_single_target_degs_scatter.py %s %s --output_directory %s --cds %s --versus_gene NONTARGETING --separator_column treatment' % (pdata, command_file, output_directory, cds)
pipeline.add(command, name=GENERATE_DEG_JOBS, dependencies=[PREPROCESS_CDS], outputs=[command_file], memory='5G')

# Informative
# Mock
output_directory = 'temp_data/pairwise_deg/informative/singles/cropseq'
easygrid.mkdir(output_directory)

#mock_informative_cds = 'temp_data/barcode_enrichment/mock.aggregated_cds.barcode_enrichment.rds'
#mock_informative_pdata = 'temp_data/barcode_enrichment/mock.aggregated_cds_pdata.barcode_enrichment.txt'
mock_informative_cds = 'temp_data/barcode_enrichment/aggregated_cds.barcode_enrichment.rds'
mock_informative_pdata = 'temp_data/barcode_enrichment/aggregated_cds.barcode_enrichment.pdata.txt'

mock_informative_command_file = 'pairwise_deg_cropseq.mock_informative.jobs'
command = 'python all_single_target_degs_scatter.py %s %s --output_directory %s --cds %s --versus_gene NONTARGETING --separator_column treatment --informative --combo_limit 1' % (mock_informative_pdata, mock_informative_command_file, output_directory, mock_informative_cds)
pipeline.add(command, name=GENERATE_DEG_JOBS, dependencies=[PREPROCESS_CDS], outputs=[mock_informative_command_file], memory='5G')

# Dox
#dox_informative_cds = 'temp_data/barcode_enrichment/dox_100nm.aggregated_cds.barcode_enrichment.rds'
#dox_informative_pdata = 'temp_data/barcode_enrichment/dox_100nm.aggregated_cds_pdata.barcode_enrichment.txt'
dox_informative_cds = mock_informative_cds
dox_informative_pdata = mock_informative_pdata
dox_informative_command_file = 'pairwise_deg_cropseq.dox_100nm_informative.jobs'
command = 'python all_single_target_degs_scatter.py %s %s --output_directory %s --cds %s --versus_gene NONTARGETING --separator_column treatment --informative --combo_limit 1' % (dox_informative_pdata, dox_informative_command_file, output_directory, dox_informative_cds)
pipeline.add(command, name=GENERATE_DEG_JOBS, dependencies=[PREPROCESS_CDS], outputs=[dox_informative_command_file], memory='5G')

# Finally these were some analyses requested by reviewers
command = 'Rscript response_to_reviewers.R'
outputs = ['reviewer_figures/cropseq_screen.enrichment.tsne.genotypes.png']
pipeline.add(command, name=REVIEWER_FIGURES, dependencies=[PREPROCESS_CDS], outputs=outputs, memory='50G')

# Run the pipeline here because this generates some commands that we need
success = pipeline.run(dry=False)
del pipeline

if not success:
    sys.exit('First part of pipeline had failures. Not running part2.')

# Resume and run DEG tests
pipeline = easygrid.JobManager('.easygrid2')

deg_commands = []
for command_f in [command_file, mock_informative_command_file, dox_informative_command_file]:
    deg_commands.extend([command.strip() for command in open(command_f)])


for command in deg_commands:
    pipeline.add(command, name=PAIRWISE_DEG, memory='100G')

pipeline.run(dry=True)
