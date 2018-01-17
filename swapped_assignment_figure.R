library(dplyr)
library(ggplot2)
library(reshape2)
source('helper_functions.R')
source("paths.R")

plot_simulation = function(deg_simulation_files, x_step=0.1, label_size=3.2, label_pattern="%1.1f-fold reduction") {
	degs_by_swap_rate = do.call(rbind, lapply(Sys.glob(deg_simulation_files), read.delim))
	deg_counts = degs_by_swap_rate %>% dplyr::filter(qval <= 0.05) %>% dplyr::group_by(swap_rate, seed) %>% dplyr::summarize(n = n())
	deg_counts$group_color = ifelse(ceiling((deg_counts$swap_rate*100)%%10) == 0, "GROUP1", "GROUP2")

	original_average = mean(subset(deg_counts, swap_rate == 0)$n)
	swap_rate_50_average = mean(subset(deg_counts, swap_rate == 0.5)$n)
	fold_reduction = original_average / swap_rate_50_average

	set.seed(0)
	ggplot(deg_counts, aes(swap_rate, n)) +
		annotate("rect", xmin = 0.475, xmax = 0.525, ymin = -Inf, ymax = Inf, alpha = .2, fill="#ED4D3C") +
		geom_linerange(aes(group=as.factor(swap_rate), ymax = n), ymin = -Inf, alpha=0.1, color="#d3d3d3", linetype="dotted") +
		geom_vline(xintercept=0.475, color="black", linetype="dashed") +
		geom_vline(xintercept=0.525, color="black", linetype="dashed") +
		geom_segment(x=0.525, xend=0.9, y=swap_rate_50_average, yend=swap_rate_50_average, linetype="dashed", color="#d3d3d3") +
		geom_label(label=sprintf(label_pattern, fold_reduction), size=label_size, x=0.9, y=swap_rate_50_average) +
		geom_jitter(aes(color=group_color), size=0.85) +
		geom_smooth(color="#ED4D3C", se=FALSE) +
		scale_x_continuous(breaks = seq(0, 1, by = x_step), limits=c(-0.05, 1.05)) +
		scale_color_manual(values = c("GROUP1"="#53585F", "GROUP2"="#d3d3d3")) +
		guides(color=FALSE) +
		theme_cfg(grid_lines=FALSE) +
		ylim(0, max(deg_counts$n) + 10) +
		xlab('simulated swap rate') +
		ylab('differentially\nexpressed genes')
}

deg_simulation_files.cropseq = Sys.glob("temp_data/swap_rate_simulations/*dox_100nm*")

deg_simulation_files.adamson_upr = Sys.glob("temp_data/swap_rate_simulations/*adamson.upr.txt")

# TODO this is actually 6000 cells so these files are named incorrectly, not a big deal...
deg_simulation_files.adamson_upr.downsampled = Sys.glob("temp_data/swap_rate_simulations/*adamson.upr.5000cells.txt")
deg_simulation_files.adamson_tf_pilot = Sys.glob("temp_data/swap_rate_simulations/*adamson.tf_pilot*")
deg_simulation_files.datlinger = Sys.glob("temp_data/swap_rate_simulations/*stimulated*")

plot_simulation(deg_simulation_files.cropseq) +
	ggsave('supplemental_figures/swap_rate_simulation.pdf', height=2.0, width=4.5)

plot_simulation(deg_simulation_files.adamson_upr) +
	ggsave('supplemental_figures/swap_rate_simulation.adamson_upr.pdf', height=2.0, width=4.5)

plot_simulation(deg_simulation_files.adamson_upr.downsampled) +
	ggsave('supplemental_figures/swap_rate_simulation.adamson_upr.downsampled.pdf', height=2.0, width=4.5)

plot_simulation(deg_simulation_files.adamson_tf_pilot, x_step=0.2, label_size=2.4, label_pattern="%1.1f-fold lower") +
	ylab('diff exprs genes') +
	ggsave('figures/swap_rate_simulation.adamson_tf_pilot.pdf', height=2.75, width=2.75)
