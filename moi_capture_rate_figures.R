library(dplyr)
library(ggplot2)
library(viridis)
source("paths.R")
source("helper_functions.R")

log_likelihoods = read.delim('temp_data/moi_capture_rate_parameter_sweep/crop_seq_log_likelihoods.txt')
log_likelihoods = subset(log_likelihoods, moi <= 3)


plot_capture_moi_loglikelihood = function(log_likelihoods, MOI_LABEL_OFFSET=1.5, CAPTURE_LABEL_OFFSET=-0.1) {
	maximum = log_likelihoods %>%
		group_by(sample) %>%
		arrange(-log_likelihood) %>%
		slice(1) %>%
		mutate(label=paste0('MOI: ', moi, ', ', 'capture: ', capture_probability)) %>%
		ungroup()

	plot = ggplot(log_likelihoods, aes(moi, capture_probability)) +
		geom_raster(aes(fill=log_likelihood)) +
		geom_contour(colour='white', size=0.35, alpha=0.3, bins=50, aes(z=log_likelihood)) +
		geom_point(data=maximum) +
		geom_segment(data=maximum, aes(yend=capture_probability + CAPTURE_LABEL_OFFSET, xend=moi + MOI_LABEL_OFFSET)) +
		geom_label(data=maximum, aes(label=label, y=capture_probability + CAPTURE_LABEL_OFFSET, x=moi + MOI_LABEL_OFFSET), size=3, fontface='bold') +
		scale_fill_viridis(option='magma', name='log likelihood') +
		xlab('multiplicity of infection') +
		ylab('capture rate') +
		theme_cfg(grid_lines=FALSE)

	return(plot)
}

plot_capture_moi_loglikelihood(log_likelihoods) +
	facet_wrap(~sample, ncol=4) +
	ggsave('supplemental_figures/crop_seq_capture_rate_moi_estimate.png', width=6, height=3)
