library(ggplot2)
library(viridis)
source('helper_functions.R')
source('paths.R')

# Observed rates of contamination of blue barcodes into green
# sorted sample and vice versa from FACS adjusted
# for the mixed cell control contamination rate
observed_blue_contamination_in_green = 0.224 - 0.018
observed_green_contamination_in_blue = 0.333 - 0.011

# Simulate what rates of contaminations you would expect for blue into green
# and green into blue barcode contamination at different swap rates
# and fraction of green plasmids (to account for non-50:50 ratio)
swap_rate = seq(0.1, 0.7, 0.01)
fraction_green_plasmid = seq(0.25, 0.75, 0.005)

simulations <<- data.frame(color=c(), swap_rate=c(), expected_value=c())

for(rate in swap_rate) {
	for(green_fraction in fraction_green_plasmid) {
		expected_green_contamination_in_blue = rate * green_fraction
		expected_blue_contamination_in_green = rate * (1 - green_fraction)

		total_error = (expected_blue_contamination_in_green - observed_blue_contamination_in_green) ^ 2 + (expected_green_contamination_in_blue - observed_green_contamination_in_blue) ^ 2

		results = data.frame("swap_rate"=c(rate), "fraction_green"=c(green_fraction), "total_error"=c(total_error))
		simulations <<- rbind(simulations, results)
	}
}

# Get the best result from sim
minimum = subset(simulations, total_error == min(total_error))

# Generate heatmap over search space, annotating the best parameters
LABEL_X_OFFSET = -0.175
LABEL_Y_OFFSET = 0.15
ggplot(simulations, aes(swap_rate, fraction_green)) +
      geom_raster(aes(fill=total_error)) +
      geom_contour(colour='#d3d3d3', size=0.3, alpha=0.25, bins=100, aes(z=total_error)) +
      scale_fill_viridis(option='magma', name="error") +
	  geom_point(data=minimum, color="white") +
      geom_segment(data=minimum, color="white", aes(yend=swap_rate + LABEL_Y_OFFSET, xend=fraction_green + LABEL_X_OFFSET)) +
	  geom_label(data=minimum, size=3, aes(label=paste("Fraction Green: ", fraction_green, "; Swap Rate: ", swap_rate, sep=""), y=swap_rate + LABEL_Y_OFFSET, x=fraction_green + LABEL_X_OFFSET)) +
      theme_cfg() +
      xlab('lentivirus swap rate') +
      ylab('fraction green plasmid') +
      scale_x_continuous(breaks=seq(0.0, 0.7, 0.1)) +
      scale_y_continuous(breaks=seq(0.25, 0.75, 0.1)) +
      ggsave(file.path(SUPPLEMENT_FIGURE_DIR, "optimal_swap_rate_parameter_sweep_matrix.pdf"), width=4.5, height=3)


# Also make a plot of the error if we fix the proportion as the value inferred from FACS data
# and then find the best fitting swap rate with green fraction fixed
FRACTION_GREEN_CELLS_FROM_FACS = 4.59 / (4.59 + 2.85)

## Subset out the simulation with closest fraction of green cells to FACS data
facs_measured_ratio_bfp_gfp = subset(simulations, abs(fraction_green - FRACTION_GREEN_CELLS_FROM_FACS) == min(abs(fraction_green - FRACTION_GREEN_CELLS_FROM_FACS)))

## Get the best solution
optimal_swap_rate = subset(facs_measured_ratio_bfp_gfp, total_error == min(total_error))

## Generate a plot of the error with most likely swap rate annotated
LABEL_X_OFFSET = 0
LABEL_Y_OFFSET = 0.05

ggplot(facs_measured_ratio_bfp_gfp, aes(swap_rate, total_error)) +
	geom_point(color="#53585F") +
	geom_line(color="#53585F") +
	geom_segment(data=optimal_swap_rate, size=0.5, linetype="dashed", color="black", aes(yend=total_error + LABEL_Y_OFFSET, xend=swap_rate + LABEL_X_OFFSET)) +
	geom_label(data=optimal_swap_rate, size=2.9, aes(label=paste("rate: ", swap_rate, sep=""), y=total_error + LABEL_Y_OFFSET, x=swap_rate + LABEL_X_OFFSET)) +
	theme_cfg() +
	xlab('lentivirus swap rate') +
	ylab('squared error from expected') +
	ggsave(file.path(FIGURE_DIR, "optimal_swap_rate_sweep_matrix.fixed_facs_green_proportion.pdf"), height=2.5, width=1.75)

