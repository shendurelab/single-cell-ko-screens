import scipy.stats
import math
import numpy as np
import collections
import argparse

# Adapted from code provided by Altray Dixit
def get_log_likelihood(observed_counts, capture_probability, moi, moi_max=10):
    	#initialize possion distribution with current guess    
	pdf=scipy.stats.poisson.pmf(range(moi_max), moi)
    
	#Zero truncation and renormalization
	pdf[0]=0
	pdf=np.divide(pdf, np.sum(pdf))

	#get probabilities after convolving with binomial distribution
	zibpdf=np.zeros((moi_max,1))
	for k in range(moi_max):
		pf=0
		for j in np.arange(k,moi_max):
			pf+=pdf[j]*scipy.stats.binom.pmf(k,j,capture_probability)
		zibpdf[k]=pf

	#evaluate log likelihood after multiplying with observed values
	ll=1.0
	for k in observed_counts:
		if k >= moi_max:
			continue # discount cells with more than max moi being considered
			
		ll+=observed_counts[k]*np.log(zibpdf[k])

	return float(ll)


def load_guide_counts(file_object, sample_column=None):

	guide_counts = {}

	colnames = next(file_object).strip().split('\t')

	for line in file_object:
		entries = dict(zip(colnames, line.strip().split('\t')))

		# Get guide count and fill NAs with zero
		guide_count = entries['guide_count']

		if guide_count == 'NA':
			guide_count = 0
		else:
			guide_count = int(guide_count)


		# Get sample
		if sample_column and sample_column not in entries:
			raise ValueError('Column specified as sample column not present in pData table' % sample_column)

		if sample_column:
			sample = entries[sample_column]
		else:
			sample = 'NA'

		# Track counts in dict
		if sample not in guide_counts:
			guide_counts[sample] = collections.Counter()

		guide_counts[sample].update([guide_count])

	return guide_counts



if __name__ == '__main__':
	parser = argparse.ArgumentParser('Script to generate a file with log-likelihoods for a set of capture rates and mois for each sample in pData file.')
	parser.add_argument('pdata', help='pData file must include a guide_count column (observed guide count).')
	parser.add_argument('output_file', help='Log format file with the log likelihood for each combination of MOI, capture rate, and sample.')
	parser.add_argument('--sample_column', help='Specify a column from pData that identifies sets of cells for which the moi-capture rate estimates should be estimated separately. Different treatments or different samples, for example.')
	parser.add_argument('--capture_probability_min', help='Min capture probability to consider.', default=0, type=float)
	parser.add_argument('--capture_probability_max', help='Max capture probability to consider.', default=1, type=float)
	parser.add_argument('--capture_probability_step', help='Step size for capture probability sweeping.', default=0.01, type=float)
	parser.add_argument('--moi_min', help='Min MOI to consider.', default=0.05, type=float)
	parser.add_argument('--moi_max', help='Max MOI to consider.', default=10, type=int)
	parser.add_argument('--moi_step', help='Step size for MOI sweeping.', default=0.1, type=float)
	args = parser.parse_args()

	guide_counts = load_guide_counts(open(args.pdata), sample_column=args.sample_column)

	with open(args.output_file, 'w') as output_file:
		# Write header
		output_file.write('\t'.join(['sample', 'capture_probability', 'moi', 'log_likelihood']) + '\n')

		for sample in guide_counts:

			observed_counts = guide_counts[sample]

			for moi in list(np.arange(args.moi_min, args.moi_max, args.moi_step)):
				for capture_probability in list(np.arange(args.capture_probability_min, args.capture_probability_max, args.capture_probability_step)):
					log_likelihood = get_log_likelihood(observed_counts, capture_probability, moi, args.moi_max)

					output_file.write('\t'.join([sample, str(capture_probability), str(moi), str(log_likelihood)]) + '\n')

