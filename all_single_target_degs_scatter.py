import subprocess
import argparse
import os
from collections import Counter

def parse_pdata(file_object):
	columns = None
	for line in file_object:

		entries = line.strip().split('\t')

		if columns is None:
			columns = entries
			continue

		yield dict(zip(columns, entries))

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Script to scatter DEG calculations for all "conditions" and "gene" values in pData of CFG CDS on cluster and save to files.')
	parser.add_argument('pdata', help='pData file from CDS')
	parser.add_argument('output_file', help='Output file with the commands.')
	parser.add_argument('--output_directory', help='Output directory for DEG results files.')
	parser.add_argument('--cds', required=True, help='CDS object to use for tests.')
	parser.add_argument('--versus_gene', default="NONTARGETING", help='Gene ID to test against.')
	parser.add_argument('--separator_column', required=True, help='Column in pData indicating which samples from a CDS should be tested as a subgroup.')
	parser.add_argument('--informative', action="store_true", help='Gene ID to test against.')
	parser.add_argument('--combo_limit', type=int, default=1, help='Limit to the order of combos to test (2 = pairwise).')
	parser.add_argument('--combo_number_min', default=3, type=int, help='Limit to the order of combos to test (2 = pairwise).')
	args = parser.parse_args()

	deg_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'get_degs_script.R')

	pdata_entries = [entry for entry in parse_pdata(open(args.pdata))]

	if args.informative:
		combos = Counter([(x[args.separator_column], x['gene']) for x in pdata_entries if x['gene'] != 'NA' and x['informative'] == 'TRUE' and (x['gene'] != args.versus_gene)])
	else:
		combos = Counter([(x[args.separator_column], x['gene']) for x in pdata_entries if x['gene'] != 'NA' and (x['gene'] != args.versus_gene)])

	if not os.path.exists(args.output_directory):
		os.mkdir(args.output_directory)

	with open(args.output_file, 'w') as output_file:
		for condition, gene_id in combos:

			# Exclude combos that occur too few times to warrant testing
			if combos[(condition, gene_id)] <= args.combo_number_min:
				continue

			# Don't test combos that have  the versus gene (NTC) in it
			if args.versus_gene in gene_id and '_' in gene_id:
				continue

			# Don't run anything against itself and don't test combos with an order over the limit
			if gene_id != args.versus_gene and len(gene_id.split('_')) <= args.combo_limit:
				output_file_name = os.path.join(args.output_directory, '%s-%s-%s.deg.txt' % (condition, gene_id, args.versus_gene))

				if args.informative:
					command = 'Rscript %s %s %s %s %s --column_name %s --column_value %s --informative' % (deg_script, args.cds, gene_id, args.versus_gene, output_file_name, args.separator_column, condition)
				else:
					command = 'Rscript %s %s %s %s %s --column_name %s --column_value %s' % (deg_script, args.cds, gene_id, args.versus_gene, output_file_name, args.separator_column, condition)
				output_file.write(command + '\n')
