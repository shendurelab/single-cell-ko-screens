import argparse
import gzip

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Script to generate count of GFP and BFP barcodes from FACS experiments.')
	parser.add_argument('--fastq_list', nargs='+', help='List of FASTQ files.')
	parser.add_argument('--output_file', help='Output file with counts for GFP and BFP barcodes.')
	args = parser.parse_args()

	# Sequences that match reverse complement of barcode and preceeding sequence to be stringent
	gfp_search_seq = 'GTCTTAAAGGCTCGA' + 'TAAGGGAGACAAGGC'
	bfp_search_seq = 'GTCTTAAAGGCTCGA' + 'AGAGCGGGAAGTGAA'

	with open(args.output_file, 'w') as output_file:
		output_file.write('\t'.join(['fastq_file', 'bfp_count', 'bfp_proportion', 'gfp_count', 'gfp_proportion']) + '\n')

		for fastq in args.fastq_list:
			bfp_counts = 0.0
			gfp_counts = 0.0

			for line in gzip.open(fastq):
				if bfp_search_seq in line:
					bfp_counts += 1
				elif gfp_search_seq in line:
					gfp_counts += 1

			bfp_proportion = bfp_counts / (bfp_counts + gfp_counts)
			gfp_proportion = gfp_counts / (bfp_counts + gfp_counts)

			output_file.write('\t'.join([fastq, str(bfp_counts), str(bfp_proportion), str(gfp_counts), str(gfp_proportion)]) + '\n')
