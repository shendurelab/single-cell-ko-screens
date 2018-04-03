from __future__ import print_function
import argparse
import pysam
import collections
import itertools
import pickle


def get_total_reads(bamfile):
    """
    Get the total number of reads from a BAM file.
    """
    read_counts = 0
    index_stats = pysam.idxstats(bamfile)

    # Versions of pysam differ in their output so check to make sure have a list
    if not isinstance(index_stats, list):
        index_stats = index_stats.split('\n')

    for line in index_stats:
        entries = line.strip().split('\t')

        if len(entries) == 4:
            read_counts += int(entries[2]) + int(entries[3])

    return read_counts


def load_whitelist(file_object):
    """
    Load mutation barcode whitelist into a set.

    Args:
            file_object (file): opened file correponding to whitelist file with one barcode per line
    Returns:
            set of string: set of whitelist barcodes in all uppercase

    Raises:
            ValueError if there are any duplicate barcodes in whitelist.
    """

    barcodes = [line.strip().upper() for line in file_object]

    if len(set(barcodes)) != len(barcodes):
        raise ValueError('Provided whitelist contains duplicate barcode entries! All whitelist barcodes must be unique.')

    return set(barcodes)


def hamming_distance(str1, str2, capdistance=None):
    """Count the # of differences between equal length strings str1 and str2. Max diff can be capped using capdistance to short-circuit full calculation when only care about a range of distances."""
    diffs = 0

    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
            if capdistance is not None and diffs >= capdistance:
                return diffs
    return diffs


def correct_barcode(barcode, barcode_whitelist, mismatch_whitelist):
    """
    Correct an observed raw barcode to a barcode in a whitelist or a whitelist of mismatches.

    Args:
            barcode (string): barcode sequence to be corrected
            barcode_whitelist (set of string): set of whitelist barcodes
            mismatch_whitelist (nested dict): dict of mismatch count to dict of mismatched sequences to real sequences

    Returns:
            string: corrected barcodes or None if barcode not correctable.

    """
    barcode = barcode.upper()

    if barcode in barcode_whitelist:
        return barcode

    for mismatch_count in sorted(mismatch_whitelist.keys()):
        if barcode in mismatch_whitelist[mismatch_count]:
            return mismatch_whitelist[mismatch_count][barcode]

    return None


def convert_to_counts(cell_barcode_dict):
    counts = {}
    for cell in cell_barcode_dict:
        counts[cell] = {}

        for barcode in cell_barcode_dict[cell]:
            stats = {'read_count': len(cell_barcode_dict[cell][barcode]),
                     'umi_count': len(set(cell_barcode_dict[cell][barcode]))
                     }
            counts[cell][barcode] = stats
    return counts


def generate_mismatches(sequence, num_mismatches):
    """
    Generate a list of mimatched sequences to a given sequence. Must only contain ATGC.
    This is heavily based on a biostars answer.

    Args:
        sequence (str): The sequence must contain only A, T, G, and C
        num_mismatches (int): number of mismatches to generate sequences for

    Yield:

    """
    letters = 'ACGT'
    sequence = sequence.upper()
    mismatches = []

    for locs in itertools.combinations(range(len(sequence)), num_mismatches):
        sequence_list = [[char] for char in sequence]
        for loc in locs:
            orig_char = sequence[loc]
            sequence_list[loc] = [l for l in letters if l != orig_char]

        for poss in itertools.product(*sequence_list):
            mismatches.append(''.join(poss))

    return mismatches


def construct_mismatch_to_whitelist_map(whitelist, edit_distance):
    """
    Constructs a precomputed set of all mimatches within a specified edit distance and the barcode whitelist.

    Args:
        whitelist (set of str): set of whitelist sequences
        edit_distance (int): max edit distance to consider

    Returns:
        dict: mapping of mismatched sequences to their whitelist sequences
    """

    mismatch_to_whitelist_map = {}

    for mismatch_count in range(1, edit_distance + 1):
        mismatch_to_whitelist_map[mismatch_count] = {}

        for sequence in whitelist:
            sequence = sequence.upper()

            # Generate all possible mismatches in range
            mismatches = generate_mismatches(sequence, num_mismatches=mismatch_count)

            # Construct a mapping to the intended sequences
            for mismatch in mismatches:
                mismatch_to_whitelist_map[mismatch_count][mismatch] = sequence

    return mismatch_to_whitelist_map


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to generate a file with counts of all mutation barcodes found within each cell in a 10X allelic series profiling experiment.')
    parser.add_argument('--input_bams', '-i', nargs='+', help='Position sorted BAM (or list of bams) from 10X pipestance.')
    parser.add_argument('--output_file', '-o', help='Tab delimited file with cell, mutation barcode, read count, umi count. All observed barcodes correctable to a whitelist are reported.')
    parser.add_argument('--whitelist', required=False, default=None, help='Optional mutation barcode whitelist.')
    parser.add_argument('--search_seq', required=True, help='Sequence to search for immediately upstream of the mutation barcode in transcript.')
    parser.add_argument('--barcode_length', required=False, type=int, help='If you are not providing a whitelist, you must provide the length of your barcodes in bp.')
    parser.add_argument('--chimeric_threshold', type=float, default=0.2, help='Threshold for calling a UMI non-chimeric.')
    parser.add_argument('--no_swalign', action='store_true', help='Flag to turn off smith waterman alignment which otherwise requires skbio package / python 3')
    parser.add_argument('--force_correction', type=int, help='Force correction to a specified edit distance. WARNING: if this is greater than or equal to the min edit distance between guides, corrections might reflect cross-talk between guide/barcode sequences.')
    parser.add_argument('--save_raw_counts', help='Save a pickle formatted file of the raw mutation barcode count data. Useful if want to mess around with the raw data rather than processed read and UMI counts.')

    args = parser.parse_args()

    # Try to load alignment library
    try:
        from skbio.alignment import StripedSmithWaterman
    except ImportError:
        if not args.no_swalign:
            raise RuntimeError('skbio library was not found. Please install (requires python 3) or rerun with --no_swalign option.')

    # Try to load progress bar library
    show_progress_bar = False
    try:
        import progressbar
        show_progress_bar = True
    except ImportError:
        print('Install progressbar2 package to show progress bar... proceeding without.')

    # Normalize inputs
    args.search_seq = args.search_seq.upper()

    # Load whitelist
    if args.whitelist:
        barcode_whitelist = load_whitelist(open(args.whitelist))
    else:
        barcode_whitelist = None

    if not barcode_whitelist and args.barcode_length is None:
        # Need to know length somehow...
        raise ValueError('When not providing a whitelist file with expected barcodes with --whitelist, you must specify the length of your barcodes with --barcode_length. Please use one of these two options.')
    elif barcode_whitelist is not None:
        # Look for conflicts and if non found just set automatically based on whitelist
        barcode_lengths = [len(x) for x in barcode_whitelist]

        if len(set(barcode_lengths)) > 1:
            raise ValueError('Whitelist contains sequences of differing length. All barcodes must be same length (can just pad with next base for example).')

        args.barcode_length = barcode_lengths[0]

    # Sanity check on hamming distance within whitelist
    hamming_distances = []
    for barcode in barcode_whitelist:
        for other_barcode in barcode_whitelist:
            if barcode == other_barcode:
                continue
            hamming_distances.append(hamming_distance(barcode, other_barcode))
    min_hamming_distance = min(hamming_distances)

    print('[INFO]: min hamming distance in barcode whitelist: %s' % min_hamming_distance)

    # Determine edit distance for correction
    if args.force_correction is not None:
        if args.force_correction >= min_hamming_distance:
            print('WARNING: edit distance specified by --force_correction (%s) is >= the min hamming distance between guides. Lowest edit distances: %s' % (args.force_correction, sorted(hamming_distances)[1:10]))
        correction_edit_distance = max(int(min_hamming_distance / 2), args.force_correction)
    else:
        correction_edit_distance = int(min_hamming_distance / 2)

    # Precompute sets of potential mismatches within one and two bases
    mismatch_to_whitelist_map = construct_mismatch_to_whitelist_map(barcode_whitelist, correction_edit_distance)

    # Set up alignment if possible
    if not args.no_swalign:
        query = StripedSmithWaterman(args.search_seq)
    else:
        query = None

    search_seq_length = len(args.search_seq)

    # Tracking of UMIs observed for cells for each guide
    mutation_barcodes = {}

    for bam in args.input_bams:

        if show_progress_bar:
            try:
                total_reads = get_total_reads(bam)
                bar = progressbar.ProgressBar(maxval=total_reads, widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ']).start()
            except:
                # if BAM not indexed or some other issue, don't show progress bar
                print('Could not get total reads (BAM may not be indexed), skipping progress bar...')
                show_progress_bar = False
        read_number = 0

        for read in pysam.Samfile(bam):
            # Progress bar
            if show_progress_bar and (read_number % 100000 == 0 or read_number == total_reads - 1):
                bar.update(read_number)

            read_number += 1

            # Ignore mapped reads (can't belong to guide transcripts...)
            if not read.is_unmapped:
                continue

            seq = read.seq.upper()
            tags = dict(read.tags)
            cell = tags.get('CB', None)
            umi = tags.get('UB', None)

            if not cell or not umi:
                continue

            # Given search seq, find barcode in seq keep track of all UMIs for that barcode/cell combo
            search_seq_start = seq.rfind(args.search_seq)
            if search_seq_start >= 0:
                # If there was a perfect match, just use that
                barcode_start = search_seq_start + len(args.search_seq)
            elif not args.no_swalign:
                # Fall back on alignment when no perfect match found
                alignment = query(seq)

                # Skip reads that have a large number of mismatches to search seq
                if (alignment.query_end - alignment.query_begin) < search_seq_length - 3:
                    continue

                barcode_start = alignment.target_end_optimal + (len(args.search_seq) - alignment.query_end)
            else:
                # No match found and no SW backup, so skip read
                continue

            barcode = seq[barcode_start: barcode_start + args.barcode_length]
            barcodes_in_cell = mutation_barcodes.get(cell, dict())

            if barcode_whitelist and barcode not in barcode_whitelist:
                corrected_barcode = correct_barcode(barcode, barcode_whitelist, mismatch_to_whitelist_map)
            else:
                corrected_barcode = barcode

            if not corrected_barcode:
                # Store the barcode with five bases flanking on either side for context
                corrected_barcode = 'unprocessed_%s' % seq[barcode_start - 5: barcode_start + args.barcode_length + 5]

            if corrected_barcode not in barcodes_in_cell:
                barcodes_in_cell[corrected_barcode] = []

            barcodes_in_cell[corrected_barcode].append(umi)
            mutation_barcodes[cell] = barcodes_in_cell

    # Write output file
    if args.save_raw_counts:
        pickle.dump(mutation_barcodes, file=open(args.save_raw_counts, 'wb'))

    with open(args.output_file, 'w') as output_file:
        output_file.write('\t'.join(['cell', 'barcode', 'read_count', 'umi_count']) + '\n')

        for cell in mutation_barcodes:

            # Get a global count of each UMI across all observed barcodes
            all_umis = []
            for barcode in mutation_barcodes[cell]:
                all_umis.extend(mutation_barcodes[cell][barcode])

            umi_counts_all_barcodes = collections.Counter(all_umis)

            for barcode in mutation_barcodes[cell]:

                # Get a count of UMIs within a particular barcode
                umi_counts_barcode = collections.Counter(mutation_barcodes[cell][barcode])

                # Get a list of UMIs that fall below the threshold for being likely non-chimeric
                chimeric_umis = set()

                for umi in umi_counts_barcode:
                    if float(umi_counts_barcode[umi]) / umi_counts_all_barcodes[umi] < args.chimeric_threshold:
                        chimeric_umis.add(umi)

                # Filter any chimeric UMIs before outputing final stats
                filtered_umis = [umi for umi in mutation_barcodes[cell][barcode] if umi not in chimeric_umis]

                read_count = len(filtered_umis)
                umi_count = len(set(filtered_umis))

                output_file.write('\t'.join([cell, barcode, str(read_count), str(umi_count)]) + '\n')
