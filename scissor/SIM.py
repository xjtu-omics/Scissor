#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2020/2/14
'''
import numpy as np
import random
import pyfaidx
from intervaltree import *
import os
import logging
from plot import Plot
import sys
from utils import merge_fasta

# contigs used to simulate rearrangements, read from chrom size file.
ALLOWED_CONTIGS = []

def run(args):

    """ The main function for simulating rearrangements """
    ## Debug usage
    # alts_type = '/mnt/d/Scissor/simulation/wgs/alts_type.txt'
    # ref_genome_fa = '/mnt/d/data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'
    # exclude_file = '/mnt/d/Scissor/simulation/grch38.excludeRegions.bed'
    # out_dir = '/mnt/d/Scissor/simulation/wgs'
    # chrom_size = '/mnt/d/Scissor/simulation/wgs/grch38.sizes.tsv'
    # min_size = 500
    # max_size = 10000


    # Initialize parameters
    alts_type = args.alts
    min_size = args.min_size
    max_size = args.max_size
    exclude_file = args.exclude
    chrom_size = args.chromsize
    ref_genome_fa = args.reference
    out_dir = args.output

    logging.basicConfig(filename=os.path.abspath(out_dir + '/Scissor.log'), filemode='w', level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    print('Initialized .log file ' + os.path.abspath(out_dir + '/Scissor.log\n'))

    try:

        with open(os.path.abspath(ref_genome_fa), 'r') as file:

            assert (file.readline().startswith('>'))  # genome .file starts with '>'

    except:
        logging.error('Reference file does not exist, is not readable or is not a valid .fasta file')
        sys.exit(1)

    alt_info = out_dir + '/gr_info.bed'
    writer = open(alt_info, 'w')
    load_allowed_contigs(chrom_size)
    alts_dict_by_chrom = initial_var_dict()

    alt_type_dict = load_alts_type_file(alts_type)

    # Prepare data for each rearrangement
    events_info = create_random_regions(alts_type, min_size, max_size, exclude_file, chrom_size)

    exclude_region_tree_dict = load_exclude_regions(exclude_file)

    logging.info("Data prepared for {0} rearrangements".format(len(events_info)))

    immutable_ref = pyfaidx.Fasta(os.path.abspath(ref_genome_fa))

    for event, info in events_info.items():
        chrom = info[0]
        this_chrom = immutable_ref[chrom]
        chrom_sequence = this_chrom[0:len(this_chrom)].seq
        ref, alt = alt_type_dict[info[3]]
        ref_tokens = ref.split(",")
        alt_tokens = alt.split(",")
        template = chrom_sequence[int(info[1]):int(info[2])]
        seg_pos_dict = dict(zip(ref_tokens[:-1], info[4].tolist()))

        print("Modify sequence at {0}:{1}-{2} of {3}".format(chrom, info[1], info[2], info[3]))

        logging.info("Modify sequence at {0}:{1}-{2} of {3}".format(chrom, info[1], info[2], info[3]))

        alt_segment_info, alt_sequence, extra_segment_info = create_variant_info(chrom, ref_tokens, alt_tokens, template, int(info[1]), seg_pos_dict, immutable_ref, exclude_region_tree_dict, min_size, max_size)

        # rearranged alt sequence of shared segments between reference and variation genome
        alts_dict_by_chrom[chrom].append((info[1], info[2], alt_sequence))


        # adding alt sequence of variation unique segments, usually a deletion caused by cut-paste.
        if len(extra_segment_info) > 0:
            for extra_seg, extra_info in extra_segment_info.items():
                chrom = extra_info[0]
                alts_dict_by_chrom[chrom].append((extra_info[3], extra_info[4], ''))

        this_alt = "{0}\t{1}\t{2}".format(info[0], info[1], info[2])

        writer.write("{0}\t{1}\t{2}\n".format(this_alt, ref + "," + alt, segment_info_to_string(alt_segment_info)))
        # region_writer = open(out_dir + "{0}_{1}_{2}.bed".format(info[0], info[1], info[2]), 'w')

        # alt_end = info[1] + len(alt_sequence)

        # prepare sequence for dotplot validation
        ref_seq_for_dotplot = chrom_sequence[int(info[1])- 10000:int(info[2]) + 10000]
        alt_seq_for_dotplot = chrom_sequence[int(info[1])- 10000:int(info[1])] + alt_sequence + chrom_sequence[int(info[2]):int(info[2]) + 10000]

        Plot.run("{0}_{1}_{2}".format(info[0], info[1], info[2]), ref_seq_for_dotplot, alt_seq_for_dotplot, args)

        # print("Original sequence length {0}, modified sequence length {1}".format(int(info[2]) - int(info[1]), len(alt_sequence)))

    for chrom in ALLOWED_CONTIGS:
        ref_genome_seq = immutable_ref[chrom]
        alt_sequence = concatenate_sequence(alts_dict_by_chrom[chrom], ref_genome_seq)
        logging.info("Modified {0}: {1}, original: {2}".format(chrom, len(alt_sequence), len(immutable_ref[chrom])))
        write_sequence(out_dir, chrom, alt_sequence)

    merge_fasta(args.output, args.haploid)

def create_random_regions(alts_type, min_size, max_size, exclude_file, chrom_size):
    """ Create random regions for rearrangements """

    excluded_regions = load_exclude_regions(exclude_file)
    chrom_size_dict = {}
    events_info_dict = {}
    for line in open(chrom_size, 'r'):
        chrom_size_dict[line.split("\t")[0]] = int(line.split("\t")[1])

    event_count = 1
    for line in open(alts_type, 'r'):
        tmp = line.strip().split("\t")
        ref_tokens = tmp[1].split(",")[:-1]
        num_events = int(tmp[3])
        for i in range(num_events):
            ref_seg_sizes = create_segment_size(ref_tokens, min_size, max_size)
            # Find a non-overlapping position for this rearrangement
            chrom, start, end = random_region_of_size(sum(ref_seg_sizes), chrom_size_dict, excluded_regions)
            # add this region to the current chrom interval tree
            excluded_regions[chrom][start:end] = (start, end)
            id = "{0}_{1}_{2}_GR{3}".format(chrom, start, end, event_count)
            seg_info = np.zeros([len(ref_tokens), 2]).astype(int)
            seg_info[0][0] = start
            seg_info[0][1] = ref_seg_sizes[0]
            # collect info for each segments
            for k in range(1, len(ref_tokens)):
                seg_info[k][0] = seg_info[k - 1][0] + ref_seg_sizes[k - 1]
                seg_info[k][1] = ref_seg_sizes[k]
            # save each event [chrom, start, end, type, ref_segment_info] to a dictionary for further usage
            events_info_dict[id] = (chrom, start, end, tmp[0], seg_info)

            event_count += 1

    return events_info_dict

def random_region_of_size(size, chromsizes, intervals):
    """ Create a non-overlapping region on the genome of specific size """

    chrom_idx = np.random.randint(0, len(chromsizes))
    chrom = list(chromsizes.keys())[chrom_idx]
    chrom_size = chromsizes[chrom]
    interval_tree_this_chrom = intervals[chrom]

    start_pos = np.random.randint(50000, chrom_size - size)
    overlapped = interval_tree_this_chrom.overlap(start_pos, start_pos + size)

    while overlapped:
        start_pos = np.random.randint(0, chrom_size - size)
        overlapped = interval_tree_this_chrom.overlap(start_pos, start_pos + size)

    return chrom, start_pos, start_pos + size

def create_variant_info(this_chrom, ref_tokens, alt_tokens, sequence, sequence_start, segment_pos_dict, ref_genome, exclude_region_dict, min_size, max_size):
    # Reference position of the altered segments
    alt_extra_segments = {}
    alt_segments_pos = np.zeros(len(alt_tokens)).astype(int)
    exclude_region_tree = exclude_region_dict[this_chrom]
    alt_segments_pos[0] = sequence_start

    # only one segment
    if len(alt_tokens) > 1:
        alt_segments_pos[-1] = sequence_start + len(sequence)

    segments_info = {}
    alt_sequence = ''
    for i in range(0, len(alt_tokens)):
        alt_segment = alt_tokens[i]
        # reverse sequence of original segment
        if "^" in alt_segment:
            tmp_segment = alt_segment.replace("^","")
            if tmp_segment in segment_pos_dict:
                ref_segment_start_pos = segment_pos_dict[tmp_segment][0] - sequence_start
                # update altered segment postion
                if i > 0:
                    alt_segments_pos[i] = alt_segments_pos[i - 1] + segment_pos_dict[tmp_segment][1]
                # add this segment to the alt_sequence
                if ref_tokens[ref_tokens.index(tmp_segment) + 1] == '#':
                    next_ref_segment_start_pos = len(sequence)
                else:
                    next_ref_segment_start_pos = segment_pos_dict[ref_tokens[ref_tokens.index(tmp_segment) + 1]][0] - sequence_start
                altered_segment_sequence = reverse(sequence[ref_segment_start_pos:next_ref_segment_start_pos])
                alt_sequence += altered_segment_sequence
                # save this altered segments
                if alt_segment in segments_info:
                    segments_info[alt_segment].append((this_chrom, alt_segments_pos[i], alt_segments_pos[i] + len(altered_segment_sequence), altered_segment_sequence, 'alt', 'h1'))
                else:
                    segments_info[alt_segment] = [(this_chrom, alt_segments_pos[i], alt_segments_pos[i] + len(altered_segment_sequence), altered_segment_sequence, 'alt', 'h1')]

            # a new segment in alt sequence
            else:
                new_segment_info = create_jump_segment(ref_genome, exclude_region_tree, min_size, max_size, 'rev')
                from_chrom = new_segment_info[0]
                from_chrom_sequence, from_chrom_sequence_start = new_segment_info[2], new_segment_info[3]
                alt_sequence += from_chrom_sequence

                if new_segment_info[1] == 'cut':
                    alt_extra_segments[alt_segment] = new_segment_info

                # update new altered segment position
                if i > 0:
                    alt_segments_pos[i] = alt_segments_pos[i - 1] + len(from_chrom_sequence)

                # update exclude regions
                if from_chrom == this_chrom:
                    exclude_region_tree[from_chrom_sequence_start: from_chrom_sequence_start + len(from_chrom_sequence)] = (from_chrom_sequence_start, from_chrom_sequence_start + len(from_chrom_sequence))
                else:
                    if from_chrom in exclude_region_dict:
                        exclude_region_dict[from_chrom][from_chrom_sequence_start:from_chrom_sequence_start + len(from_chrom_sequence)] = (from_chrom_sequence_start, from_chrom_sequence_start + len(from_chrom_sequence))
                    else:
                        exclude_region_tree[from_chrom_sequence_start: from_chrom_sequence_start + len(from_chrom_sequence)] = (from_chrom_sequence_start, from_chrom_sequence_start + len(from_chrom_sequence))

        # not a reverse segment
        else:
            if alt_segment in segment_pos_dict:
                ref_segment_start_pos = segment_pos_dict[alt_segment][0] - sequence_start
                # update altered segment postion
                if i > 0:
                    alt_segments_pos[i] = alt_segments_pos[i - 1] + segment_pos_dict[alt_segment][1]

                if ref_tokens[ref_tokens.index(alt_segment) + 1] == '#':
                    next_ref_segment_start_pos = len(sequence)
                else:
                    next_ref_segment_start_pos = segment_pos_dict[ref_tokens[ref_tokens.index(alt_segment) + 1]][0] - sequence_start

                altered_segment_sequence = sequence[ref_segment_start_pos:next_ref_segment_start_pos]
                alt_sequence += altered_segment_sequence
                # save this altered segments
                if alt_segment in segments_info:
                    segments_info[alt_segment].append((this_chrom, alt_segments_pos[i], alt_segments_pos[i] + len(altered_segment_sequence), altered_segment_sequence, 'alt', 'h1'))
                else:
                    segments_info[alt_segment] = [(this_chrom, alt_segments_pos[i], alt_segments_pos[i] + len(altered_segment_sequence), altered_segment_sequence, 'alt', 'h1')]
            else:
                new_segment_info = create_jump_segment(ref_genome, exclude_region_tree, min_size, max_size, 'fwd')
                from_chrom = new_segment_info[0]
                from_chrom_sequence, from_chrom_sequence_start = new_segment_info[2], new_segment_info[3]
                if new_segment_info[1] == 'cut':
                    alt_extra_segments[alt_segment] = new_segment_info
                alt_sequence += from_chrom_sequence

                # update new altered segment position
                if i > 0:
                    alt_segments_pos[i] = alt_segments_pos[i - 1] + len(from_chrom_sequence)

                # update exclude regions
                if from_chrom == this_chrom:
                    # add source region to the interval tree
                    exclude_region_tree[from_chrom_sequence_start: from_chrom_sequence_start + len(from_chrom_sequence)] = (from_chrom_sequence_start, from_chrom_sequence_start + len(from_chrom_sequence))
                else:
                    if from_chrom in exclude_region_dict:
                        exclude_region_dict[from_chrom][from_chrom_sequence_start: from_chrom_sequence_start + len(from_chrom_sequence)] = (from_chrom_sequence_start, from_chrom_sequence_start + len(from_chrom_sequence))
                    else:
                        exclude_region_tree[from_chrom_sequence_start: from_chrom_sequence_start + len(from_chrom_sequence)] = (from_chrom_sequence_start, from_chrom_sequence_start + len(from_chrom_sequence))



    # save segment ref info
    for seg in segments_info.keys():
        segments_info[seg].append((this_chrom, segment_pos_dict[seg.replace("^","")][0], -1, '', 'ref'))

    # adding novel segment
    for seg_symbol, info in alt_extra_segments.items():
        alt_pos = alt_segments_pos[alt_tokens.index(seg_symbol)]
        segments_info[seg_symbol] = [(this_chrom, alt_pos, info[3] + len(info[2]), info[2], 'alt'), (info[0], info[3], info[3] + len(info[2]), '', 'ref')]

    return segments_info, alt_sequence, alt_extra_segments

def create_jump_segment(ref_genome, exclude_regions, min_size, max_size, is_reverse):
    """ Simulate a jumping segment, either cut-paste or copy-paste """

    from_chrom = ALLOWED_CONTIGS[np.random.randint(0, len(ALLOWED_CONTIGS))]
    cut_copy = random.choice(["cut", "copy"])
    from_chrom_sequence, from_chrom_sequence_start = get_sequence_from(ref_genome, from_chrom, min_size, max_size, exclude_regions)
    if is_reverse == 'rev':
        from_chrom_sequence = reverse(from_chrom_sequence)
    return (from_chrom, cut_copy, from_chrom_sequence, from_chrom_sequence_start, from_chrom_sequence_start + len(from_chrom_sequence))

def load_exclude_regions(exclude_file):
    """ Load exclude regions of the genome. """

    exclude_region_dict = {}
    for line in open(exclude_file, 'r'):
        tmp = line.strip().split("\t")
        chrom = tmp[0]
        if chrom in exclude_region_dict:
            exclude_region_dict[chrom][int(tmp[1]) : int(tmp[2])] = (int(tmp[1]), int(tmp[2]))
        else:
            exclude_region_dict[chrom] = IntervalTree()
            exclude_region_dict[chrom][int(tmp[1]) : int(tmp[2])] = (int(tmp[1]), int(tmp[2]))
    return exclude_region_dict

def load_allowed_contigs(chrom_size):
    """ Load contig names from the chrom size file """
    for line in open(chrom_size, 'r'):
        tmp = line.strip().split('\t')
        ALLOWED_CONTIGS.append(tmp[0])

def segment_info_to_string(alt_segment_dict):
    """ Convert the segment information of current chromosome to string for output """

    out_str = ""
    for segment, vars_info in alt_segment_dict.items():
        alt_infos = []
        alt_info = vars_info[0]
        ref_info = vars_info[-1]
        out_str += "SEG={0},{1},ALT={2},{3},REF={4},{5};".format(segment, len(alt_info[3]), alt_info[0], alt_info[1], ref_info[0], ref_info[1])

    return out_str[:-1]

def create_segment_size(ref_symbol, min_size, max_size):
    """ Randomly assign size between [min_size, max_size] to a segment """

    segment_num = len(ref_symbol)
    sizes = []
    sum_len = 0
    while segment_num > 0:
        size = random.randrange(min_size, max_size)
        sizes.append(size)
        sum_len += size
        segment_num -= 1

    return sizes

def get_sequence_from(ref_genome, from_chrom, min_size, max_size, restrict_region_tree):
    """ Get sequence from a specific chromosome of size in range [min_size, max_size] """

    from_chrom_ref = ref_genome[from_chrom]
    from_chrom_seq = from_chrom_ref[0:len(from_chrom_ref)].seq

    sequence_start = np.random.randint(10000, len(from_chrom_seq) - 10000)
    sequence_size = np.random.randint(min_size, max_size)

    while restrict_region_tree.overlap(sequence_start, sequence_start + sequence_size):
        sequence_start = np.random.randint(10000, len(from_chrom_seq) - 10000)
        sequence_size = np.random.randint(1000, 50000)

    return from_chrom_seq[sequence_start:sequence_start + sequence_size], sequence_start

def load_alts_type_file(alts_type):
    alt_type_dict = {}
    try:
        for line in open(alts_type, 'r'):
            tmp = line.strip().split("\t")
            ref = tmp[1]
            alt = tmp[2]
            alt_type_dict[tmp[0]] = (ref, alt)
    except:
        logging.error("Cannot find the rearrangement type file!")
        sys.exit(1)

    return alt_type_dict

def initial_var_dict():
    var_dict = {}
    for chrom in ALLOWED_CONTIGS:
        var_dict[chrom] = []
    return var_dict

def reverse(sequence):
    """ Create a reverse complementary sequence """
    trans = str.maketrans('ATGC', 'TACG')
    new_seq = sequence[::-1].translate(trans)

    return new_seq

def concatenate_sequence(var_list, ref_genome_seq):
    """ Build the whole variation genome with rearrangements in var_list """

    if var_list == []:
        return ref_genome_seq[0:len(ref_genome_seq)].seq

    alt_seq = ''
    sorted_vars = sorted(var_list, key=lambda v:v[1])
    alt_seq += ref_genome_seq[0:sorted_vars[0][0]].seq
    alt_seq += sorted_vars[0][2]
    for i in range(1, len(sorted_vars)):
        # add unmodified sequence between two events
        alt_seq += ref_genome_seq[sorted_vars[i - 1][1]:sorted_vars[i][0]].seq
        alt_seq += sorted_vars[i][2]

    alt_seq += ref_genome_seq[sorted_vars[-1][1]:len(ref_genome_seq)].seq
    return alt_seq

def write_sequence(out_dir, chrom, alt_sequence):
    out_fasta = '{0}{1}_alt.fa'.format(out_dir, chrom)
    writer = open(out_fasta, 'w')
    writer.write('>' + chrom + '\n' + alt_sequence)
    writer.close()

# if __name__ == '__main__':
#     run('')


