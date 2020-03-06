#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2020/3/2
'''

from plot.Sequence import Sequence
from plot.HashAligner import HashAligner
from plot.PlotSigleImg import PlotSingleImg
import logging


def run(var_id, ref, seq, args):
    logging.info("Generate Dotplot for rearrangement {0}".format(var_id))

    k = args.kmer
    min_match = args.min_match
    repeat_thresh = args.repeat
    mismatch = args.mismatch

    logging.info("Running parameters: K-mer {0}, window size {1}, repeat_thresh {2}".format(k, min_match, repeat_thresh))

    ref = Sequence(ref)
    read = Sequence(seq)

    aligner_ref = HashAligner(k, min_match, mismatch, repeat_thresh)
    aligner_ref.run(ref, ref)

    diff_segs = []
    if args.denoise:
        diff_segs = aligner_ref.get_self_diff_segments()

    avoid_kmers = aligner_ref.get_avoid_kmer()
    y_hashvalues = aligner_ref.get_hash_values()
    aligner_merge = HashAligner(k, min_match, mismatch, repeat_thresh)
    aligner_merge.run(read, ref, diff_segs, y_hashvalues, avoid_kmers)
    segments = aligner_merge.get_merge_segments()

    read_length = read.length()
    ref_length = ref.length()

    # plot imgs
    img_out = args.output + '/{0}.png'.format(var_id)
    PlotSingleImg(segments, read_length, ref_length, img_out, var_id)


# if __name__ == '__main__':
#     ref = 'AATTCTGAAACACTCTTTCAGAGGGTCTGCAAGTGGATATTTTAGAGCTTTGGGACAATTGTGGAAAAGTAAATATCTTCACATAGAAACTACACGGAAGCATTCTGAGAAACTTCTTTGGAGGTGTGCATTCAACTCACAGAGTTGAACCTATCTTTTCATTGAGCAGTTTTGAATCTCTCTTTTTGTAGACTCTGCTTGCAGATATTTGGAGAGCTTTGAGGCCTATTGTGGAAAAGGGAATATGTTCACATAAAAACACACAGAAGCACTCTGAGAAACTTCTTTGTGAGGTGTGCATTCAACTCACAGAGTTGAACCTATCTTTTGATGGAGAAGTCTTGAATCTCTCTTTTTGTAGAAGCTGCATGTGGATATTTGGAGACGTTTGTGGCCTATGGTAGAAAAGGATATATCTTCAAATAAAAACTAGACAGAAGCATTTTGAGAAAATTCTCTGTGCTGTGTGCATTCATATCACATGGTTGAAACTACCTTTTGATTGAGCAGTTTCGAGTCTCTCTGTTTGTACCATCTGCAATGGATATTTGGAGCCCTTTGTGGTCTGTGGTGGAAAAGGAACTATCCTCAAATAAAAACTACACGGAAGTATTCTGAGAAACTTCTTTGTGATGTGTGCATTTATCTCACAGAGTTGAACCTTTGGTTTGATTGAGCAGTTTTGAGATAATCTTTCCATAGAATCTGGAAGTGAATACTTGGATAACTTTGAGATCTATTTTGGAGAAGGAGATATCTTTATATAAAAACTGCACAGAAGCATTCTGAGAAACATCTTTGTGAGGTGTGCAATGAAGTCACAGAGTTGAAACTATCTTTTGATTCAGCAGTTTTGAGTCTCTCTTTTTGCAGAATCTGCGAGTGGATATCTGGAGAACGTTGAGGCCTACTTGGAAAAGGAAATATCTTCACATAAAAACTACGCAGAAGCATTTTGAGATACTTCTTTGTGAGGTGTGCATTCAACTCACAGAGTTGAACTTATCTTTCCATGGAGCACTTTCATATCTCTTTTTTTGTGGAATCTGCAAGTGGATATTTGGAGCTCTTTGCACCCTGTGGTGGAAAGGGAAATATCTTCATATAAAAACTACAAAGAAGCATTCAGAGAAACTTCTTTGTGATGAATGCATTCCTCACACAGAAGTTGAAGCCTTTCTTTTTATTGAGCAGTATTGAAACGCTCCTTTTGCAGAATCACCAAGTGGATATTTGGAGAGCTTTGGGGCCTGATTTGGAAAATGAAATATCTTCAAAGTAAAACTACACAGAACCATTCTGAGAAACTTCTTTATGATG'
#
#     read = 'AATTCTGAAACACTCTTTCAGAGGGTCTGCAAGTGGATATTTTAGAGCTTTGGGACAATTGTGGAAAAGTAAATATCTTCACATAGAAACTACACGGAAGCATTCTGAGAAACTTCTTTGGAGGTGTGCATTCAACTCACAGAGTTGAACCTATCTTTTCATTGAGCAGTTTTGAATCTCTCTTTTTGTAGACTCTGCTTGCAGATATTTGGAGAGCTTTGAGGCCTATTGTGGAAAAGGGAATATGTTCACATAAAAACACACAGAAGCACTCTGAGAAACTTCTTTGTGAGGTGTGCATTCAACTCACAGAGTTGAACCTATCTTTTGATGGAGAAGTCTTGAATCTCTCTTTTTGTAGAAGCTGCATGTGGATATTTGGAGACGTTTGTGGCCTATGGTAGAAAAGGATATATCTTCAAATAAAAACTAGACAGAAGCATTTTGAGAAAATTCTCTGTGCTGTGTGCATTCATATCATGATATGAATGCACACAGCACAGAGAATTTTCTCAAAATGCTTCTGTCTAGTTTTTATTTGAAGATATATCCTTTTCTACCATAGGCCACAAACGTCTCCAAATATCCACATGCAGCTTCTACAAAAAGAGAGATTCAAGACTTCTCCATCAAAAGATAGGTTCAACTCTGTGAGTTGAATGATATGAATGCACACAGCACAGAGAATTTTCTCAAAATGCTTCTGTCTAGTTTTTATTTGAAGATATATCCTTTTCTACCATAGGCCACAAACGTCTCCAAATATCCACATGCAGCTTCTACAAAAAGAGAGATTCAAGACTTCTCCATCAAAAGATAGGTTCAACTCTGTGAGTTGAATGATATGAATGCACACAGCACAGAGAATTTTCTCAAAATGCTTCTGTCTAGTTTTTATTTGAAGATATATCCTTTTCTACCATAGGCCACAAACGTCTCCAAATATCCACATGCAGCTTCTACAAAAAGAGAGATTCAAGACTTCTCCATCAAAAGATAGGTTCAACTCTGTGAGTTGAATGATATGAATGCACACAGCACAGAGAATTTTCTCAAAATGCTTCTGTCTAGTTTTTATTTGAAGATATATCCTTTTCTACCATAGGCCACAAACGTCTCCAAATATCCACATGCAGCTTCTACAAAAAGAGAGATTCAAGACTTCTCCATCAAAAGATAGGTTCAACTCTGTGAGTTGAACATGGTTGAAACTACCTTTTGATTGAGCAGTTTCGAGTCTCTCTGTTTGTACCATCTGCAATGGATATTTGGAGCCCTTTGTGGTCTGTGGTGGAAAAGGAACTATCCTCAAATAAAAACTACACGGAAGTATTCTGAGAAACTTCTTTGTGATGTGTGCATTTATCTCACAGAGTTGAACCTTTGGTTTGATTGAGCAGTTTTGAGATAATCTTTCCATAGAATCTGGAAGTGAATACTTGGATAACTTTGAGATCTATTTTGGAGAAGGAGATATCTTTATATAAAAACTGCACAGAAGCATTCTGAGAAACATCTTTGTGAGGTGTGCAATGAAGTCACAGAGTTGAAACTATCTTTTGATTCAGCAGTTTTGAGTCTCTCTTTTTGCAGAATCTGCGAGTGGATATCTGGAGAACGTTGAGGCCTACTTGGAAAAGGAAATATCTTCACATAAAAACTACGCAGAAGCATTTTGAGATACTTCTTTGTGAGGTGTGCATTCAACTCACAGAGTTGAACTTATCTTTCCATGGAGCACTTTCATATCTCTTTTTTTGTGGAATCTGCAAGTGGATATTTGGAGCTCTTTGCACCCTGTGGTGGAAAGGGAAATATCTTCATATAAAAACTACAAAGAAGCATTCAGAGAAACTTCTTTGTGATGAATGCATTCCTCACACAGAAGTTGAAGCCTTTCTTTTTATTGAGCAGTATTGAAACGCTCCTTTTGCAGAATCACCAAGTGGATATTTGGAGAGCTTTGGGGCCTGATTTGGAAAATGAAATATCTTCAAAGTAAAACTACACAGAACCATTCTGAGAAACTTCTTTATGATG'
#
#     run('test', ref, read)