#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2020/2/29
'''

import argparse
from argparse import RawTextHelpFormatter

def parse_arguments():


    parser = argparse.ArgumentParser(prog='Scissor', description='''Complex genome rearrangement simulator.''',
                                     epilog='Author: Jiadong Lin, Songbo Wang \nContact: jiadong324@gmail.com',
                                     formatter_class=RawTextHelpFormatter)

    subparsers = parser.add_subparsers(title='Modules', dest='command', metavar='sim, short, long')

    parser_sim = subparsers.add_parser('sim', help='Simulate genome in FASTA format containing input rearrangement type.')

    required = parser_sim.add_argument_group('Required arguments')

    required.add_argument('-g', dest='reference', help='Template reference genome', metavar='FASTA', required=True)
    required.add_argument('-s', dest='chromsize', help='Dimension of the reference genome', metavar='FASTA', required=True)
    required.add_argument('-t', dest='alts',help='Rearrangement type represented by symbols (tab separated file)', metavar='TXT',  required=True)
    required.add_argument('-x', dest='exclude', help='Regions to exclude for rearrangement scissor', metavar='BED', required=True)
    required.add_argument('-o', dest='output', help='Output folder', metavar='FOLDER', required=True)

    optional = parser_sim.add_argument_group('Optional parameters')
    optional.add_argument('-i', dest='haploid', help='Which halpotype is used for CGRs [h1, h2, default=h1]', type=str, default='h1')
    optional.add_argument('-l', dest='min_size', help='Minimum size of sequence segment to modify', type=int, default=500)
    optional.add_argument('-u', dest='max_size',help='Maximum size of sequence segment to modify', type=int, default=10000)

    plot = parser_sim.add_argument_group('Dotplot parameters')
    plot.add_argument('-k', dest='kmer', help='Size of kmer [default=32]', default=32, type=int)
    plot.add_argument('-w', dest='min_match', help='Minimum length of matched segment through kmer [defaut=50]', default=50, type=int)
    plot.add_argument('-r', dest='repeat', help='Threshold for repeat kmer comparison [defualt=10]', default=10, type=int)
    plot.add_argument('-m', dest='mismatch', help='Mismatch number allowed for each segment [defualt=0]', default=0, type=int)
    plot.add_argument('-d', dest='denoise', help='Plot without repeats from reference sequence [defualt=False]', default=False, action='store_true')

    # short read sequencing
    parser_short = subparsers.add_parser('short', help='Short read sequencing of variation genome.')
    required = parser_short.add_argument_group('Required arguments')

    required.add_argument('-g', dest='reference', help='Template reference genome', metavar='FASTA', required=True)
    required.add_argument('-v', dest='variation', help='Simulated variation genome', metavar='FOLDER', required=True)
    required.add_argument('-o', dest='output', help='Output folder for simulated FASTQ', metavar='FOLDER',
                          required=True)
    required.add_argument('-f', dest='config', help='A configure file for sequencing the alteration genome',
                          metavar='BED', required=True)
    required.add_argument('-n', dest='prefix', help='The prefix of aligned file name', metavar='STRING', required=True)

    wgsim = parser_short.add_argument_group('wgsim parameters')
    wgsim.add_argument('-t', dest='threads', help='Number of cores used for alignment [4]', metavar='', default=4,
                       type=int)
    wgsim.add_argument('-c', dest='coverage', help='Mean coverage for the simulated region [30.0]', metavar='',
                       default=30.0, type=float)
    wgsim.add_argument('-e', dest='error', help='Base error rate [0.010]', metavar='', default=0.010, type=float)
    wgsim.add_argument('-l', dest='read_length', help='Length of reads [150]', metavar='', default=150, type=int)
    wgsim.add_argument('-i', dest='indels', help='Fractions of indels [0.000000001]', metavar='', default=0.000000001,
                       type=float)
    wgsim.add_argument('-p', dest='probability', help='Probability an indel is extended [0.000000001]', metavar='',
                       default=0.000000001, type=float)
    wgsim.add_argument('-is', dest='insertsize', help='0uter distance between the two ends [500]', metavar='',
                       default=500, type=int)
    wgsim.add_argument('-sd', dest='standarddev', help='Standard deviation for insert size [50]', metavar='',
                       default=50, type=int)

    # long read sequencing
    parser_long = subparsers.add_parser('long', help='Long read sequencing of variation genome.')
    required = parser_long.add_argument_group('Required arguments')

    required.add_argument('-g', dest='reference', help='Template reference genome', metavar='FASTA', required=True)
    required.add_argument('-v', dest='variation', help='Simulated variation genome', metavar='FASTA', required=True)
    required.add_argument('-o', dest='output', help='Output folder for simulated FASTQ', metavar='FOLDER',required=True)
    required.add_argument('-f', dest='config', help='A configure file for sequencing the alteration genome',metavar='BED', required=True)
    required.add_argument('-n', dest='prefix', help='The prefix of aligned file name', metavar='STRING', required=True)

    pbsim = parser_long.add_argument_group('pbsim parameters for FASTQ simulations')
    pbsim.add_argument('-t', dest='seq', help='Long read (clr/ccs) sequencing technology [clr]', metavar='', default='clr')
    pbsim.add_argument('-c', dest='coverage', help='Mean coverage for the simulated region [20]', metavar='', default=20.0,type=float)
    pbsim.add_argument('-a', dest='accuracy', help='Mean for simulated reads [0.90]', metavar='', default=0.90,type=float)
    pbsim.add_argument('-l', dest='length', help='Mean length for simulated reads [8000]', metavar='', default=8000,type=int)
    pbsim.add_argument('-r', dest='ratio', help='substitution:insertion:deletion ratio [30:30:40]', metavar='',default='30:30:40', type=str)



    return parser.parse_args()


