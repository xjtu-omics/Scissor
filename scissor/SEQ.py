#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2020/2/29
'''

import pyfaidx
import os
import subprocess


def sim_short(ref_genome_fa, alt_genome_fa, args):

    for line in open(args.config, 'r'):
        tmp = line.strip().split("\t")
        chrom = tmp[0]
        start = int(tmp[1])
        end = int(tmp[2])
        vaf = float(tmp[3])
        short_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, start, end, vaf, args)

    # combine all chroms of short reads
    with open(args.output + '/sr.r1.fq', 'w') as out:
        subprocess.call(['cat', os.path.abspath(args.output + '/*.r1.fq')], stdout=out, stderr=open(os.devnull, 'wb'))

    with open(args.output + '/sr.r2.fq', 'w') as out:
        subprocess.call(['cat', os.path.abspath(args.output + '/*.r2.fq')], stdout=out, stderr=open(os.devnull, 'wb'))

    os.remove(args.output + '/*.r1.fq')
    os.remove(args.output + '/*.r2.fq')

def sim_long(ref_genome_fa, alt_genome_fa, args):

    for line in open(args.config, 'r'):
        tmp = line.strip().split("\t")
        chrom = tmp[0]
        start = int(tmp[1])
        end = int(tmp[2])
        vaf = float(tmp[3])
        long_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, start, end, vaf, args)

    # combine all chroms of short reads
    with open(args.output + '/lr.r1.fq', 'w') as out:
        subprocess.call(['cat', os.path.abspath(args.output + '/*.r1.fq')], stdout=out, stderr=open(os.devnull, 'wb'))

    with open(args.output + '/lr.r2.fq', 'w') as out:
        subprocess.call(['cat', os.path.abspath(args.output + '/*.r2.fq')], stdout=out, stderr=open(os.devnull, 'wb'))

    os.remove(args.output + '/*.r1.fq')
    os.remove(args.output + '/*.r2.fq')

def short_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, start, end, vaf, args):
    """ Create short reads for each chromosome with wgsim """

    # read altered genome of one chromosome
    alt_fa = pyfaidx.Fasta(os.path.abspath(args.input + "/" + alt_genome_fa))
    alt_seq = alt_fa[0:len(alt_fa[chrom])].seq
    alt_Ns = alt_seq.count('N')

    # altered genome is smaller than the ref genome
    if len(alt_seq) < end - start:
        num_reads = round((args.cov * (len(alt_seq) - alt_Ns)) / args.read_length) / 2
    else:
        num_reads = round((args.cov * (end - start - alt_Ns)) / args.read_length) / 2

    if vaf != 100:
        # Extract reference genome of this chromosome
        with open(os.path.abspath(args.input + '/reference_chrom.fa'), 'w') as out:
            subprocess.call(['samtools', 'faidx', ref_genome_fa, chrom + ':' + str(start) + '-' + str(end)], stdout=out, stderr=open(os.devnull, 'wb'))

        alt_reads = round((num_reads / 100) * vaf)
        ref_reads = num_reads - alt_reads

        # generate reads for alt genome
        subprocess.call(['wgsim', '-e', str(args.error), '-d', str(args.insertsize), '-s', str(args.standarddev), '-N', str(alt_reads), '-1',
             str(args.read_length), '-2', str(args.read_length), '-R', str(args.indels), '-X', str(args.probability),
             os.path.abspath(args.input + "/" + alt_genome_fa), os.path.abspath(args.output + '/alt.1.fq'),
             os.path.abspath(args.output + '/alt.2.fq')], stderr=open(os.devnull, 'wb'),
            stdout=open(os.devnull, 'wb'))

        # generate reads for ref genome
        subprocess.call(['wgsim', '-e', str(args.error), '-d', str(args.insertsize), '-s', str(args.standarddev), '-N', str(ref_reads), '-1',
                         str(args.read_length), '-2', str(args.read_length), '-R', str(args.indels), '-X',
                         str(args.probability),os.path.abspath(args.input + '/reference_chrom.fa'),os.path.abspath(args.output + '/ref.1.fq'),
                         os.path.abspath(args.output + '/ref.1.fq')], stderr=open(os.devnull, 'wb'),
                        stdout=open(os.devnull, 'wb'))


        # combine R1 fastq
        with open(args.output + '/{0}.r1.fq'.format(chrom), 'w') as out:
            subprocess.call(['cat', os.path.abspath(args.output + '/alt.1.fq'), os.path.abspath(args.output + '/ref.1.fq')], stdout=out, stderr=open(os.devnull, 'wb'))

        # remove tmp fastq files
        os.remove(args.output + '/alt.1.fq')
        os.remove(args.output + '/ref.1.fq')

        # combine R2 fastq
        with open(args.output + '/{0}.r2.fq'.format(chrom), 'w') as out:
            subprocess.call(['cat', os.path.abspath(args.output + '/alt.2.fq'), os.path.abspath(args.output + '/ref.2.fq')], stdout=out, stderr=open(os.devnull, 'wb'))

        # remove tmp fastq files
        os.remove(args.output + '/alt.2.fq')
        os.remove(args.output + '/ref.2.fq')

        # remove the tmp reference fastq file
        os.remove(args.input + '/reference_chrom.fa')

    else:
        subprocess.call(['wgsim', '-e', str(args.error), '-d', str(args.insertsize), '-s', str(args.standarddev), '-N',
                         str(num_reads), '-1', str(args.read_length), '-2', str(args.read_length), '-R', str(args.indels), '-X',
                         str(args.probability), os.path.abspath(args.input + "/" + alt_genome_fa), os.path.abspath(args.output + '/{0}.r1.fq'.format(chrom)),
                         os.path.abspath(args.output + '/{0}.r2.fq'.format(chrom))], stderr=open(os.devnull, 'wb'),
                        stdout=open(os.devnull, 'wb'))


def long_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, start, end, vaf, args):

    if vaf != 100:
        # Extract reference genome of this chromosome
        with open(os.path.abspath(args.input + 'reference_chrom.fa'), 'w') as out:
            subprocess.call(['samtools', 'faidx', ref_genome_fa, chrom + ':' + str(start) + '-' + str(end)], stdout=out, stderr=open(os.devnull, 'wb'))

        alt_cov = (args.cov / 100) * vaf
        ref_cov = args.cov - alt_cov

        subprocess.call(['pbsim', '--model_qc', args.model_qc, '--prefix', args.output + 'simref', '--length-mean', str(args.mean),
                         '--accuracy-mean', str(args.accuracy), '--difference-ratio', args.ratio, '--depth', str(ref_cov),
                         os.path.abspath(args.input + 'reference_chrom.fa')], stderr=open(os.devnull, 'wb'),
                        stdout=open(os.devnull, 'wb'))

        # remove tmp files created by pbsim
        os.remove(os.path.abspath(args.output + 'reference_chrom.fa'))
        os.remove(os.path.abspath(args.output + 'simref_0001.ref'))
        os.remove(os.path.abspath(args.output + 'simref_0001.maf'))

        subprocess.call(['pbsim', '--model_qc', args.model_qc, '--prefix', args.output + 'simalt', '--length-mean', str(args.mean),
             '--accuracy-mean', str(args.accuracy), '--difference-ratio', args.ratio, '--depth', str(ref_cov),
             os.path.abspath(args.input + alt_genome_fa)], stderr=open(os.devnull, 'wb'),
            stdout=open(os.devnull, 'wb'))

        # remove tmp files created by pbsim
        os.remove(os.path.abspath(args.output + 'simalt_0001.ref'))
        os.remove(os.path.abspath(args.output + 'simalt_0001.maf'))

def run(args):
    if args.seq == 'short':
        sim_short(args.reference, args.variation, args)

    elif args.seq == 'long':
        sim_long(args.reference, args.variation, args)