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
import glob


def run(args, seqtype):

    label = 1
    for file in os.listdir(args.variation):
        # find each hacked haploid genome in the directory

        if file.endswith("h1.fa") or file.endswith("h2.fa"):
            if seqtype == 'short':
                sim_short(args.reference, args.variation + file, label, args)
            elif seqtype == 'long':
                sim_long(args.reference, args.variation + file, label, args)
            label += 1

    # Merge bams
    bams = [y for x in os.walk(os.path.abspath(args.output)) for y in glob.glob(os.path.join(x[0], '*.srt.bam'))]
    if len(bams) > 1:
        with open(os.path.abspath(args.output + 'bamstomerge.txt'), 'w') as bamstomerge:

            for bam in bams:
                bamstomerge.write(bam + '\n')

        subprocess.call(['samtools', 'merge', '-@', str(args.threads - 1), '-b', os.path.abspath(args.output + '/bamstomerge.txt'),
             os.path.abspath(args.output + args.prefix + '.tmp.srt.bam')], stderr=open(os.devnull, 'wb'))
        subprocess.call(['samtools', 'index', os.path.abspath(args.output + args.prefix + '.srt.bam')],
                        stderr=open(os.devnull, 'wb'))

    else:
        subprocess.call(['mv', args.output + '*.tmp.srt.bam', args.output + args.prefix + '.srt.bam'], stderr=open(os.devnull, 'wb'))
        subprocess.call(['samtools', 'index', os.path.abspath(args.output + args.prefix + '.srt.bam')], stderr=open(os.devnull, 'wb'))

def sim_short(ref_genome_fa, alt_genome_fa, label, args):

    for line in open(args.config, 'r'):
        tmp = line.strip().split("\t")
        chrom = tmp[0]
        start = int(tmp[1])
        end = int(tmp[2])
        vaf = float(tmp[3])
        short_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, start, end, vaf, args)

    # combine all chroms of short reads
    merge_fq(args.output, 'short')

    os.remove(os.path.abspath(args.output + '*.tmp1.fq'))
    os.remove(os.path.abspath(args.output + '*.tmp2.fq'))

    # Do the alignment
    with open(os.path.abspath(args.output + 'tmp.sam'), 'w') as samout:
        subprocess.call(['bwa', 'mem', '-t', str(args.threads), ref_genome_fa, os.path.abspath(args.output + 'short.r1.fq'), os.path.abspath(args.output + 'short.r2.fq')], stdout=samout, stderr=open(os.devnull, 'wb'))

    convert_sam(label, args.threads, args.output)


def sim_long(ref_genome_fa, alt_genome_fa, label, args):

    for line in open(args.config, 'r'):
        tmp = line.strip().split("\t")
        chrom = tmp[0]
        start = int(tmp[1])
        end = int(tmp[2])
        vaf = float(tmp[3])
        long_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, start, end, vaf, args)

    # combine all chroms of long reads
    merge_fq(args.output, 'long')
    os.remove(os.path.abspath(args.output + '*.sim_0001.fq'))

    # Do the alignment
    with open(os.path.abspath(args.output + 'tmp.sam'), 'w') as samout:
        subprocess.call(['ngmlr', '-t', str(args.threads), '-r', ref_genome_fa, '-q', args.output + 'lr.fq'], stdout=samout, stderr=open(os.devnull, 'wb'))

    convert_sam(label, args.cores, args.output)

def short_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, start, end, vaf, args):
    """ Create short reads for each chromosome with wgsim """

    # read altered genome of one chromosome
    with open(os.path.abspath(args.output + 'alt_chrom.fa'), 'w') as regionout:
        subprocess.call(['samtools', 'faidx', alt_genome_fa, chrom], stdout=regionout, stderr=open(os.devnull, 'wb'))

    this_alt_chrom_fa = pyfaidx.Fasta(os.path.abspath(args.output + 'alt_chrom.fa'))
    this_chrom_alt_seq = this_alt_chrom_fa[0:len(this_alt_chrom_fa[chrom])].seq
    this_chrom_alt_Ns = this_chrom_alt_seq.count('N')

    # altered genome is smaller than the ref genome
    if len(this_chrom_alt_seq) < end - start:
        num_reads = round((args.cov * (len(this_chrom_alt_seq) - this_chrom_alt_Ns)) / args.read_length) / 2
    else:
        num_reads = round((args.cov * (end - start - this_chrom_alt_Ns)) / args.read_length) / 2

    if vaf != 100:
        # Extract reference genome of this chromosome
        with open(os.path.abspath(args.output + 'reference_chrom.fa'), 'w') as out:
            subprocess.call(['samtools', 'faidx', ref_genome_fa, chrom], stdout=out, stderr=open(os.devnull, 'wb'))

        alt_reads = round((num_reads / 100) * vaf)
        ref_reads = num_reads - alt_reads

        # generate reads for alt genome
        subprocess.call(['wgsim', '-e', str(args.error), '-d', str(args.insertsize), '-s', str(args.standarddev), '-N', str(alt_reads), '-1',
             str(args.read_length), '-2', str(args.read_length), '-R', str(args.indels), '-X', str(args.probability),
             os.path.abspath(args.output + "alt_chrom.fa"), os.path.abspath(args.output + 'alt.1.fq'),
             os.path.abspath(args.output + 'alt.2.fq')], stderr=open(os.devnull, 'wb'),
            stdout=open(os.devnull, 'wb'))

        # generate reads for ref genome
        subprocess.call(['wgsim', '-e', str(args.error), '-d', str(args.insertsize), '-s', str(args.standarddev), '-N', str(ref_reads), '-1',
                         str(args.read_length), '-2', str(args.read_length), '-R', str(args.indels), '-X',
                         str(args.probability),os.path.abspath(args.input + 'reference_chrom.fa'),os.path.abspath(args.output + 'ref.1.fq'),
                         os.path.abspath(args.output + 'ref.1.fq')], stderr=open(os.devnull, 'wb'),
                        stdout=open(os.devnull, 'wb'))


        # combine R1 fastq
        with open(args.output + '{0}.tmp1.fq'.format(chrom), 'w') as out:
            subprocess.call(['cat', os.path.abspath(args.output + 'alt.1.fq'), os.path.abspath(args.output + 'ref.1.fq')], stdout=out, stderr=open(os.devnull, 'wb'))

        # remove tmp fastq files
        os.remove(args.output + 'alt.1.fq')
        os.remove(args.output + 'ref.1.fq')

        # combine R2 fastq
        with open(args.output + '{0}.tmp2.fq'.format(chrom), 'w') as out:
            subprocess.call(['cat', os.path.abspath(args.output + 'alt.2.fq'), os.path.abspath(args.output + 'ref.2.fq')], stdout=out, stderr=open(os.devnull, 'wb'))

        # remove tmp fastq files
        os.remove(args.output + 'alt.2.fq')
        os.remove(args.output + 'ref.2.fq')

        # remove the tmp reference fastq file
        os.remove(args.input + 'reference_chrom.fa')

    else:
        subprocess.call(['wgsim', '-e', str(args.error), '-d', str(args.insertsize), '-s', str(args.standarddev), '-N',
                         str(num_reads), '-1', str(args.read_length), '-2', str(args.read_length), '-R', str(args.indels), '-X',
                         str(args.probability), os.path.abspath(args.output + "alt_chrom.fa"), os.path.abspath(args.output + '{0}.tmp1.fq'.format(chrom)),
                         os.path.abspath(args.output + '{0}.tmp2.fq'.format(chrom))], stderr=open(os.devnull, 'wb'),
                        stdout=open(os.devnull, 'wb'))


def long_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, start, end, vaf, args):

    # read altered genome of one chromosome
    with open(os.path.abspath(args.output + 'alt_chrom.fa'), 'w') as regionout:
        subprocess.call(['samtools', 'faidx', alt_genome_fa, chrom], stdout=regionout, stderr=open(os.devnull, 'wb'))

    if vaf != 100:
        # Extract reference genome of this chromosome
        with open(os.path.abspath(args.output + 'reference_chrom.fa'), 'w') as out:
            subprocess.call(['samtools', 'faidx', ref_genome_fa, chrom], stdout=out, stderr=open(os.devnull, 'wb'))

        alt_cov = (args.cov / 100) * vaf
        ref_cov = args.cov - alt_cov

        subprocess.call(['pbsim', '--model_qc', args.model_qc, '--prefix', args.output + 'simref', '--length-mean', str(args.mean),
                         '--accuracy-mean', str(args.accuracy), '--difference-ratio', args.ratio, '--depth', str(ref_cov),
                         os.path.abspath(args.output + 'reference_chrom.fa')], stderr=open(os.devnull, 'wb'),
                        stdout=open(os.devnull, 'wb'))

        # remove tmp files created by pbsim
        os.remove(os.path.abspath(args.output + 'reference_chrom.fa'))
        os.remove(os.path.abspath(args.output + 'simref_0001.ref'))
        os.remove(os.path.abspath(args.output + 'simref_0001.maf'))

        subprocess.call(['pbsim', '--model_qc', args.model_qc, '--prefix', args.output + 'simalt', '--length-mean', str(args.mean),
             '--accuracy-mean', str(args.accuracy), '--difference-ratio', args.ratio, '--depth', str(ref_cov),
             os.path.abspath(args.output + 'alt_chrom.fa')], stderr=open(os.devnull, 'wb'),
            stdout=open(os.devnull, 'wb'))

        # remove tmp files created by pbsim
        os.remove(os.path.abspath(args.output + 'simalt_0001.ref'))
        os.remove(os.path.abspath(args.output + 'simalt_0001.maf'))

        with open(os.path.abspath(args.output + '{0}.sim_0001.fastq'.format(chrom)), 'w') as regionout:

            subprocess.call(['cat', os.path.abspath(args.output + 'simref_0001.fastq'), os.path.abspath(args.output + 'simalt_0001.fastq')],stdout=regionout, stderr=open(os.devnull, 'wb'))

    else:

        subprocess.call(['pbsim', '--model_qc', args.model_qc, '--prefix', args.output + '{0}.sim'.format(chrom), '--length-mean', str(args.mean),
             '--accuracy-mean', str(args.accuracy), '--difference-ratio', args.ratio, '--depth', str(args.cov),
             os.path.abspath(args.output + 'reference_chrom.fa')], stderr=open(os.devnull, 'wb'),
            stdout=open(os.devnull, 'wb'))

        os.remove(os.path.abspath(args.output + 'alt_chrom.fa'))
        os.remove(os.path.abspath(args.output + '{0}.sim_0001.ref'.format(chrom)))
        os.remove(os.path.abspath(args.output + '{0}.sim_0001.maf'.format(chrom)))

def merge_fq(output, seqtype):
    if seqtype == 'short':
        merged_r1_fq = open(os.path.join(output, '{0}.r1.fq').format(seqtype), 'w')
        merged_r2_fq = open(os.path.join(output, '{0}.r2.fq').format(seqtype), 'w')
        for file in os.listdir(output):
            if file.endswith('.tmp1.fq'):
                cmd = 'cat %s >> %s'%(os.path.join(output, file), merged_r1_fq)
                subprocess.call(cmd, shell=True)
            elif file.endswith('.tmp2.fq'):
                cmd = 'cat %s >> %s' % (os.path.join(output, file), merged_r2_fq)
                subprocess.call(cmd, shell=True)
    elif seqtype == 'long':
        merged_fq = open(os.path.join(output, '{0}.fq').format(seqtype), 'w')
        for file in os.listdir(output):
            if file.endswith('.sim_0001.fq'):
                cmd = 'cat %s >> %s' % (os.path.join(output, file), merged_fq)
                subprocess.call(cmd, shell=True)


def convert_sam(label, threads, output):

    with open(os.path.abspath(output + 'tmp.bam'), 'w') as bamout:
        subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + 'tmp.sam')], stdout=bamout, stderr=open(os.devnull, 'wb'))

    os.remove(os.path.abspath(output + 'tmp.sam'))

    with open(os.path.abspath(output + label + '.srt.bam'), 'w') as srtbamout:
        subprocess.call(['samtools', 'sort', '-@', str(threads - 1), os.path.abspath(output + 'tmp.bam')], stdout=srtbamout, stderr=open(os.devnull, 'wb'))

    os.remove(os.path.abspath(output + 'tmp.bam'))

    subprocess.call(['samtools', 'index', os.path.abspath(output + label + '.tmp.srt.bam')], stderr=open(os.devnull, 'wb'))
