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
import logging
import sys
from shutil import which

def run(args, seqtype):

    logging.basicConfig(filename=os.path.abspath(args.output + '/Scissor_SEQ.log'), filemode='w', level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    print('Initialized .log file ' + os.path.abspath(args.output + '/Scissor_SEQ.log'))

    external_tools = ['wgsim', 'ngmlr', 'pbsim', 'bwa', 'samtools']

    for tool in external_tools:
        if which(tool) is None:
            logging.error(tool + ' was not found in the current environment.')
    try:
        with open(os.path.abspath(args.config), 'r') as config:
            assert (config.readline().split("\t"))
    except:
        logging.error("Scissor config file has to be TAB separated")
        sys.exit(1)

    logging.info("Template genome: {0}".format(os.path.abspath(args.reference)))
    logging.info("Variation genome: {0}".format(os.path.abspath(args.variation)))
    logging.info("Output directory: {0}".format(os.path.abspath(args.output)))

    label = 1
    for file in os.listdir(args.variation):
        # find each hacked haploid genome in the directory
        if file.endswith("h1.fa") or file.endswith("h2.fa"):
            fai_index(file)
            if seqtype == 'short':
                logging.info("Short read sequencing with wgsim")
                short_read_sequencing(args.reference, args.variation + file, label, args)

            elif seqtype == 'long':
                logging.info("Long read sequencing with pbsim")

                long_read_sequencing(args.reference, args.variation + file, label, args)

            label += 1

    # Merge bams
    bams = [y for x in os.walk(os.path.abspath(args.output)) for y in glob.glob(os.path.join(x[0], '*.srt.bam'))]
    if len(bams) > 1:
        rm_cmd = ['rm']
        with open(os.path.abspath(args.output + 'bamstomerge.txt'), 'w') as bamstomerge:

            for bam in bams:
                bamstomerge.write(bam + '\n')
                rm_cmd.append(bam)

        subprocess.call(['samtools', 'merge', '-@', str(args.threads - 1), '-b', os.path.abspath(args.output + '/bamstomerge.txt'),
             os.path.abspath(args.output + args.prefix + '.srt.bam')], stderr=open(os.devnull, 'wb'))
        subprocess.call(['samtools', 'index', os.path.abspath(args.output + args.prefix + '.srt.bam')],
                        stderr=open(os.devnull, 'wb'))

        rm_cmd.append(os.path.join(args.output, 'bamstomerge.txt'))
        subprocess.call(rm_cmd)

    else:
        subprocess.call(['mv', args.output + '1.srt.bam', args.output + args.prefix + '.srt.bam'], stderr=open(os.devnull, 'wb'))
        subprocess.call(['samtools', 'index', os.path.abspath(args.output + args.prefix + '.srt.bam')], stderr=open(os.devnull, 'wb'))

def short_read_sequencing(ref_genome_fa, alt_genome_fa, label, args):

    for line in open(args.config, 'r'):
        tmp = line.strip().split("\t")
        chrom = tmp[0]
        start = int(tmp[1])
        end = int(tmp[2])
        vaf = float(tmp[3])
        short_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, start, end, vaf, args)

    # combine all chroms of short reads
    merge_fq(args.output, 'short')

    # Do the alignment
    logging.info("Start align short.r1.fq and short.r2.fq to reference {0}".format(os.path.abspath(ref_genome_fa)))

    if not os.path.exists(os.path.abspath(ref_genome_fa + '.sa')):
        try:

            logging.info('Creating bwa index for ' + os.path.abspath(ref_genome_fa))
            bwa_index(os.path.abspath(ref_genome_fa))

        except:
            logging.error('It was not possible to generate bwa index for ' + os.path.abspath(ref_genome_fa))
            sys.exit(1)

    with open(os.path.abspath(args.output + 'tmp.sam'), 'w') as samout:
        subprocess.call(['bwa', 'mem', '-t', str(args.threads), ref_genome_fa, os.path.abspath(args.output + 'short.r1.fq'), os.path.abspath(args.output + 'short.r2.fq')], stdout=samout, stderr=open(os.devnull, 'wb'))

    convert_sam(str(label), args.threads, args.output)

def long_read_sequencing(ref_genome_fa, alt_genome_fa, label, args):

    model_qc = os.path.abspath(os.path.dirname(__file__) + '/model_qc_clr')
    if args.seq == 'ccs':
        model_qc = os.path.abspath(os.path.dirname(__file__) + '/model_qc_ccs')

    logging.info("Sequence reads with model {0}".format(model_qc))

    for line in open(args.config, 'r'):
        tmp = line.strip().split("\t")
        chrom = tmp[0]
        vaf = float(tmp[3])
        long_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, vaf, args, model_qc)

    # combine all chroms of long reads
    merge_fq(args.output, 'long')

    # Do the alignment
    logging.info("Start align {0} to reference {1}".format("long.fq", os.path.abspath(ref_genome_fa)))

    if not os.path.exists(os.path.abspath(ref_genome_fa + '.sa')):
        try:

            logging.info('Creating bwa index for ' + os.path.abspath(ref_genome_fa))
            bwa_index(os.path.abspath(ref_genome_fa))

        except:
            logging.error('It was not possible to generate bwa index for ' + os.path.abspath(ref_genome_fa))
            sys.exit(1)

    with open(os.path.abspath(args.output + 'tmp.sam'), 'w') as samout:
        subprocess.call(['ngmlr', '-t', str(args.threads), '-r', ref_genome_fa, '-q', args.output + 'long.fq'], stdout=samout, stderr=open(os.devnull, 'wb'))

    convert_sam(str(label), args.threads, args.output)

def short_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, start, end, vaf, args):
    """ Create short reads for each chromosome with wgsim """

    # read altered genome of one chromosome
    with open(os.path.abspath(args.output + 'alt_chrom.fa'), 'w') as regionout:
        subprocess.call(['samtools', 'faidx', alt_genome_fa, chrom], stdout=regionout, stderr=open(os.devnull, 'wb'))

    this_alt_chrom_fa = pyfaidx.Fasta(os.path.abspath(args.output + 'alt_chrom.fa'))

    this_alt_chrom = this_alt_chrom_fa[chrom]
    this_chrom_alt_seq = this_alt_chrom[0:len(this_alt_chrom)].seq
    this_chrom_alt_Ns = this_chrom_alt_seq.count('N')

    # altered genome is smaller than the ref genome
    if len(this_chrom_alt_seq) < end - start:
        num_reads = round((args.coverage * (len(this_chrom_alt_seq) - this_chrom_alt_Ns)) / args.read_length) / 2
    else:
        num_reads = round((args.coverage * (end - start - this_chrom_alt_Ns)) / args.read_length) / 2

    if vaf != 100:
        # Extract reference genome of this chromosome
        with open(os.path.abspath(args.output + 'reference_chrom.fa'), 'w') as out:
            subprocess.call(['samtools', 'faidx', ref_genome_fa, chrom], stdout=out, stderr=open(os.devnull, 'wb'))

        with open(os.path.abspath(args.output + 'reference_chrom.fa.fai'), 'w') as idxout:
            subprocess.call(['samtools', 'faidx', args.output + 'reference_chrom.fa'], stdout=idxout, stderr=open(os.devnull, 'wb'))

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
                         str(args.probability),os.path.abspath(args.output + 'reference_chrom.fa'),os.path.abspath(args.output + 'ref.1.fq'),
                         os.path.abspath(args.output + 'ref.2.fq')], stderr=open(os.devnull, 'wb'),
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
        os.remove(args.output + 'reference_chrom.fa')
        os.remove(args.output + 'reference_chrom.fa.fai')
        os.remove(args.output + 'alt_chrom.fa')
        os.remove(args.output + 'alt_chrom.fa.fai')

    else:
        subprocess.call(['wgsim', '-e', str(args.error), '-d', str(args.insertsize), '-s', str(args.standarddev), '-N',
                         str(num_reads), '-1', str(args.read_length), '-2', str(args.read_length), '-R', str(args.indels), '-X',
                         str(args.probability), os.path.abspath(args.output + "alt_chrom.fa"), os.path.abspath(args.output + '{0}.tmp1.fq'.format(chrom)),
                         os.path.abspath(args.output + '{0}.tmp2.fq'.format(chrom))], stderr=open(os.devnull, 'wb'),
                        stdout=open(os.devnull, 'wb'))

        os.remove(args.output + 'alt_chrom.fa')
        os.remove(args.output + 'alt_chrom.fa.fai')

def long_reads_single_chrom(ref_genome_fa, alt_genome_fa, chrom, vaf, args, model_qc):

    # read altered genome of one chromosome
    with open(os.path.abspath(args.output + 'alt_chrom.fa'), 'w') as regionout:
        subprocess.call(['samtools', 'faidx', alt_genome_fa, chrom], stdout=regionout, stderr=open(os.devnull, 'wb'))

    if vaf != 100:
        # Extract reference genome of this chromosome
        with open(os.path.abspath(args.output + 'reference_chrom.fa'), 'w') as out:
            subprocess.call(['samtools', 'faidx', ref_genome_fa, chrom], stdout=out, stderr=open(os.devnull, 'wb'))

        alt_cov = (args.coverage / 100) * vaf
        ref_cov = args.coverage - alt_cov

        subprocess.call(['pbsim', '--model_qc', model_qc, '--prefix', args.output + 'simref', '--length-mean', str(args.mean),
                         '--accuracy-mean', str(args.accuracy), '--difference-ratio', args.ratio, '--depth', str(ref_cov),
                         os.path.abspath(args.output + 'reference_chrom.fa')], stderr=open(os.devnull, 'wb'),
                        stdout=open(os.devnull, 'wb'))

        # remove tmp files created by pbsim
        os.remove(os.path.abspath(args.output + 'reference_chrom.fa'))
        # os.remove(os.path.abspath(args.output + 'reference_chrom.fa.fai'))
        os.remove(os.path.abspath(args.output + 'simref_0001.ref'))
        os.remove(os.path.abspath(args.output + 'simref_0001.maf'))

        subprocess.call(['pbsim', '--model_qc', model_qc, '--prefix', args.output + 'simalt', '--length-mean', str(args.mean),
             '--accuracy-mean', str(args.accuracy), '--difference-ratio', args.ratio, '--depth', str(ref_cov),
             os.path.abspath(args.output + 'alt_chrom.fa')], stderr=open(os.devnull, 'wb'),
            stdout=open(os.devnull, 'wb'))

        with open(os.path.abspath(args.output + '{0}.sim_0001.fastq'.format(chrom)), 'w') as regionout:

            subprocess.call(['cat', os.path.abspath(args.output + 'simref_0001.fastq'), os.path.abspath(args.output + 'simalt_0001.fastq')],stdout=regionout, stderr=open(os.devnull, 'wb'))

        # remove the tmp reference fastq file
        os.remove(os.path.join(args.output, 'alt_chrom.fa'))
        # os.remove(os.path.join(args.output, 'alt_chrom.fa.fai'))
        os.remove(os.path.abspath(args.output + 'simalt_0001.ref'))
        os.remove(os.path.abspath(args.output + 'simalt_0001.maf'))

    else:
        subprocess.call(['pbsim', '--model_qc', model_qc, '--prefix', args.output + '{0}.sim'.format(chrom), '--length-mean', str(args.mean),
             '--accuracy-mean', str(args.accuracy), '--difference-ratio', args.ratio, '--depth', str(args.coverage),
             os.path.abspath(args.output + 'alt_chrom.fa')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

        os.remove(os.path.join(args.output, 'alt_chrom.fa'))
        # os.remove(os.path.join(args.output, 'alt_chrom.fa.fai'))
        os.remove(os.path.abspath(args.output + '{0}.sim_0001.ref'.format(chrom)))
        os.remove(os.path.abspath(args.output + '{0}.sim_0001.maf'.format(chrom)))

def merge_fq(output, seqtype):
    if seqtype == 'short':
        r1_files = []
        r2_files = []
        for file in os.listdir(output):
            if file.endswith('.tmp1.fq'):
                r1_files.append(os.path.join(output, file))
            elif file.endswith('.tmp2.fq'):
                r2_files.append(os.path.join(output, file))

        with open(os.path.join(output, '{0}.r1.fq').format(seqtype), 'w') as merged_r1_fq:
            cmd = ['cat'] + r1_files
            subprocess.call(cmd, stdout=merged_r1_fq)

        with open(os.path.join(output, '{0}.r2.fq').format(seqtype), 'w') as merged_r2_fq:
            cmd = ['cat'] + r2_files
            subprocess.call(cmd, stdout=merged_r2_fq)

        logging.info("All sequenced reads merged to {0},{1}".format(os.path.abspath(os.path.join(output, "short.r1.fq")),
                                                                    os.path.abspath(os.path.join(output, "short.r2.fq"))))

        rm_cmd = ['rm'] + r1_files + r2_files
        subprocess.call(rm_cmd)

    elif seqtype == 'long':
        fq_files = []
        for file in os.listdir(output):
            if file.endswith('.sim_0001.fastq'):
                fq_files.append(os.path.join(output, file))

        with open(os.path.join(output, '{0}.fq').format(seqtype), 'w') as merged_fq:
            cmd = ['cat'] + fq_files
            subprocess.call(cmd, stdout=merged_fq)

        rm_cmd = ['rm'] + fq_files
        subprocess.call(rm_cmd)

        logging.info("All sequenced reads merged to {0}".format(os.path.abspath(os.path.join(output, "long.fq"))))

def bwa_index(fasta):
    subprocess.call(['bwa', 'index', os.path.abspath(fasta)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

def fai_index(fasta):
    subprocess.call(['samtools', 'faidx', os.path.abspath(fasta)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

def convert_sam(label, threads, output):

    with open(os.path.join(output, 'tmp.bam'), 'w') as bamout:
        subprocess.call(['samtools', 'view', '-b', os.path.join(output, 'tmp.sam')], stdout=bamout, stderr=open(os.devnull, 'wb'))

    os.remove(os.path.join(output, 'tmp.sam'))

    with open(os.path.abspath(output + label + '.srt.bam'), 'w') as srtbamout:
        subprocess.call(['samtools', 'sort', '-@', str(threads - 1), os.path.join(output, 'tmp.bam')], stdout=srtbamout, stderr=open(os.devnull, 'wb'))
    os.remove(os.path.abspath(output + 'tmp.bam'))

    logging.info("Alignment sorted")