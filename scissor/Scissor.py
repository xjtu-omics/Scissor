#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2020/3/2
'''

from scissor.arguments import *
from scissor import SIM, SEQ
import sys
import logging
import os


def main():

    options = parse_arguments(sys.argv[1:])

    # external_tools = ['wgsim', 'pbsim']
    #
    # for tools in external_tools:
    #     if which(tools) is None:
    #         logging.error(tools + ' was not found in the current environment.')

    # validate template reference genome
    try:
        with open(os.path.abspath(options.reference), 'r') as file:
            assert (file.readline().startswith('>'))
    except:
        logging.error('Reference file does not exist, is not readable or is not a valid .fasta file')


    if options.command == 'sim':
        logging.info("Simulate rearrangement")
        logging.info("Template reference genome: {0}".format(options.reference))
        logging.info("Output directory: {0}".format(options.output))

        SIM.run(options)


    elif options.command == 'short':
        logging.info("Short read sequencing with wgsim")
        logging.info("Template reference genome: {0}".format(options.reference))
        logging.info("Variation genome: {0}".format(options.variation))
        logging.info("Output directory: {0}".format(options.output))

        SEQ.run(options, 'short')

    elif options.command == 'long':
        logging.info("Long read sequencing with pbsim")
        logging.info("Template reference genome: {0}".format(options.reference))
        logging.info("Variation genome: {0}".format(options.variation))
        logging.info("Output directory: {0}".format(options.output))

        SEQ.run(options, 'long')

if __name__ == '__main__':
    main()