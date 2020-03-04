#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2020/3/2
'''

class Sequence():
    def __init__(self, sequence):
        self.bases = sequence
        self.segSequence = []

    def add_sequence(self, sequence):
        sequence = sequence.upper()
        self.bases += sequence

    def clear(self):
        self.bases = ""

    def length(self):
        return len(self.bases)

    def getBases(self):
        return self.bases

    def getReverseComplementBases(self):
        inv_seq = ""
        for i in range(len(self.bases) - 1, -1, -1):
            bp = self.bases[i]
            inv_bp = ''
            if bp == 'A':
                inv_bp = 'T'
            elif bp == 'T':
                inv_bp = 'A'
            elif bp == 'C':
                inv_bp = 'G'
            elif bp == 'G':
                inv_bp = 'C'
            else:
                inv_bp = 'N'

            inv_seq += inv_bp

        return inv_seq