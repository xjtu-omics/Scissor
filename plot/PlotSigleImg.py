#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2020/3/2
'''

from matplotlib import pyplot as plt
import os


class PlotSingleImg():

    def __init__(self, segments_orignal, refLength, readLength, outDir, title):

        self.figure_size = 8
        self.ratio = 2
        self.title = title
        self.segments_orignal = segments_orignal

        self.readLength = int(readLength / self.ratio)
        self.refLength = int(refLength / self.ratio)

        self.outDir = outDir
        self.plot()

    def plot(self):

        figsize = self.figure_size, self.figure_size

        fig, ax1 = plt.subplots(1, 1, figsize=figsize)

        for seg in self.segments_orignal:

            if seg.forward():
                ax1.plot([seg.yStart(), seg.yEnd()], [seg.xStart(), seg.xEnd()], color='b', linewidth=2)
            else:
                ax1.plot([seg.yStart(), seg.yEnd()], [seg.xStart(), seg.xEnd()], color='r', linewidth=2)


        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)

        ax1.xaxis.set_ticks_position('top')
        ax1.yaxis.set_ticks_position('left')

        plt.xticks([])
        plt.yticks([])

        font = {
                 'weight': 'normal',
                 'size': 10,
                 }

        ax1.set_ylabel('VARIATION', font)
        ax1.set_xlabel('REFERENCE', font)
        ax1.set_title(self.title, font)
        ax1.xaxis.set_ticks_position('top')
        ax1.invert_yaxis()

        # plt.show()
        plt.savefig(os.path.join(self.outDir))