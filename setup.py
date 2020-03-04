#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2020/3/4
'''

from setuptools import setup,find_packages

setup(name='Scissor',
  version=1.0,
  description='Simualte complex genome rearrangements',
  url='',
  requires=['python (>= 3.6)'],
  author='Jiadong Lin, Songbo Wang',
  author_email='jiadong324@gmail.com',
  license='LICENSE.txt',
  install_requires=['pyfaidx >= 0.5.5.2', 'intervaltree == 3.0.2', 'matplotlib == 3.1.1', 'numpy >= 1.15.1'],
  zip_safe=False,
  packages=find_packages(),
  include_package_data=True,
  entry_points={'console_scripts': ['VISOR=VISOR.VISOR:main']}
)
