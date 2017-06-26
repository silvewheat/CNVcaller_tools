# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 21:32:33 2017

@author: YudongCai
"""

import os
import numpy as np
import pandas as pd
from collections import defaultdict


def loadcnvrfile(cnvfile, chromorder=None):
    """
    parse raw cnvcaller output format
    the first 8 columns are [chr, start, end, number, gap, repeat, gc, kmer]
    the last two columns are [average, sd]
    the in-between columns are each individual
    when chromorder was state, chromsome will be sorted as the given order, contings that not states will be remove, you can use this option to exclude some contigs
    """
    if chromorder:
        chromorder = loadonecol(chromorder)
    df = pd.read_csv(cnvfile, sep='\t', low_memory=False)
    df['chr'] = df['chr'].astype('str')  # chr is always str type
    df['chr'] = df['chr'].astype('category', categories=chromorder, ordered=True)
    if chromorder:
        print('original contigs number is %s' % df.shape[0])
        df.dropna(inplace=True)  # contig that not listed in chromorder will be removed
        df = df.sort_values('chr')
        print('after sort and filter, %s remains' % df.shape[0])
    return df


def loadgtfile(gtfile, chromorder=None):
    """
    parse HYH genotype file
    file not contain original cnv data
    the first 8 columns are [chr, start, end, number, gap, repeat, gc, kmer]
    the last two columns are [aa, Aa, AA, AB, BB, BC, M, average1, average2, average3, sd1, sd2, sd3]
    the in-between columns are each individual's genotype
    when chromorder was state, chromsome will be sorted as the given order, contings that not states will be remove, you can use this option to exclude some contigs
    """
    if chromorder:
        chromorder = loadonecol(chromorder)
    df = pd.read_csv(gtfile, sep='\t', low_memory=False)
    df['chr'] = df['chr'].astype('str')  # chr is always str type
    df['chr'] = df['chr'].astype('category', categories=chromorder, ordered=True)
    if chromorder:
        print('original contigs number is %s' % df.shape[0])
        df.dropna(inplace=True)  # contig that not listed in chromorder will be removed
        df = df.sort_values('chr')
        print('after sort and filter, %s remains' % df.shape[0])
    return df


def loadann(annfile):
    """
    load ann result
    file format [chr, start, end, type, genename]
    """
    df = pd.read_csv(annfile, sep='\t', low_memory=False)
    df['chr'] = df['chr'].astype('str')  # chr is always str type
    df['chr'] = df['chr'].astype('category', categories=None, ordered=True)
    return df



def loadcorrlist(corrlistfile):
    """
    load corr file list to distinguish male and female
    1 is female, 2 is male
    """
    sexdict = {}
    with open(corrlistfile) as f:
        for line in f:
            line = line.strip()
            basename = os.path.basename(line)
            sex = basename.split('_')[-1]
            smid = '_'.join(basename.split('_')[:-6])
            sexdict[smid] = sex
    return sexdict



def loadtwocol(infile):
    """
    load two columns file,
    save to a dict,
    with the first as keys the second as values
    """
    kvdict = {}
    with open(infile) as f:
        for nline, line in enumerate(f):
            line = line.strip()
            if line: # to exclude blank line
                k, v = line.split('\t')
                kvdict[k] = v
    return kvdict


def loadtwocol_dlist(infile):
    """
    load two columns file,
    save to a dict, value is list, the second col as k
    """
    kvdict = defaultdict(list)
    with open(infile) as f:
        for line in f:
            line = line.strip()
            if line:
                k, v = line.split('\t')
                kvdict[v].append(k)
    return kvdict



def loadonecol(infile):
    """
    load one columns file,
    save as a list
    """
    slist = []
    with open(infile) as f:
        for line in f:
            line = line.strip()
            if line: # exclude blank line
                slist.append(line)
    return slist




