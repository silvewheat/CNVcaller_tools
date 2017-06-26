# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 15:34:30 2017

@author: Caiyd
"""

import numpy as np
import pandas as pd


def addann(targetdf, anndf):
    """
    anndf is read from func loadann
    targetdf must contain [chr, start, end]
    """
    newdf = pd.merge(targetdf, anndf, on=['chr', 'start', 'end'], how='left')
    assert newdf.shape[0] == targetdf.shape[0]
    return newdf


def gt2pedmapfile(gtdf, outprefix, samplelist, chromlist):
    """
    genotype file to map and ped file
    samplelist to fliter sample
    chromlist to filter contigs
    """
    cols = ['chr', 'start', 'end']
    cols += samplelist
    chromset = set(chromlist)
    gtdf = gtdf.loc[gtdf['chr'].isin(chromset), cols]
    maploc = []
    pedgtlist = []
    with open('%s.ped' % outprefix, 'w') as f_ped:
        for nline, line in enumerate(gtdf[cols].values):
            gts = []
            for gt in line[3:]:
                gt = list(gt) if gt != 'M' else ['M', 'M']
                gts += gt
            if len(set(gts)) <= 2:
                loc = '%s:%s-%s' % tuple(line[:3])
                maploc.append(loc)
                allele1_list = []
                allele2_list = []
                for allele1, allele2 in zip(gts[::2], gts[1::2]):
                    allele1_list.append(allele1)
                    allele2_list.append(allele2)
                pedgtlist.append(allele1_list)
                pedgtlist.append(allele2_list)
        for smid, gts in zip(samplelist, np.array(pedgtlist).T):
            f_ped.write('{FID}\t{IID}\t{PID}\t{MID}\t{Sex}\t{Pheno}\t'.format(
                    FID=smid,
                    IID=smid,
                    PID=0,
                    MID=0,
                    Sex=0,
                    Pheno=0))
            ogt = {}
            ogt[gts[0]] = 'G'
            gts = [ogt.get(x, 'T') for x in gts] # swtich to G and T
            f_ped.write('\t'.join(gts))
            f_ped.write('\n')

    with open('%s.map' % outprefix, 'w') as f_map:
        for loc in maploc:
            f_map.write('{chrom}\t{varidentifier}\t{morganspos}\t{basepairpos}\n'.format(
                    chrom=loc.split(':')[0],
                    varidentifier=loc,
                    morganspos=0,
                    basepairpos=loc.split(':')[1].split('-')[0]))




