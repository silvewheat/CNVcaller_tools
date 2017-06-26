# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 22:46:59 2017

@author: Caiyd
"""

def filtersample(cnvdf, smlist):
    """
    filter sample
    and recalculate avg and sd
    """
    allsm = set(cnvdf.columns.tolist())
    losssm = set(smlist) - allsm
    if losssm:
        print('%s sample loss in input cnvfile' % len(losssm))
        print('\n'.join(list(losssm)))
        smlist = [x for x in smlist if x in allsm]
    else:
        print('all sample found')
    cols = ['chr', 'start', 'end', 'number', 'gap', 'repeat', 'gc', 'kmer']
    cols += smlist
    new_df = cnvdf.loc[:, cols]
    new_df['average'] = cnvdf.loc[:, smlist].mean(axis = 1)
    new_df['sd'] = cnvdf.loc[:, smlist].std(axis=1, ddof=1)
    return new_df


def filtergenotype(gtdf, majorgtcount):
    """
    cnv with major genotype count larger than majorgtcount will be exclude
    """
    typelist = ['aa', 'Aa', 'AA', 'AB', 'BB', 'BC', 'M']
    gtdf = gtdf.loc[(gtdf[typelist].max(axis=1) <= majorgtcount), :]
    return gtdf


def filterpar():
    pass

