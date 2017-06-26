# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 14:02:22 2017

@author: YudongCai
"""

import click
from cnvcallertools.utilities.parse import *
from cnvcallertools.utilities.groupdescribe import *
from cnvcallertools.utilities.filter import *
from cnvcallertools.utilities.ann_annovar import *
from cnvcallertools.utilities.reshape import *


@click.group()
def main():
    """
    tools for CNVcaller
    """



@main.command('groupdescribe', short_help='describe info by each group')
@click.option('--cnvfile', help='cnvfile generated from CNVcaller')
@click.option('--groupfile', help='groupfile to decribe sample group, the first col is sample id, the second col is group id')
@click.option('--outfile', help='output file name')
@click.option('--chromorder', default=None, help='file to decibre chrom order, default is empty')
def gdescribe(cnvfile, groupfile, outfile, chromorder):
    """
    decribe for each group
    """
    cnvdf = loadcnvrfile(cnvfile, chromorder)
    gdict = loadtwocol_dlist(groupfile)
    newdf = groupdescribe(cnvdf, gdict)
    newdf.to_csv(outfile, sep='\t', index=False, float_format='%.2f', na_rep='NA')



@main.command('filtersample', short_help='filter sub group from an cnvfile')
@click.option('--cnvfile', help='cnvfile generated from CNVcaller')
@click.option('--samplelistfile', help='one sample id per row')
@click.option('--outfile', help='output file name')
@click.option('--chromorder', default=None, help='file to decibre chrom order, default is empty')
def filtersm(cnvfile, samplelistfile, outfile, chromorder):
    """
    filter sample and recalculate avg and sd
    also can be used to sort sample order
    or sort chromsome order(use --chtomorder option)
    """
    cnvdf = loadcnvrfile(cnvfile, chromorder)
    samplelist = loadonecol(samplelistfile)
    newdf = filtersample(cnvdf, samplelist)
    newdf.to_csv(outfile, sep='\t', index=False, float_format='%.2f', na_rep='NA')


@main.command('filtergenotype', short_help='filter genotype result')
@click.option('--gtfile', help='gtfile is same as cnvfile, just replace copy number to genotype')
@click.option('--majorgtcount', help='major genotype count cutoff', type=int)
@click.option('--outfile', help='output file name')
@click.option('--chromorder', default=None, help='file to decibre chrom order, default is empty')
def filtergt(gtfile, majorgtcount, outfile, chromorder):
    """
    cnv with major genotype count larger than majorgtcount will be exclude
    """
    gtdf = loadgtfile(gtfile, chromorder)
    newdf = filtergenotype(gtdf, majorgtcount)
    newdf.to_csv(outfile, sep='\t', index=False, float_format='%.2f', na_rep='NA')


@main.command('annannovar', short_help='use annovar to annotate cnv interval')
@click.option('--cnvfile', help='cnvfile generated from CNVcaller')
@click.option('--outprefix', help='output prefix')
@click.option('--annovar', help='annotate_variation.pl file and path')
@click.option('--buildver', help='{buildver}_refGene.txt, only the part in brackets')
@click.option('--buildpath', help='the path contain the {buildver}_refGene.txt')
def annannovar(cnvfile, outprefix, annovar, buildver, buildpath):
    fetchinterval(cnvfile, outprefix)
    ann(annovar, buildver, buildpath, outprefix)
    reshapeout(outprefix)


@main.command('addanninfo', short_help='add annotaion information to table')
@click.option('--intable', help='input table file, add ann to this')
@click.option('--annfile', help='annfile produced from annannovar')
@click.option('--outtable', help='output table file')
def addanninfo(intable, annfile, outtable):
    """
    add annotaion information produced from annannovar to other table file
    table must contain chr, start, end. e.g. raw cnv file or genotype file
    """
    targetdf = loadcnvrfile(intable, chromorder=None)
    anndf = loadann(annfile)
    newdf = addann(targetdf, anndf)
    newdf.to_csv(outtable, sep='\t', index=False, float_format='%.2f', na_rep='NA')


@main.command('gt2pedmap', short_help='convert genotype file to ped map file')
@click.option('--gtfile', help='gtfile is same as cnvfile, just replace copy number to genotype')
@click.option('--samplelistfile', help='sample included')
@click.option('--chromlistfile', help='chrom that included')
@click.option('--outprefix', help='output prefix')
def gt2pedmap(gtfile, samplelistfile, chromlistfile, outprefix):
    gtdf = loadgtfile(gtfile)
    samplelist = loadonecol(samplelistfile)
    chromlist = loadonecol(chromlistfile)
    gt2pedmapfile(gtdf, outprefix, samplelist, chromlist)


if __name__ == '__main__':
    main()