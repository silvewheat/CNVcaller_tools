# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 14:02:22 2017

@author: YudongCai
"""

import click
from cnvcallertools.utilities.parse import *
from cnvcallertools.utilities.groupdescribe import *


@click.group()
def cli():
    """
    tools for CNVcaller
    """



@cli.command('groupdescribe', short_help='describe info by each group')
@click.option('--cnvfile', help='cnvfile generated from CNVcaller')
@click.option('--groupfile', help='groupfile to decribe sample group, the first col is sample id, the second col is group id')
@click.option('--oufile', help='output file name')
@click.option('--chromorder', default=None, help='file to decibre chrom order, default is empty')
def groupdescribe(cnvfile, groupfile, outfile, chromorder):
    cnvdf = loadcnvrfile(cnvfile, chromorder)
    gdict = loadtwocol_dlist(groupfile)
    newdf = groupdescribe(cnvdf, gdict)
    newdf.to_csv(outfile, sep='\t', index=False)

