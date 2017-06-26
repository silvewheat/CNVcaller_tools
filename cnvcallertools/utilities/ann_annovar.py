# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:33:04 2017

@author: Caiyd
"""

import os
import numpy as np
import pandas as pd


def fetchinterval(cnvfile, outprefix):
    with open('%s.region' % outprefix, 'w') as f_write:
        with open(cnvfile) as f_read:
            header = f_read.readline()
            for line in f_read:
                line = line.strip().split('\t')
                f_write.write('%s\t%s\t%s\t0\t0\n' % tuple(line[:3]))


def ann(annovar, buildver, buildpath, outprefix):
    """
    annotate_variation.pl must be add in $PATH
    """
    cmd = 'perl {annv} --outfile {out} --buildver {build} {out}.region {buildpath}'.format(annv=annovar, out=outprefix, build=buildver, buildpath=buildpath)
    print(cmd)
    os.system(cmd)


def reshapeout(outprefix):
    os.remove('%s.exonic_variant_function' % outprefix)
    os.remove('%s.log' % outprefix)
    os.remove('%s.region' % outprefix)
    os.system(""" awk 'BEGIN{print "chr\tstart\tend\ttype\tgenename"};{print $3"\t"$4"\t"$5"\t"$1"\t"$2}' %s.variant_function > %s.ann """ % (outprefix, outprefix))
    os.remove('%s.variant_function' % outprefix)
