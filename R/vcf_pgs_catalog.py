# Author: SK, Originally written by Mark Chaffin

"""
This script prepares a vcf.gz file for scoring via PLINK. Takes as input a vcf.gz and score file and harmonizes the two. The vcf.gz should
follow standard vcf format, and the score file should be a typical PLINK score file with first 3 columns Variant, Effect Allele, and Effect Weight.
Variants in the score file are expected to be denoted as chromosome:position:a1:a2. Order of alleles a1 and a2 do not matter in the score file as 
both orders will be checked for in the vcf.gz. Output is saved as a .bcf file for quicker computation downstream with plink.

Usage:
python vcf_prep_for_score_V2.py --vcf=VCF --score=SCORE --header --out=OUT

Ex)
python vcf_prep_for_score_V2.py --vcf /path/to/yourvcf.vcf.gz \
--score CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt \
--header \
--out /path/to/output.bcf

VCF: The path and name of the vcf.gz file to prepare
SCORE: The score file in traditional PLINK format. The first 3 columns should be variant, effect allele, and effect weight. Any rows beginning with #
       in the header of the file will be ignored. If a header is present, use the --header flag, otherwise i 
OUT: The path and name of the output file of interest (to be saved as a .bcf)

"""
from __future__ import division
import os
import re
import sys
import gzip
import time, sys, traceback, argparse
import pandas as pd
import numpy as np
from pysam import VariantFile


pd.options.mode.chained_assignment = None
np.set_printoptions(precision=3)

__version__ = '0.0.1'

def writeVCF():
      for en,rec in enumerate(vcf_in.fetch()):
          chrom = rec.chrom
          try:
            test = int(chrom)
          except ValueError:
            continue 
          if rec.id is not None:
            id1 = str(rec.id) + ":" + str(rec.ref)
            id2 = str(rec.id) + ":" + ''.join(rec.alts)
            if id1 in inIDs:
                vcf_out.write(rec)
                inIDs = list(filter(lambda x: x != id1, inIDs))
            elif id2 in inIDs:
                vcf_out.write(rec)
                inIDs = list(filter(lambda x: x != id2, inIDs))
            if len(inIDs) == 0:
              break


def extract_vcf(vcf, score, out, header):
    '''Stream the vcf, keeping only rows of interest and re-naming variant IDs'''
    import datetime 
    # print(str(datetime.datetime.now()))
    if header:
        scorefile = pd.read_csv(score, sep=' ')
        if (scorefile.columns[1] in ['A','C','G','T']):
            raise ValueError('Looks like the score file does not have a header? Try removing the --header flag')
    else:
        scorefile = pd.read_csv(score, sep=' ', comment='#', header=None)
        if (not scorefile.iloc[0].tolist()[1] in ['A','C','G','T']):
            raise ValueError('Looks like the score file has a header? Try adding the --header flag')
    #store the variants in dataframe by chromosome
    needCols = ['id', 'effect_allele']
    rsCols = ['rsID', 'effect_allele']
    if set(rsCols).issubset(list(scorefile)):
      # inIDs = scorefile['rsID'].map(str)+':'+scorefile['effect_allele'].map(str)
      inIDs = scorefile['rsID']
      inIDs = inIDs.tolist()
      # inIDs = dict((k,None) for k in inIDs)
      inIDs = dict.fromkeys(inIDs)
      vcf_in = VariantFile(vcf)  # auto-detect input format
      vcf_out = VariantFile(out, 'wb', header=vcf_in.header)
      testDat = 0 
      for en,rec in enumerate(vcf_in.fetch()):
          chrom = rec.chrom
          try:
            test = int(chrom)
          except ValueError:
            continue
          if rec.id is not None:
              splitIds = rec.id.split(':')
          else:
              continue
          idList = []
          for i in splitIds:
              idList.append(i+":" + str(rec.ref))
              idList.append(i+":" + ''.join(rec.alts))
          # id1 = str(rec.id) + ":" + str(rec.ref)
          # id2 = str(rec.id) + ":" + ''.join(rec.alts)
          for i in idList:
              if i in inIDs:
                  rec.id = i
                  vcf_out.write(rec)
                  del inIDs[i]

          # if id1 in inIDs:
              # rec.id = id1
              # vcf_out.write(rec)
              # # inIDS = [x for x in inIDs if id1 not in inIDs] 
              # # inIDs = list(filter(lambda x: x != id1, inIDs))
              # # inIDs = [inID for inID in inIDs if inID not in id1]
              # del inIDs[id1]
          # elif id2 in inIDs:
              # rec.id = id2
              # vcf_out.write(rec)
              # # inIDs = list(filter(lambda x: x != id2, inIDs))
              # del inIDs[id2]
              # inIDs = [inID for inID in inIDs if inID not in id2]
          if len(inIDs) == 0:
            break
    elif set(needCols).issubset(list(scorefile)):
      inChrom = re.sub(r'.+_', '', re.sub(r'_DBN.vcf.gz', '', vcf))
      scorefile = scorefile.loc[scorefile['chr_name'].map(str) == inChrom]
      # scorefile['chr_name'] = [int(x.split(':')[0]) for x in scorefile[scorefile.columns[0]]]
      scorefile['id'] = scorefile['chr_name'].map(str)+':'+scorefile['chr_position'].map(str)+':'+scorefile['effect_allele'].map(str)
      ##Flaw if reference_allele is not reference_allele
      # scorefile['oid'] = scorefile['chr_name'].map(str)+':'+scorefile['chr_position'].map(str)+':'+scorefile['alternate_allele'].map(str)+":"+scorefile['effect_allele'].map(str)
      inIDs = scorefile['id'].tolist()
      inIDs = dict.fromkeys(inIDs)
      # inIDs = dict((k,None) for k in inIDs)
      # for i in range(1,23,1):
          # goodvariants = None
          # goodvariants = set(scorefile[scorefile['chr_name'].map(int)==i]['id'].tolist())
          # vardict[i] = goodvariants
      #use pysam to stream through the VCF
      vcf_in = VariantFile(vcf)  # auto-detect input format
      vcf_out = VariantFile(out, 'wb', header=vcf_in.header)
      for en,rec in enumerate(vcf_in.fetch()):
          chrom = rec.chrom

          # try:
            # test = int(chrom)
          # except ValueError:
            # continue 
          id1 = rec.chrom + ":" + str(rec.pos) + ":" + rec.ref
          id2 = rec.chrom + ":" + str(rec.pos) + ":" + ''.join(rec.alts)

          if id1 in inIDs:
              rec.id = id1
              # vardict[int(chrom)].remove(id1)
              del inIDs[id1]
              ## TODO Find a way  to check if set of set is empty
              vcf_out.write(rec)
          elif id2 in inIDs:
          # elif id2 in vardict[int(chrom)]:
              rec.id = id2
              # vardict[int(chrom)].remove(id2)
              del inIDs[id2]
              vcf_out.write(rec)
          else:
              if len(inIDs) == 0:
              # if len(vardict) == 0:
                break


parser = argparse.ArgumentParser()
parser.add_argument('--vcf', type=str, help='Path to the vcf.gz to convert. Should be of typical VCF format, in gzipped form and must have tabix-index file in same location. Multi-allelic variants must have separate rows for each variant.', required=True)
parser.add_argument('--score', type=str, help='Score file for score, formatted as input for PLINK. First 3 columns should be Variant, Effect Allele, and Effect Weight, respectively. Assumes the file is tab-delimited. The variant ID should be denoted as chromosome:position:a1:a2, where the order of a1 and a2 does not matter. Any commented lines at the top with \"#\" are ignored.', required=True)
parser.add_argument('--header', help='Flag to indicate that the score file has a header row (ignoring lines beginning with \"#\").', action='store_true')
parser.add_argument('--out', type=str, help='The name of the output bcf file. Default is out.bcf in current directory.', default="out.bcf")

if __name__ == "__main__":
    #Check of the input to make sure everything necessary is provided
    args = parser.parse_args()
    if args.vcf is None:
        raise ValueError('--vcf is required')
    if args.score is None:
        raise ValueError('--score is required')
    if args.out is None:
        raise ValueError('--out is required')

    #Check if your variant is present in the two data sources
    print('######################################################\n' + \
          '#      RUNNING POLYGENIC SCORE PREPARATION FILE      #\n' + \
          '#    Note: This is a beta script, use at own risk    #\n' + \
          '######################################################')
    extract_vcf(args.vcf, args.score, args.out, args.header)


