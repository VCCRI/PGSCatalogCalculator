import argparse
import pandas as pd
import numpy as np
from cyvcf2 import VCF, Writer


pd.options.mode.chained_assignment = None
np.set_printoptions(precision=3)

__version__ = '0.0.1'

# def writeVCF():
      # for en,rec in enumerate(vcf_in.fetch()):
          # chrom = rec.chrom
          # try:
            # test = int(chrom)
          # except ValueError:
            # continue 
          # if rec.id is not None:
            # id1 = str(rec.id) + ":" + str(rec.ref)
            # id2 = str(rec.id) + ":" + ''.join(rec.alts)
            # if id1 in inIDs:
                # vcf_out.write(rec)
                # inIDs = list(filter(lambda x: x != id1, inIDs))
            # elif id2 in inIDs:
                # vcf_out.write(rec)
                # inIDs = list(filter(lambda x: x != id2, inIDs))
            # if len(inIDs) == 0:
              # break

def checkRec(inIDs, rec):
      id1 = str(rec.ID) + ":" + str(rec.REF)
      id2 = str(rec.ID) + ":" + ''.join(rec.ALT)
      test = list(set(a) & set(b))
      return(id1, id2)

def removeId(inId, exId):
      inIDs = [inID for inID in inIDs if inID not in exId]
      return(inIDs)

def writeVCF(vcf, inIDs, out):
      vcf_in = VCF(vcf)
      vcf_out = Writer(out, vcf_in)
      # vcf_out = VariantFile(out, 'wb', header=vcf_in.header)
      for rec in vcf_in:
      # for en,rec in enumerate(vcf_in.fetch()):
          chrom = rec.CHROM
          try:
            test = int(chrom)
          except ValueError:
            continue 
          id1 = str(rec.ID) + ":" + str(rec.REF)
          id2 = str(rec.ID) + ":" + ''.join(rec.ALT)
          # recChang = list(set.intersection(*map(set,[[id1, id2], inIDs])))
          # if len(recChang) != 0:
            # vcf_out.write_record(rec)
            # inIDs = [inID for inID in inIDs if inID not in recChang[0]]
          # recChang = list(set(id1) & set(id2) & set(inIDS))
          if id1 in inIDs:
              #modify id
              rec.ID = id1
              vcf_out.write_record(rec)
              # inIDS = [x for x in inIDs if id1 not in inIDs] 
              # inIDs = list(filter(lambda x: x != id1, inIDs))
              inIDs = [inID for inID in inIDs if inID not in id1]
          elif id2 in inIDs:
              rec.ID = id2
              #modify id
              # vcf_out.write(rec)
              vcf_out.write_record(rec)
              # inIDs = list(filter(lambda x: x != id2, inIDs))
              inIDs = [inID for inID in inIDs if inID not in id2]
          if len(inIDs) == 0:
            break

def extract_vcf(vcf, score, out, header):
    '''Stream the vcf, keeping only rows of interest and re-naming variant IDs'''
    if header:
        scorefile = pd.read_csv(score, sep=' ')
        scorefile = scorefile.head(20000)
        if (scorefile.columns[1] in ['A','C','G','T']):
            raise ValueError('Looks like the score file does not have a header? Try removing the --header flag')
    else:
        scorefile = pd.read_csv(score, sep=' ', comment='#', header=None)
        if (not scorefile.iloc[0].tolist()[1] in ['A','C','G','T']):
            raise ValueError('Looks like the score file has a header? Try adding the --header flag')
    #store the variants in dataframe by chromosome
    vardict = {}
    needCols = ['id', 'effect_allele']
    rsCols = ['rsID', 'effect_allele']
    if set(needCols).issubset(list(scorefile)):
      scorefile['chr_name'] = [int(x.split(':')[0]) for x in scorefile[scorefile.columns[0]]]
      # scorefile['id'] = scorefile['chr'].map(str)+':'+scorefile['chr_pos'].map(str)+':'+":"+scorefile['effect_allele'].map(str)
      ##Flaw if reference_allele is not reference_allele
      # scorefile['oid'] = scorefile['chr_name'].map(str)+':'+scorefile['chr_position'].map(str)+':'+scorefile['alternate_allele'].map(str)+":"+scorefile['effect_allele'].map(str)
      import pdb
      for i in range(1,23,1):
          goodvariants = None
          goodvariants = set(scorefile[scorefile['chr_name'].map(int)==i]['id'].tolist())
          vardict[i] = goodvariants
      vcf_in = VCF(vcf)
      vcf_out = Writer(out, vcf_in)
      #use pysam to stream through the VCF
      for rec in VCF(vcf_in):
      # for en,rec in enumerate(vcf_in.fetch()):
          chrom = rec.chrom
          try:
            test = int(chrom)
          except ValueError:
            continue 
          id1 = rec.chrom + ":" + str(rec.pos) + ":" + rec.ref
          id2 = rec.chrom + ":" + str(rec.pos) + ":" + ''.join(rec.alts)

          if id1 in vardict[int(chrom)]:
              rec.ID = id1
              vardict[int(chrom)].remove(id1)
              ## TODO Find a way  to check if set of set is empty
              # vcf_out.write(rec)
              vcf_out.write_record(rec)
          elif id2 in vardict[int(chrom)]:
              rec.ID = id2
              vardict[int(chrom)].remove(id2)
              # vcf_out.write(rec)
              vcf_out.write_record(rec)
          else:
              if len(vardict) == 0:
                break
    elif set(rsCols).issubset(list(scorefile)):
      inIDs = scorefile['rsID'].map(str)+':'+scorefile['effect_allele'].map(str)
      inIDs = inIDs.tolist()
      writeVCF(vcf, inIDs,out)
      # vcf_in = VariantFile(vcf)  # auto-detect input format

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', type=str, help='Path to the vcf.gz to convert. Should be of typical VCF format, in gzipped form and must have tabix-index file in same location. Multi-allelic variants must have separate rows for each variant.', required=True)
parser.add_argument('--score', type=str, help='Score file for score, formatted as input for PLINK. First 3 columns should be Variant, Effect Allele, and Effect Weight, respectively. Assumes the file is tab-delimited. The variant ID should be denoted as chromosome:position:a1:a2, where the order of a1 and a2 does not matter. Any commented lines at the top with \"#\" are ignored.', required=True)
parser.add_argument('--header', help='Flag to indicate that the score file has a header row (ignoring lines beginning with \"#\").', action='store_true')
parser.add_argument('--out', type=str, help='The name of the output bcf file. Default is out.bcf in current directory.', default="out.bcf")
if __name__ == "__main__":
    #Check of the input to make sure everything necessary is provided
    args = parser.parse_args()
    import cProfile
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
