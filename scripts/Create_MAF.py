#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_MAF.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create MAF file for OncoKB annotations
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# Usage:
"""
python3 scripts/Create_MAF.py -i {input.Mutations} -o {output.MAF}
"""
#
# TODO:
# 1) 
#
# History:
#  20-08-2024: File creation, write code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
import argparse
import pandas as pd
#-------------------------------------------------------------------------------
# 1.1 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 Create_MAF.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Create MAF file for OncoKB annotations  ')
    parser.add_argument('-i', help='path to selected mutations file',
                        dest='input',
                        type=str)
    parser.add_argument('-s', help='subtype',
                        dest='subtype',
                        type=str)
    parser.add_argument('-o', help='path to output file',
                        dest='output',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

"""
args.input = 'output/KappaHyperExome/Selected_mutations_NSCLC.txt'
"""
#-------------------------------------------------------------------------------
# 2.1 Read data
#-------------------------------------------------------------------------------
# read Mutations
Mutations = pd.read_csv(args.input, sep = "\t")

#-------------------------------------------------------------------------------
# 3.1 Select and translate column names for MAF
#-------------------------------------------------------------------------------
# Select mutations which are present
Mutations = Mutations[Mutations.Mutation_present]
Mutations['NCBI_Build'] = 'GRCh37'
Mutations['Variant_Classification'] = Mutations['exonic.func'].fillna(Mutations.func)
Mutations['Chromosome'] = [chr.replace('chr','') for chr in Mutations.chr]
Mutations['SampleID'] = [sample+'-'+region for sample,region in zip(Mutations.SampleID,Mutations.Region)]
# Rename columns to MAF format
Mutations = Mutations.rename(columns={"SampleID" : "Tumor_Sample_Barcode", "AAChange" : "HGVSp_Short" ,'AAChange' : 'HGVSP','start' : 'Start_Position','stop':'End_Position','NucleotideChange':'HGVSg', 'ref':'Reference_Allele','var' : 'Tumor_Seq_Allele1','var':'Tumor_Seq_Allele2','VAF':'Variant_Allele_Frequency'})

Mutations['HGVSp_Short'] = Mutations['HGVSP']
Mutations['HGVSp'] = Mutations['HGVSP']
Mutations['Tumor_Seq_Allele1'] = Mutations['Tumor_Seq_Allele2']

#-------------------------------------------------------------------------------
# 3.2 Recode Variant classifcation column
#-------------------------------------------------------------------------------
# Define translation dict (see https://github.com/mskcc/vcf2maf/blob/main/vcf2maf.pl)
Translation_dict = {'intronic' : 'Intron', 'nonsynonymous SNV' : 'Missense', 'UTR3' : "3'UTR", 'stopgain SNV' :  'Nonsense_Mutation', 'ncRNA_intronic' :'Intron','synonymous SNV' : 'Silent','ncRNA_exonic' : 'RNA','intergenic' : 'Intergenic', 'upstream;downstream' : '','upstream' : "5'Flank",'frameshift substitution' : 'Frame_Shift_DEl','unknown' : '','ncRNA_UTR3' : "3'UTR",'stoploss SNV' : 'Nonstop_Mutation', 'nonsynonymous' : 'Missense','UTR5' : "5'UTR",'splicing' : "Splice_site", 'downstream' : "3'Flank",'frameshift insertion' : 'Frame_Shift_Ins','immediate-stopgain' : 'Nonsense_Mutation','nonframeshift substitution' : "In_Frame_Del",'ncRNA_UTR5' : "5'UTR", 'nonframeshift insertion' : 'In_Frame_Ins','UTR5;UTR3' : '', 'UNKNOWN' : '','synonymous' : 'Silent', 'ncRNA_splicing' : "Splice_site"}

# Recode Variant Classification
Mutations['Variant_Classification'] = [Translation_dict[x] for x in Mutations['Variant_Classification'] ]

# Select columns
MAF = Mutations.loc[:,['NCBI_Build','Hugo_Symbol','Variant_Classification','Tumor_Sample_Barcode','HGVSp_Short','HGVSp','HGVSg','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','Variant_Allele_Frequency']]

MAF['Cancer_Type'] = args.subtype
# drop duplicates
MAF.drop_duplicates(inplace=True)

#-------------------------------------------------------------------------------
# 3.1 Write to file
#-------------------------------------------------------------------------------
MAF.to_csv(args.output, sep = '\t', index = False)
