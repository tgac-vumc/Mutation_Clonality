#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Fetch_FoundationOne_regions.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Obtain FoundationOne panel regions based on techinical supplement:
# https://info.foundationmedicine.com/hubfs/FMI%20Labels/FoundationOne_CDx_Label_Technical_Info.pdf
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) 
#
# History:
#  19-12-2023: File creation, write code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(dplyr)))

#-------------------------------------------------------------------------------
# 1.1 Fetch genes for which all exons are in panel
#-------------------------------------------------------------------------------
# Genes in Table 2. to take exon regions from
ExonGenes <- c('ABL1','BRAF','CDKN1A','EPHA3','FGFR4','IKZF1','MCL1','NKX2-1','PMS2','RNF43','TET2','ACVR1B','BRCA1','CDKN1B','EPHB1','FH','INPP4B','MDM2','NOTCH1','POLD1','ROS1','TGFBR2','AKT1','BRCA2','CDKN2A','EPHB4','FLCN','IRF2','MDM4','NOTCH2','POLE','RPTOR','TIPARP','AKT2','BRD4','CDKN2B','ERBB2','FLT1','IRF4','MED12','NOTCH3','PPARG','SDHA','TNFAIP3','AKT3','BRIP1','CDKN2C','ERBB3','FLT3','IRS2','MEF2B','NPM1','PPP2R1A','SDHB','TNFRSF14', 'ALK','BTG1','CEBPA','ERBB4','FOXL2','JAK1','MEN1','NRAS','PPP2R2A','SDHC','TP53','ALOX12B','BTG2','CHEK1','ERCC4','FUBP1','JAK2','MERTK','NT5C2','PRDM1','SDHD','TSC1','AMER1','BTK','CHEK2','ERG','GABRA6','JAK3','MET','NTRK1','PRKAR1A','SETD2','TSC2','APC','C11orf30','CIC','ERRFI1','GATA3','JUN','MITF','NTRK2','PRKCI','SF3B1','TYRO3','AR','CALR','CREBBP','ESR1','GATA4','KDM5A','MKNK1','NTRK3','PTCH1','SGK1','U2AF1','ARAF','CARD11','CRKL','EZH2','GATA6','KDM5C','MLH1','P2RY8','PTEN','SMAD2','VEGFA','ARFRP1','CASP8','CSF1R','FAM46C','GID4','(C17orf39)','KDM6A','MPL','PALB2','PTPN11','SMAD4','VHL','ARID1A','CBFB','CSF3R','FANCA','GNA11','KDR','MRE11A','PARK2','PTPRO','SMARCA4','WHSC1','ASXL1','CBL','CTCF','FANCC','GNA13','KEAP1','MSH2','PARP1','QKI','SMARCB1','WHSC1L1','ATM','CCND1','CTNNA1','FANCG','GNAQ','KEL','MSH3','PARP2','RAC1','SMO','WT1','ATR','CCND2','CTNNB1','FANCL','GNAS','KIT','MSH6','PARP3','RAD21','SNCAIP','XPO1','ATRX','CCND3','CUL3','FAS','GRM3','KLHL6','MST1R','PAX5','RAD51','SOCS1','XRCC2','AURKA','CCNE1','CUL4A','FBXW7','GSK3B','KMT2A','MLL','MTAP','PBRM1','RAD51B','SOX2','ZNF217','AURKB','CD22','CXCR4','FGF10','H3F3A', 'KMT2D','MLL2','MTOR','PDCD1','RAD51C','SOX9','ZNF703','AXIN1','CD274','CYP17A1','FGF12','HDAC1','KRAS','MUTYH','PDCD1LG2','RAD51D','SPEN','AXL','CD70','DAXX','FGF14','HGF','LTK','MYC','PDGFRA','RAD52','SPOP','BAP1','CD79A','DDR1','FGF19','HNF1A','LYN','MYCL','PDGFRB','RAD54L','SRC','BARD1','CD79B','DDR2','FGF23','HRAS','MAF','MYCN','PDK1','RAF1','STAG2','BCL2','CDC73','DIS3','FGF3','HSD3B1','MAP2K1','MYD88','PIK3C2B','RARA','STAT3','BCL2L1','CDH1','DNMT3A','FGF4','ID3','MAP2K2','NBN','PIK3C2G','RB1','STK11','BCL2L2','CDK12','DOT1L','FGF6','IDH1','MAP2K4','NF1','PIK3CA','RBM10','SUFU','BCL6','CDK4','EED','FGFR1','IDH2','MAP3K1','NF2','PIK3CB','REL','SYK','BCOR','CDK6','EGFR','FGFR2','IGF1R','MAP3K13','NFE2L2','PIK3R1','RET','TBX3','BCORL1','CDK8','EP300','FGFR3','IKBKE','MAPK1','NFKBIA','PIM1','RICTOR','TEK')

#-------------------------------------------------------------------------------
# 1.2 Fetch regions with intronic targets
#-------------------------------------------------------------------------------
# Initialize panel df
panel <- data.frame()
# Read gencode annotation data (source: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz)
gencode <- read.table('/net/beegfs/cfg/tgac/jjanssen4/Resources/gencode/gencode.v19.annotation.REORDERED.gtf',sep =  "\t", col.names = c("chromosome_name","annotation_source","feature_type","genomic_start_location","genomic_end_location","score","genomic_strand","genomic_phase","info")) %>%
    tidyr::separate(info,
                    into = c("gene_id","transcript_id","gene_type","gene_status","gene_name","transcript_type","transcript_status","transcript_name","exon_number","exon_id","level"),
                    sep = ";") %>%
    dplyr::mutate(gene_name = gsub(" gene_name ","",gene_name)) 

# genes and intronic regions from Table 3 to fetch
Intron_regions <- data.frame(
    gene = c(
        c('ALK','ALK'),
        c('BRCA1','BRCA1','BRCA1','BRCA1','BRCA1','BRCA1','BRCA1'),
        c('ETV4'),
        c('EZR','EZR','EZR'),
        c('KIT'),
        c('MYC'),
        c('NUTM1'),
        c('RET','RET','RET','RET','RET'),
        c('SLC34A2'),
        c('BCL2'),
        c('BRCA2'),
        c('ETV5','ETV5'),
        c('FGFR1','FGFR1','FGFR1'),
        c('KMT2A','KMT2A','KMT2A','KMT2A','KMT2A','KMT2A'),
        c('NOTCH2'),
        c('PDGFRA','PDGFRA','PDGFRA'),
        c('ROS1','ROS1','ROS1','ROS1','ROS1'),
        c('TERC'),
        c('BCR','BCR','BCR'),
        c('CD74','CD74','CD74'),
        c('ETV6','ETV6'),
        c('FGFR2','FGFR2'),
        c('MSH2'),
        c('NTRK1','NTRK1','NTRK1','NTRK1'),
        c('RAF1','RAF1','RAF1','RAF1','RAF1'),
        c('RSPO2'),
        c('TERT'),
        c('BRAF','BRAF','BRAF','BRAF'),
        c('EGFR','EGFR','EGFR','EGFR','EGFR','EGFR'),
        c('EWSR1','EWSR1','EWSR1','EWSR1','EWSR1','EWSR1','EWSR1'),
        c('FGFR3'),
        c('MYB'),
        c('NTRK2'),
        c('RARA'),
        c('SDC4'),
        c('TMPRSS2','TMPRSS2','TMPRSS2')
        ),
    region = c(
        c('intron-18','intron-19'),
        c('intron-2','intron-7','intron-8','intron-12','intron-16','intron-19','intron-20'),
        c('intron-8'),
        c('intron-9','intron-10','intron-11'),
        c('intron-16'),
        c('intron-1'),
        c('intron-1'),
        c('intron-7','intron-8','intron-9','intron-10','intron-11'),
        c('intron-4'),
        c('3-UTR'),
        c('intron-2'),
        c('intron-6','intron-7'),
        c('intron-1','intron-5','intron-17'),
        c('intron-6','intron-7','intron-8','intron-9','intron-10','intron-11'),
        c('itnron-26'),
        c('intron-7','intron-9','intron-11'),
        c('intron-31','intron-32','intron-33','intron-34','intron-35'),
        c('ncRNA'),
        c('intron-8','intron-13','intron-14'),
        c('intron-6','intron-7','intron-8'),
        c('intron-5','intron-6'),
        c('intron-1','intron-17'),
        c('intron-5'),
        c('intron-8','intron-9','intron-10','intron-11'),
        c('intron-4','intron-5','intron-6','intron-7','intron-8'),
        c('intron-1'),
        c('promoter'),
        c('intron-7','intron-8','intron-9','intron-10'),
        c('intron-7','intron-15','intron-24','intron-25','intron-26','intron-27'),
        c('intron-7','intron-8','intron-9','intron-10','intron-11','intron-12','intron-13'),
        c('intron-17'),
        c('intron-14'),
        c('intron-12'),
        c('intron-2'),
        c('intron-2'),
        c('intron-1','intron-2','intron-3')
        )
)

# Iterate over intron regions and add panel regions
# For introns, take exon-X ends, and exon-X+1 starts as beginning and end respectively
for(i in seq(nrow(Intron_regions))){
    gene <- Intron_regions$gene[i]
    region <- Intron_regions$region[i]
    if(grepl('intron',region)){
        intron_number <- as.integer(strsplit(region,'-')[[1]][2])
        intron_region <-
            gencode %>%
            filter(gene_name == gene, c(exon_number == paste0(' exon_number ',intron_number) | exon_number == paste0(' exon_number ',intron_number+1))) %>% 
            # Select (the first) transcript with the two exon nunbers present
            group_by(transcript_id) %>%
            mutate(number_of_exons = n_distinct(exon_number)) %>%
            ungroup() %>%
            filter(number_of_exons > 1) %>% 
            filter(transcript_id == first(transcript_id)) %>% 
            select(chromosome_name,genomic_start_location,genomic_end_location,exon_number,genomic_strand) %>% 
            unique()
 
        if(unique(intron_region$genomic_strand) == '+'){
            start <- 'genomic_end_location'
            end <- 'genomic_start_location'
        }else{
            start <- 'genomic_start_location'
            end <- 'genomic_end_location'
        }
                            
        
        panel <-
            rbind(panel,
                  data.frame(
                      chromosome_name = unique(intron_region$chromosome_name), #Chromosome
                      genomic_start_location = intron_region[intron_region$exon_number == paste0(' exon_number ',intron_number), start], # start location
                      genomic_end_location = intron_region[intron_region$exon_number == paste0(' exon_number ',intron_number+1), end] # end location
                  ))
    }else if(region == 'promoter' & gene == 'TERT'){
        # Taken from genome browser (ENSMBL, bounds up)
        panel <- rbind(panel,
                       c('chr5',1293602,1296199))

    }else if(region == 'ncRNA'  & gene == 'TERC'){
        panel <-
            rbind(panel,
                  gencode %>%
                  filter(gene_name == gene, feature_type == 'gene') %>%
                  select(chromosome_name,genomic_start_location,genomic_end_location) %>%
                  unique()
                  )
    }else if (region == '3-UTR' & gene == 'BCL2'){
        panel <-
            rbind(panel,
                  gencode %>%
                  filter(gene_name == 'BCL2',grepl('UTR',feature_type), !grepl('5_UTR',level)) %>%
                  filter(genomic_end_location == max(genomic_end_location)) %>%
                  select(chromosome_name,genomic_start_location,genomic_end_location) %>%
                  unique()
                  )
        
    }
}

panel <- unique(panel)

# Add exon regions
panel <- rbind(panel,
               gencode %>% filter(gene_name %in% ExonGenes, feature_type == 'exon')  %>%
               select(chromosome_name,genomic_start_location,genomic_end_location) %>%
               unique())

colnames(panel) <- c('chr','chromStart','chromEnd')



panel <- panel %>%
    mutate(
        min = case_when(
            chromStart > chromEnd ~ chromEnd,
            TRUE ~ chromStart
        ),
        max =  case_when(
            chromStart < chromEnd ~ chromEnd,
            TRUE ~ chromStart
        )) %>%
    select(chr,min,max)

colnames(panel) <- c('chr','chromStart','chromEnd')



