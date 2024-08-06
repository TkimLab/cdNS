cdNS (context-depedendent dnds ratio scores) - Evolutionary selection of gene pairs with mutations in cancer genomes




#
# This code is modified from the original dndscv R pacakges (https://www.sanger.ac.uk/tool/dndscv/ by Martincorena et al 2017).
# The main purpose of the modification is (1) to run the code in a standalone mode, (2) reuse the resource to speed up the iteration
# The code has been used to calculate two dndscv scores for each gene in the list - 
# (1) dndscv of genes in genomes with mutations of the corresponding gene (context+)
# and (2) dndscv of genes in genomes without mutations of the corresponding gene (context-)
# So, cdNS scores (context-dependent dnds ratios) can be calculate for genes as ratios of them
# Three outputs will be generated (per gene) in default OUTPUT folder:
# - dndscv scores and related outputs of context= and context+ genomes: (tumortype)_OUTPUT_gene1_SELCV.txt
# - dndscv scores of context- genomes with permutation (default 100 runs): (tumortype)_OUTPUT_gene1_permW_POS.txt
# - dndscv scores of context+ genomes with permutation (default 100 runs): (tumortype)_OUTPUT_gene1_permW_NEG.txt
#

# The following files are NOT provided in Github, and will be given by request
#
# dndsoutput_ANNOTATION_totalTCGA.txt (275 Mb)
# dndsoutput_REFCDS_totalTCGA.rds (70.9 Mb)
#
