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

# The following options are default ones used by original dndscv R packages
gene_list = NULL
refdb = "hg19"
sm = "192r_3w"
kc = "cgc81"
cv = "hg19"
max_muts_per_gene_per_sample = 3
max_coding_muts_per_sample = 3000
use_indel_sites = T
min_indels = 5
maxcovs = 20
constrain_wnon_wspl = T
outp = 3
numcode = 1
outmats = F
mingenecovs = 500
dc = NULL 

# load the dataset of dndscv R packages
substmodel = readRDS("dndsoutput_substmodel.RDS")
covs = readRDS("dndsoutput_covs.RDS")
rownames(covs)[3522] = "CDKN2A"
known_cancergenes = readRDS("dndsoutput_known_cancergenes.RDS")
known_cancergenes = c(known_cancergenes, "CDKN2A")

# selection of tumor types
tumortypeTCGA = read.table("TCGA-CDR-tumortype_tcgaID.txt", header=T, sep="\t", row.names=1)
tumortype2 = tumortypeTCGA$type
tumortype2[tumortype2 == "LGG" | tumortype2 == "GBM"] = "LGGGBM"
tumortype2[tumortype2 == "COAD" | tumortype2 == "READ"] = "COADREAD"

tumortype_selected = c("ALL", "ESCA", "PAAD", "SARC", "KIRP", "CESC", "LIHC", "BLCA", "STAD", "SKCM", "PRAD", "LUSC", "THCA", "LUAD", "HNSC", "KIRC", "UCEC", "OV", "COADREAD", "BRCA", "LGGGBM")

# For PanCancer analysis, [1] ("ALL") should be designated, otherwise provide index for specific tumor types in 'tumortype_selected'
tumortype_target = tumortype_selected[1]

# TCGA mutation annotaiton
annot_total = read.table("dndsoutput_ANNOTATION_totalTCGA.txt", sep="\t", header=T)
annot_total$gene[annot_total$gene == "CDKN2A.p14arf" | annot_total$gene == "CDKN2A.p16INK4a"] = "CDKN2A"
if (tumortype_target != "ALL") {
annot_total = annot_total[substr(annot_total$sampleID, 1, 12) %in% rownames(tumortypeTCGA)[tumortype2 == tumortype_target],]	
}

# RefCDS is provided and can be reused
RefCDS = readRDS("dndsoutput_REFCDS_totalTCGA.rds")
RefCDS[[3522]]$gene_name = "CDKN2A"
genesRefCDS = sapply(RefCDS, function(x) x$gene_name)

#tempcall_forindels = substr(annot_total$ntchange, nchar(annot_total$ntchange)-9, nchar(annot_total$ntchange))
#truncating_total_call = ifelse(annot_total$impact %in% c("Essential_Splice", "Nonsense"), 1, ifelse(tempcall_forindels %in% c("delfrshift", "insfrshift"),1,0)) # 213189 (ind-fr 90440)
nonsilent_total_call = ifelse(annot_total$impact == "Synonymous", 0, 1)												# 1444822

# cdNS scores is dnds (N1) - dnds (N2) - for each mutation pair, dnds/N1 and dnds/N2 is calculated simultaneously
# refN_total for refN1 - refN2 
refN_total = array(0, dim=c(length(genesRefCDS), nrow(substmodel), 4))
triContext = paste0(annot_total$ref3_cod, ">", annot_total$mut3_cod)
syn_m = as.matrix(table(annot_total$gene[annot_total$impact == "Synonymous"], triContext[annot_total$impact == "Synonymous"]))
refN_total[match(rownames(syn_m),genesRefCDS), match(colnames(syn_m),rownames(substmodel)), 1] = syn_m
mis_m = as.matrix(table(annot_total$gene[annot_total$impact == "Missense"], triContext[annot_total$impact == "Missense"]))
refN_total[match(rownames(mis_m),genesRefCDS), match(colnames(mis_m),rownames(substmodel)), 2] = mis_m
non_m = as.matrix(table(annot_total$gene[annot_total$impact == "Nonsense"], triContext[annot_total$impact == "Nonsense"]))
refN_total[match(rownames(non_m),genesRefCDS), match(colnames(non_m),rownames(substmodel)), 3] = non_m
spl_m = as.matrix(table(annot_total$gene[annot_total$impact == "Essential_Splice"], triContext[annot_total$impact == "Essential_Splice"]))
refN_total[match(rownames(spl_m),genesRefCDS), match(colnames(spl_m),rownames(substmodel)), 4] = spl_m

# global
Lall = array(sapply(RefCDS, function(x) x$L), dim = c(192, 4, length(RefCDS)))	# 192X4X20091
L = apply(Lall, c(1, 2), sum)							# 192X4 (no sample) so global!

# for dNdSloc
fit_substmodel = function(N) {
	l = c(L)
	n = c(N)
	r = c(substmodel)
	n = n[l != 0]
	r = r[l != 0]
	l = l[l != 0]
	params = unique(base::strsplit(x = paste(r, collapse = "*"), split = "\\*")[[1]])
	indmat = as.data.frame(array(0, dim = c(length(r), length(params))))
	colnames(indmat) = params
	for (j in 1:length(r)) {
	    indmat[j, base::strsplit(r[j], split = "\\*")[[1]]] = 1
	}
	model = glm(formula = n ~ offset(log(l)) + . - 1, data = indmat, family = poisson(link = log))
	mle = exp(coefficients(model))
	ci = exp(confint.default(model))
	par = data.frame(name = gsub("`", "", rownames(ci)), mle = mle[rownames(ci)], cilow = ci[, 1], cihigh = ci[, 2])
	return(par)
}

# for dNdScv
mle_tcv = function(n_neutral, exp_rel_neutral, shape, scale) {
    tml = (n_neutral + shape - 1)/(exp_rel_neutral + (1/scale))
    if (shape <= 1) { tml = max(shape * scale, tml) }
    return(tml)
}
selfun_cv1 = function(j) {
    y = as.numeric(genemuts1[j, -1])
    exp_rel = y[5:8]/y[5]
    shape = theta1
    scale = y[9]/theta1
    indneut = 1
    opt_t = mle_tcv(n_neutral = sum(y[indneut]), exp_rel_neutral = sum(exp_rel[indneut]), shape = shape, scale = scale)
    mrfold = max(1e-10, opt_t/sum(y[5]))
    wfree = y[2:4]/y[6:8]/mrfold
    wfree[y[2:4] == 0] = 0
    if (constrain_wnon_wspl == 0) {
	return(c(wfree))
    } else {
	wmisfree = y[2]/y[6]/mrfold
	wmisfree[y[2] == 0] = 0
	wtruncfree = sum(y[3:4])/sum(y[7:8])/mrfold
	wtruncfree[sum(y[3:4]) == 0] = 0
	return(c(wmisfree, wtruncfree, wtruncfree))
    }
}
selfun_cv2 = function(j) {
    y = as.numeric(genemuts2[j, -1])
    exp_rel = y[5:8]/y[5]
    shape = theta2
    scale = y[9]/theta2
    indneut = 1
    opt_t = mle_tcv(n_neutral = sum(y[indneut]), exp_rel_neutral = sum(exp_rel[indneut]), shape = shape, scale = scale)
    mrfold = max(1e-10, opt_t/sum(y[5]))
    wfree = y[2:4]/y[6:8]/mrfold
    wfree[y[2:4] == 0] = 0
    if (constrain_wnon_wspl == 0) {
	return(c(wfree))
    } else {
	wmisfree = y[2]/y[6]/mrfold
	wmisfree[y[2] == 0] = 0
	wtruncfree = sum(y[3:4])/sum(y[7:8])/mrfold
	wtruncfree[sum(y[3:4]) == 0] = 0
	return(c(wmisfree, wtruncfree, wtruncfree))
    }
}


#
# main function - read gene list and perform dndscv analysis for context(+) genomes and context(-) genomes
#		cf. context(+) and (-) genomes are those with and without the mutations in the corresponding genes
#		- provided example is the 312 genes with high level of dnds ratios (likey to be evolutionarily selected genes)
#		- the results will be (1) dndscv ratios of context(+) and context(-) genomes 
#		- ** cdNS (context-dependent dnds ratios) scores of a gene can be calculated as dnds(context+)/dnds(context-)
#		- higher cdNS scores represent that the context mutation has conferred positive evolutionary selective pressure on the mutations of the corresponding gene
#		(2) and (3) results will be permutated results (default: 100 permutations) from which nominal P values of cdNS scores can be estimated

genes_consensus_template = read.table("INPUT_contextgenes.txt", header=T, sep="\t")
genes_consensus = unique(genes_consensus_template$contextgene)

for (context_gene_index in 1:nrow(genes_consensus_template)) {

	context_gene = genes_consensus[context_gene_index]
	print (paste("----", context_gene_index, "/", length(genes_consensus), "-", format(Sys.time(), "%a %b %d %X %Y"))); flush.console();

	context_genomes = unique(annot_total[annot_total$gene == context_gene & nonsilent_total_call,]$sampleID)

	# sel_loc
	whetherContext = annot_total$sampleID %in% context_genomes
	annot = annot_total[whetherContext,]		
	#RefCDS = RefCDS_total
	#genesRefCDS = sapply(RefCDS, function(x) x$gene_name)

	triContext = paste0(annot$ref3_cod, ">", annot$mut3_cod)
	refN1 = array(0, dim=c(length(genesRefCDS), nrow(substmodel), 4))
	syn_m = as.matrix(table(annot$gene[annot$impact == "Synonymous"], triContext[annot$impact == "Synonymous"]))
	refN1[match(rownames(syn_m),genesRefCDS), match(colnames(syn_m),rownames(substmodel)), 1] = syn_m
	mis_m = as.matrix(table(annot$gene[annot$impact == "Missense"], triContext[annot$impact == "Missense"]))
	refN1[match(rownames(mis_m),genesRefCDS), match(colnames(mis_m),rownames(substmodel)), 2] = mis_m
	non_m = as.matrix(table(annot$gene[annot$impact == "Nonsense"], triContext[annot$impact == "Nonsense"]))
	refN1[match(rownames(non_m),genesRefCDS), match(colnames(non_m),rownames(substmodel)), 3] = non_m
	spl_m = as.matrix(table(annot$gene[annot$impact == "Essential_Splice"], triContext[annot$impact == "Essential_Splice"]))
	refN1[match(rownames(spl_m),genesRefCDS), match(colnames(spl_m),rownames(substmodel)), 4] = spl_m

	refN2 = array(0, dim=c(length(genesRefCDS), nrow(substmodel), 4))
	refN2 = refN_total - refN1

	####################################################################################################
	# GLOBAL-DNDS - prepare 'genemuts' and 'mutrates'
	N1 = apply(refN1, c(2,3), sum)	# context+
	N2 = apply(refN2, c(2,3), sum)	# context-

	par1 = fit_substmodel(N1); parmle1 = setNames(par1[, 2], par1[, 1])
	genemuts1 = data.frame(gene_name = sapply(RefCDS, function(x) x$gene_name), n_syn = NA, n_mis = NA, n_non = NA, n_spl = NA, exp_syn = NA, exp_mis = NA, exp_non = NA, exp_spl = NA, stringsAsFactors = F)
	genemuts1[, 2:5] = apply(refN1, c(1,3), sum)	
	mutrates1 = sapply(substmodel[, 1], function(x) prod(parmle1[base::strsplit(x, split = "\\*")[[1]]]))	# mutrates based on "syn" mutations
	genemuts1[, 6:9] = t(sapply(RefCDS, function(x) colSums(x$L * mutrates1)))	
	numrates1 = length(mutrates1)

	par2 = fit_substmodel(N2); parmle2 = setNames(par2[, 2], par2[, 1])
	genemuts2 = data.frame(gene_name = sapply(RefCDS, function(x) x$gene_name), n_syn = NA, n_mis = NA, n_non = NA, n_spl = NA, exp_syn = NA, exp_mis = NA, exp_non = NA, exp_spl = NA, stringsAsFactors = F)
	genemuts2[, 2:5] = apply(refN2, c(1, 3), sum)	
	mutrates2 = sapply(substmodel[, 1], function(x) prod(parmle2[base::strsplit(x, split = "\\*")[[1]]]))	# mutrates based on "syn" mutations
	genemuts2[, 6:9] = t(sapply(RefCDS, function(x) colSums(x$L * mutrates2)))	
	numrates2 = length(mutrates2)

	# DNDScv
	nbrdf1 = cbind(genemuts1[, c("n_syn", "exp_syn")], covs)
	model1 = tryCatch({
	  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) + ., data = nbrdf1)
	}, warning = function(w) {
	  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf1)
	}, error = function(e) {
	  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf1)
	})
	nbrdf2 = cbind(genemuts2[, c("n_syn", "exp_syn")], covs)
	model2 = tryCatch({
	  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) + ., data = nbrdf2)
	}, warning = function(w) {
	  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf2)
	}, error = function(e) {
	  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf2)
	})
	if (all(model1$y == genemuts1$n_syn)) {genemuts1$exp_syn_cv = model1$fitted.values}
	if (all(model2$y == genemuts2$n_syn)) {genemuts2$exp_syn_cv = model2$fitted.values}
	theta1 = model1$theta; nbreg1 = model1
	theta2 = model2$theta; nbreg2 = model2

	sel_cv1 = as.data.frame(t(sapply(1:nrow(genemuts1), selfun_cv1)))
	colnames(sel_cv1) = c("wmis_cv", "wnon_cv", "wspl_cv")
	sel_cv1 = cbind(genemuts1, sel_cv1)
	sel_cv2 = as.data.frame(t(sapply(1:nrow(genemuts2), selfun_cv2)))
	colnames(sel_cv2) = c("wmis_cv", "wnon_cv", "wspl_cv")
	sel_cv2 = cbind(genemuts2, sel_cv2)

	# INDEL - context_
	indels = annot[annot$impact == "no-SNV",]
	if (nrow(indels) >= min_indels) {						#<<- indels
	    geneindels = as.data.frame(array(0, dim = c(length(RefCDS), 8)))
	    colnames(geneindels) = c("gene_name", "n_ind", "n_induniq", "n_indused", "cds_length", "excl", "exp_unif", "exp_indcv")
	    geneindels$gene_name = sapply(RefCDS, function(x) x$gene_name)
	    geneindels$n_ind = as.numeric(table(indels$gene)[geneindels[, 1]])
	    geneindels[is.na(geneindels[, 2]), 2] = 0
	    geneindels$n_induniq = as.numeric(table(unique(indels[, -1])$gene)[geneindels[, 1]])
	    geneindels[is.na(geneindels[, 3]), 3] = 0
	    geneindels$cds_length = sapply(RefCDS, function(x) x$CDS_length)
	    if (!is.null(dc)) {
		geneindels$cds_length = geneindels$cds_length * dc
	    }
	    if (use_indel_sites) {
		geneindels$n_indused = geneindels[, 3]
	    } else {
		geneindels$n_indused = geneindels[, 2]
	    }
	    geneindels$excl = (geneindels[, 1] %in% known_cancergenes)
	    min_bkg_genes = 50
	    if (sum(!geneindels$excl) < min_bkg_genes | sum(geneindels[!geneindels$excl, "n_indused"]) == 0) {
		newkc = as.vector(sel_cv$gene_name[sel_cv1$qallsubs_cv < 0.01])
		geneindels$excl = (geneindels[, 1] %in% newkc)
		if (sum(!geneindels$excl) < min_bkg_genes | sum(geneindels[!geneindels$excl, "n_indused"]) == 0) {
		  geneindels$excl = Fh
		  message("    No gene was excluded from the background indel model.")
		} else {
		  warning(sprintf("    Genes were excluded from the indel background model based on the substitution data: %s.", paste(newkc, collapse = ", ")))
		}
	    }
	    geneindels$exp_unif = sum(geneindels[!geneindels$excl, "n_indused"])/sum(geneindels[!geneindels$excl, "cds_length"]) * geneindels$cds_length
	    if (is.null(cv)) {
		nbrdf = geneindels[, c("n_indused", "exp_unif")][!geneindels[, 6], ]
		model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
		nbrdf_all = geneindels[, c("n_indused", "exp_unif")]
	    } else {
		nbrdf = cbind(geneindels[, c("n_indused", "exp_unif")], covs)[!geneindels[, 6], ]
		if (sum(!geneindels$excl) < 500) {
		  model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
		}  else {
		  model = tryCatch({
		    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) + ., data = nbrdf)
		  }, warning = function(w) {
		    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
		  }, error = function(e) {
		    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
		  })
		}
		nbrdf_all = cbind(geneindels[, c("n_indused", "exp_unif")], covs)
	    }
	    #message(sprintf("    Regression model for indels (theta = %0.3g)", model$theta))
	    theta_indels = model$theta
	    nbregind = model
	    geneindels$exp_indcv = exp(predict(model, nbrdf_all))
	    geneindels$wind = geneindels$n_indused/geneindels$exp_indcv
	    sel_cv1 = merge(sel_cv1, geneindels, by = "gene_name")[, c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_indused", "wmis_cv", "wnon_cv", "wspl_cv", "wind")]
	    colnames(sel_cv1) = c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_ind", "wmis_cv", "wnon_cv", "wspl_cv", "wind_cv")
	}
	indels = annot_total[!whetherContext & annot_total$impact == "no-SNV",]
	if (nrow(indels) >= min_indels) {						#<<- indels
	    geneindels = as.data.frame(array(0, dim = c(length(RefCDS), 8)))
	    colnames(geneindels) = c("gene_name", "n_ind", "n_induniq", "n_indused", "cds_length", "excl", "exp_unif", "exp_indcv")
	    geneindels$gene_name = sapply(RefCDS, function(x) x$gene_name)
	    geneindels$n_ind = as.numeric(table(indels$gene)[geneindels[, 1]])
	    geneindels[is.na(geneindels[, 2]), 2] = 0
	    geneindels$n_induniq = as.numeric(table(unique(indels[, -1])$gene)[geneindels[, 1]])
	    geneindels[is.na(geneindels[, 3]), 3] = 0
	    geneindels$cds_length = sapply(RefCDS, function(x) x$CDS_length)
	    if (!is.null(dc)) {
		geneindels$cds_length = geneindels$cds_length * dc
	    }
	    if (use_indel_sites) {
		geneindels$n_indused = geneindels[, 3]
	    } else {
		geneindels$n_indused = geneindels[, 2]
	    }
	    geneindels$excl = (geneindels[, 1] %in% known_cancergenes)
	    min_bkg_genes = 50
	    if (sum(!geneindels$excl) < min_bkg_genes | sum(geneindels[!geneindels$excl, "n_indused"]) == 0) {
		newkc = as.vector(sel_cv$gene_name[sel_cv2$qallsubs_cv < 0.01])
		geneindels$excl = (geneindels[, 1] %in% newkc)
		if (sum(!geneindels$excl) < min_bkg_genes | sum(geneindels[!geneindels$excl, "n_indused"]) == 0) {
		  geneindels$excl = Fh
		  message("    No gene was excluded from the background indel model.")
		} else {
		  warning(sprintf("    Genes were excluded from the indel background model based on the substitution data: %s.", paste(newkc, collapse = ", ")))
		}
	    }
	    geneindels$exp_unif = sum(geneindels[!geneindels$excl, "n_indused"])/sum(geneindels[!geneindels$excl, "cds_length"]) * geneindels$cds_length
	    if (is.null(cv)) {
		nbrdf = geneindels[, c("n_indused", "exp_unif")][!geneindels[, 6], ]
		model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
		nbrdf_all = geneindels[, c("n_indused", "exp_unif")]
	    } else {
		nbrdf = cbind(geneindels[, c("n_indused", "exp_unif")], covs)[!geneindels[, 6], ]
		if (sum(!geneindels$excl) < 500) {
		  model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
		}  else {
		  model = tryCatch({
		    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) + ., data = nbrdf)
		  }, warning = function(w) {
		    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
		  }, error = function(e) {
		    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
		  })
		}
		nbrdf_all = cbind(geneindels[, c("n_indused", "exp_unif")], covs)
	    }
	    #message(sprintf("    Regression model for indels (theta = %0.3g)", model$theta))
	    theta_indels = model$theta
	    nbregind = model
	    geneindels$exp_indcv = exp(predict(model, nbrdf_all))
	    geneindels$wind = geneindels$n_indused/geneindels$exp_indcv
	    sel_cv2 = merge(sel_cv2, geneindels, by = "gene_name")[, c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_indused", "wmis_cv", "wnon_cv", "wspl_cv", "wind")]
	    colnames(sel_cv2) = c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_ind", "wmis_cv", "wnon_cv", "wspl_cv", "wind_cv")
	}


	###################################### PERMUTATION ###########################################################################################################################
	sel_cv1_template = sel_cv1
	sel_cv2_template = sel_cv2

	permtime = 100
	wmis_matrix_permP = matrix(nrow=nrow(sel_cv1), ncol=permtime); colnames(wmis_matrix_permP) = paste0("wmisDNDScv", 1:permtime)
	wnon_matrix_permP = matrix(nrow=nrow(sel_cv1), ncol=permtime); colnames(wnon_matrix_permP) = paste0("wnonDNDScv", 1:permtime)
	wind_matrix_permP = matrix(nrow=nrow(sel_cv1), ncol=permtime); colnames(wind_matrix_permP) = paste0("windDNDScv", 1:permtime)
	wmis_matrix_permN = matrix(nrow=nrow(sel_cv1), ncol=permtime); colnames(wmis_matrix_permN) = paste0("wmisDNDScv", 1:permtime)
	wnon_matrix_permN = matrix(nrow=nrow(sel_cv1), ncol=permtime); colnames(wnon_matrix_permN) = paste0("wnonDNDScv", 1:permtime)
	wind_matrix_permN = matrix(nrow=nrow(sel_cv1), ncol=permtime); colnames(wind_matrix_permN) = paste0("windDNDScv", 1:permtime)

	for (perm in 1:permtime) {

		rm("refN1"); rm("refN2"); gc()

		print (paste(perm, "/", permtime, ":", format(Sys.time(), "%a %b %d %X %Y"))); flush.console();; flush.console()

		context_genomes = unique(annot_total$sampleID)[sample(1:length(unique(annot_total$sampleID)), length(context_genomes), replace=F)]

		# sel_loc
		whetherContext = annot_total$sampleID %in% context_genomes
		annot = annot_total[whetherContext,]		
		triContext = paste0(annot$ref3_cod, ">", annot$mut3_cod)
		refN1 = array(0, dim=c(length(genesRefCDS), nrow(substmodel), 4))

		syn_m = as.matrix(table(annot$gene[annot$impact == "Synonymous"], triContext[annot$impact == "Synonymous"]))
		refN1[match(rownames(syn_m),genesRefCDS), match(colnames(syn_m),rownames(substmodel)), 1] = syn_m
		mis_m = as.matrix(table(annot$gene[annot$impact == "Missense"], triContext[annot$impact == "Missense"]))
		refN1[match(rownames(mis_m),genesRefCDS), match(colnames(mis_m),rownames(substmodel)), 2] = mis_m
		non_m = as.matrix(table(annot$gene[annot$impact == "Nonsense"], triContext[annot$impact == "Nonsense"]))
		refN1[match(rownames(non_m),genesRefCDS), match(colnames(non_m),rownames(substmodel)), 3] = non_m
		spl_m = as.matrix(table(annot$gene[annot$impact == "Essential_Splice"], triContext[annot$impact == "Essential_Splice"]))
		refN1[match(rownames(spl_m),genesRefCDS), match(colnames(spl_m),rownames(substmodel)), 4] = spl_m

		refN2 = array(0, dim=c(length(genesRefCDS), nrow(substmodel), 4))
		refN2 = refN_total - refN1

		####################################################################################################
		# GLOBAL-DNDS - prepare 'genemuts' and 'mutrates'
		N1 = apply(refN1, c(2,3), sum)	# context+
		N2 = apply(refN2, c(2,3), sum)	# context-

		par1 = fit_substmodel(N1); parmle1 = setNames(par1[, 2], par1[, 1])
		genemuts1 = data.frame(gene_name = sapply(RefCDS, function(x) x$gene_name), n_syn = NA, n_mis = NA, n_non = NA, n_spl = NA, exp_syn = NA, exp_mis = NA, exp_non = NA, exp_spl = NA, stringsAsFactors = F)
		genemuts1[, 2:5] = apply(refN1, c(1,3), sum)	
		mutrates1 = sapply(substmodel[, 1], function(x) prod(parmle1[base::strsplit(x, split = "\\*")[[1]]]))	# mutrates based on "syn" mutations
		genemuts1[, 6:9] = t(sapply(RefCDS, function(x) colSums(x$L * mutrates1)))	
		numrates1 = length(mutrates1)

		par2 = fit_substmodel(N2); parmle2 = setNames(par2[, 2], par2[, 1])
		genemuts2 = data.frame(gene_name = sapply(RefCDS, function(x) x$gene_name), n_syn = NA, n_mis = NA, n_non = NA, n_spl = NA, exp_syn = NA, exp_mis = NA, exp_non = NA, exp_spl = NA, stringsAsFactors = F)
		genemuts2[, 2:5] = apply(refN2, c(1, 3), sum)	
		mutrates2 = sapply(substmodel[, 1], function(x) prod(parmle2[base::strsplit(x, split = "\\*")[[1]]]))	# mutrates based on "syn" mutations
		genemuts2[, 6:9] = t(sapply(RefCDS, function(x) colSums(x$L * mutrates2)))	
		numrates2 = length(mutrates2)

		# DNDScv
		nbrdf1 = cbind(genemuts1[, c("n_syn", "exp_syn")], covs)
		model1 = tryCatch({
		  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) + ., data = nbrdf1)
		}, warning = function(w) {
		  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf1)
		}, error = function(e) {
		  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf1)
		})
		nbrdf2 = cbind(genemuts2[, c("n_syn", "exp_syn")], covs)
		model2 = tryCatch({
		  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) + ., data = nbrdf2)
		}, warning = function(w) {
		  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf2)
		}, error = function(e) {
		  MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1, data = nbrdf2)
		})
		if (all(model1$y == genemuts1$n_syn)) {genemuts1$exp_syn_cv = model1$fitted.values}
		if (all(model2$y == genemuts2$n_syn)) {genemuts2$exp_syn_cv = model2$fitted.values}
		theta1 = model1$theta; nbreg1 = model1
		theta2 = model2$theta; nbreg2 = model2

		sel_cv1 = as.data.frame(t(sapply(1:nrow(genemuts1), selfun_cv1)))
		colnames(sel_cv1) = c("wmis_cv", "wnon_cv", "wspl_cv")
		sel_cv1 = cbind(genemuts1, sel_cv1)
		sel_cv2 = as.data.frame(t(sapply(1:nrow(genemuts2), selfun_cv2)))
		colnames(sel_cv2) = c("wmis_cv", "wnon_cv", "wspl_cv")
		sel_cv2 = cbind(genemuts2, sel_cv2)

		# INDEL - context_
		indels = annot[annot$impact == "no-SNV",]
		if (nrow(indels) >= min_indels) {						#<<- indels
		    geneindels = as.data.frame(array(0, dim = c(length(RefCDS), 8)))
		    colnames(geneindels) = c("gene_name", "n_ind", "n_induniq", "n_indused", "cds_length", "excl", "exp_unif", "exp_indcv")
		    geneindels$gene_name = sapply(RefCDS, function(x) x$gene_name)
		    geneindels$n_ind = as.numeric(table(indels$gene)[geneindels[, 1]])
		    geneindels[is.na(geneindels[, 2]), 2] = 0
		    geneindels$n_induniq = as.numeric(table(unique(indels[, -1])$gene)[geneindels[, 1]])
		    geneindels[is.na(geneindels[, 3]), 3] = 0
		    geneindels$cds_length = sapply(RefCDS, function(x) x$CDS_length)
		    if (!is.null(dc)) {
			geneindels$cds_length = geneindels$cds_length * dc
		    }
		    if (use_indel_sites) {
			geneindels$n_indused = geneindels[, 3]
		    } else {
			geneindels$n_indused = geneindels[, 2]
		    }
		    geneindels$excl = (geneindels[, 1] %in% known_cancergenes)
		    min_bkg_genes = 50
		    if (sum(!geneindels$excl) < min_bkg_genes | sum(geneindels[!geneindels$excl, "n_indused"]) == 0) {
			newkc = as.vector(sel_cv$gene_name[sel_cv1$qallsubs_cv < 0.01])
			geneindels$excl = (geneindels[, 1] %in% newkc)
			if (sum(!geneindels$excl) < min_bkg_genes | sum(geneindels[!geneindels$excl, "n_indused"]) == 0) {
			  geneindels$excl = Fh
			  message("    No gene was excluded from the background indel model.")
			} else {
			  warning(sprintf("    Genes were excluded from the indel background model based on the substitution data: %s.", paste(newkc, collapse = ", ")))
			}
		    }
		    geneindels$exp_unif = sum(geneindels[!geneindels$excl, "n_indused"])/sum(geneindels[!geneindels$excl, "cds_length"]) * geneindels$cds_length
		    if (is.null(cv)) {
			nbrdf = geneindels[, c("n_indused", "exp_unif")][!geneindels[, 6], ]
			model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
			nbrdf_all = geneindels[, c("n_indused", "exp_unif")]
		    } else {
			nbrdf = cbind(geneindels[, c("n_indused", "exp_unif")], covs)[!geneindels[, 6], ]
			if (sum(!geneindels$excl) < 500) {
			  model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
			}  else {
			  model = tryCatch({
			    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) + ., data = nbrdf)
			  }, warning = function(w) {
			    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
			  }, error = function(e) {
			    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
			  })
			}
			nbrdf_all = cbind(geneindels[, c("n_indused", "exp_unif")], covs)
		    }
		    #message(sprintf("    Regression model for indels (theta = %0.3g)", model$theta))
		    theta_indels = model$theta
		    nbregind = model
		    geneindels$exp_indcv = exp(predict(model, nbrdf_all))
		    geneindels$wind = geneindels$n_indused/geneindels$exp_indcv
		    sel_cv1 = merge(sel_cv1, geneindels, by = "gene_name")[, c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_indused", "wmis_cv", "wnon_cv", "wspl_cv", "wind")]
		    colnames(sel_cv1) = c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_ind", "wmis_cv", "wnon_cv", "wspl_cv", "wind_cv")
		}
		indels = annot_total[!whetherContext & annot_total$impact == "no-SNV",]
		if (nrow(indels) >= min_indels) {						#<<- indels
		    geneindels = as.data.frame(array(0, dim = c(length(RefCDS), 8)))
		    colnames(geneindels) = c("gene_name", "n_ind", "n_induniq", "n_indused", "cds_length", "excl", "exp_unif", "exp_indcv")
		    geneindels$gene_name = sapply(RefCDS, function(x) x$gene_name)
		    geneindels$n_ind = as.numeric(table(indels$gene)[geneindels[, 1]])
		    geneindels[is.na(geneindels[, 2]), 2] = 0
		    geneindels$n_induniq = as.numeric(table(unique(indels[, -1])$gene)[geneindels[, 1]])
		    geneindels[is.na(geneindels[, 3]), 3] = 0
		    geneindels$cds_length = sapply(RefCDS, function(x) x$CDS_length)
		    if (!is.null(dc)) {
			geneindels$cds_length = geneindels$cds_length * dc
		    }
		    if (use_indel_sites) {
			geneindels$n_indused = geneindels[, 3]
		    } else {
			geneindels$n_indused = geneindels[, 2]
		    }
		    geneindels$excl = (geneindels[, 1] %in% known_cancergenes)
		    min_bkg_genes = 50
		    if (sum(!geneindels$excl) < min_bkg_genes | sum(geneindels[!geneindels$excl, "n_indused"]) == 0) {
			newkc = as.vector(sel_cv$gene_name[sel_cv2$qallsubs_cv < 0.01])
			geneindels$excl = (geneindels[, 1] %in% newkc)
			if (sum(!geneindels$excl) < min_bkg_genes | sum(geneindels[!geneindels$excl, "n_indused"]) == 0) {
			  geneindels$excl = Fh
			  message("    No gene was excluded from the background indel model.")
			} else {
			  warning(sprintf("    Genes were excluded from the indel background model based on the substitution data: %s.", paste(newkc, collapse = ", ")))
			}
		    }
		    geneindels$exp_unif = sum(geneindels[!geneindels$excl, "n_indused"])/sum(geneindels[!geneindels$excl, "cds_length"]) * geneindels$cds_length
		    if (is.null(cv)) {
			nbrdf = geneindels[, c("n_indused", "exp_unif")][!geneindels[, 6], ]
			model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
			nbrdf_all = geneindels[, c("n_indused", "exp_unif")]
		    } else {
			nbrdf = cbind(geneindels[, c("n_indused", "exp_unif")], covs)[!geneindels[, 6], ]
			if (sum(!geneindels$excl) < 500) {
			  model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
			}  else {
			  model = tryCatch({
			    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) + ., data = nbrdf)
			  }, warning = function(w) {
			    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
			  }, error = function(e) {
			    MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf)
			  })
			}
			nbrdf_all = cbind(geneindels[, c("n_indused", "exp_unif")], covs)
		    }
		    #message(sprintf("    Regression model for indels (theta = %0.3g)", model$theta))
		    theta_indels = model$theta
		    nbregind = model
		    geneindels$exp_indcv = exp(predict(model, nbrdf_all))
		    geneindels$wind = geneindels$n_indused/geneindels$exp_indcv
		    sel_cv2 = merge(sel_cv2, geneindels, by = "gene_name")[, c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_indused", "wmis_cv", "wnon_cv", "wspl_cv", "wind")]
		    colnames(sel_cv2) = c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_ind", "wmis_cv", "wnon_cv", "wspl_cv", "wind_cv")
		}
		wmis_matrix_permP[,perm] = sel_cv1$wmis_cv; wnon_matrix_permP[,perm] = sel_cv1$wnon_cv; wind_matrix_permP[,perm] = sel_cv1$wind_cv; 
		wmis_matrix_permN[,perm] = sel_cv2$wmis_cv; wnon_matrix_permN[,perm] = sel_cv2$wnon_cv; wind_matrix_permN[,perm] = sel_cv2$wind_cv; 

	}

	if (FALSE) {	# maybe wrong for dndscv==0 cases!
		p_mis1 = p_trunc1 = p_ind1 = rep(NA, length(genesRefCDS))
		p_mis2 = p_trunc2 = p_ind2 = rep(NA, length(genesRefCDS))
		for (i in 1:length(genesRefCDS)) {
			if (sel_cv1_template[i,7] > 0 & sel_cv2_template[i,7] > 0) {
				r_mis = sel_cv1_template[i,7]/sel_cv2_template[i,7]
				r_mis_perm = wmis_matrix_permP[i,] / wmis_matrix_permN[i,]
				p_mis1[i] = sum(r_mis < r_mis_perm)/permtime; p_mis2[i] = sum(r_mis > r_mis_perm)/permtime
			}
			if (sel_cv1_template[i,8] > 0 & sel_cv2_template[i,8] > 0) {
				r_trunc = sel_cv1_template[i,8]/sel_cv2_template[i,8]
				r_trunc_perm = wnon_matrix_permP[i,] / wnon_matrix_permN[i,]
				p_trunc1[i] = sum(r_trunc < r_trunc_perm)/permtime; p_trunc2[i] = sum(r_trunc > r_trunc_perm)/permtime
			}
			if (sel_cv1_template[i,10] > 0 & sel_cv2_template[i,10] > 0) {
				r_ind = sel_cv1_template[i,10]/sel_cv2_template[i,10]
				r_ind_perm = wind_matrix_permP[i,] / wind_matrix_permN[i,]
				p_ind1[i] = sum(r_ind < r_ind_perm)/permtime; p_ind2[i] = sum(r_ind > r_ind_perm)/permtime
			}
		}
		#write.table(cbind(sel_cv1_template, sel_cv2_template, p_mis1, p_mis2, p_trunc1, p_trunc2, p_ind1, p_ind2), paste0("OUTPUT_517PERM_ver6_v3/OUTPUT_gene_", context_gene_index, "_contextP_perm_", permtime, ".txt"), sep="\t", quote=F)
		write.table(cbind(sel_cv1_template, sel_cv2_template, p_mis1, p_mis2, p_trunc1, p_trunc2, p_ind1, p_ind2), paste0("OUTPUT_gene_", context_gene_index, "_contextP_perm_", permtime, ".txt"), sep="\t", quote=F)
	}

	write.table(cbind(sel_cv1_template, sel_cv2_template), paste0("OUTPUT/", tumortype_target, "_OUTPUT_gene", context_gene_index, "_SELCV.txt"), sep="\t", quote=F)
	write.table(cbind(wmis_matrix_permP, wnon_matrix_permP, wind_matrix_permP), paste0("OUTPUT/", tumortype_target, "_OUTPUT_gene", context_gene_index, "_permW_POS.txt"), sep="\t", quote=F)
	write.table(cbind(wmis_matrix_permN, wnon_matrix_permN, wind_matrix_permN), paste0("OUTPUT/", tumortype_target, "_OUTPUT_gene", context_gene_index, "_permW_NEG.txt"), sep="\t", quote=F)

}
