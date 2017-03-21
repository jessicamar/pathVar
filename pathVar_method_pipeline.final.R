#load("~/pways_reactome.RData")
#load("~/pways_kegg.RData")
#load("~/dat_mat.RData")

#library("AnnotationDbi")
#library(org.Hs.eg.db)
#library(EMT)
#library(mclust)
#library(ggplot2)
#library(gridExtra)
#library(Matching)
#library(reshape)


#Kegg 291
#Reactome 1547


#######################################################
#makeDBList put your pathways text file into a list
#pway$PATHNAME is the pathway names from the file
#pway$PATHID is a vector of pathway ID numbers is there are any. Otherwise it will be a vector filled with NA
#pway$GENES is a list of vectors, where each vector are the genes for a single pathway
#######################################################
#file is a tab delimited text file, where first and second columns are pathwayID and pathway name. The third (or last column is the genes associated with each pathway, seperated by commas.
#pID is a boolean where options are TRUE or FALSE. If TRUE then the first column of the tab delimited file is expected to be path IDs. If FALSE, then the first column is expected to be pathway names.

makeDBList <- function(file, pID = TRUE) {
    pwayTable <- read.table(file, header = TRUE, sep = "\t", quote = "")
    pways <- NULL
    if (pID == TRUE) {
        pwayTable <- pwayTable[c(1, 2, ncol(pwayTable))]
        pways$PATHNAME <- as.vector(pwayTable[, 2])
        pways$PATHID <- as.vector(pwayTable[, 1])
        pways$GENES <- list(length(pways$PATHID))
        pways$SIZE <- NULL
        for (i in 1:length(pways$PATHID)) {
            geneVec <- unlist(strsplit(as.character(pwayTable[i, 3]), "[,/]"))
            if (length(as.vector(geneVec)) > 0) {
                pways$GENES[[i]] <- as.vector(geneVec)
                pways$SIZE[i] <- length(as.vector(geneVec))
            } else {
                pways$PATHID <- pways$PATHID[-i]
                pways$PATHNAME <- pways$PATHNAME[-i]
            }
        }
    } else {
        pwayTable <- pwayTable[c(1, ncol(pwayTable))]
        pways$PATHNAME <- as.vector(pwayTable[, 1])
        pways$PATHID <- rep("NA", length(pways$PATHNAME))
        pways$GENES <- list(length(pways$PATHID))
        pways$SIZE <- NULL
        for (i in 1:length(pways$PATHID)) {
            geneVec <- unlist(strsplit(as.character(pwayTable[i, 2]), "[,/]"))
            if (length(as.vector(geneVec)) > 0) {
                pways$GENES[[i]] <- as.vector(geneVec)
                pways$SIZE[i] <- length(as.vector(geneVec))
            } else {
                pways$PATHID <- pways$PATHID[-i]
                pways$PATHNAME <- pways$PATHNAME[-i]
            }
        }
    }
    names(pways$GENES) <- pways$PATHNAME
    return(pways)
}

#######################################################
#pathVarOneSample
#######################################################
# 1. Compute the statistics (sd, mad, cv or mean) for each gene. 
# 2. Classify the genes with respect to the statistics in at most 4 clusters.
# 3. For each pathway, we extract the gene in our dataset and in which cluster they belong.
# 4. For each pathway, we look how the gene counts in each category and compare it to the reference counts with all the gene from the dataset with the Chi-Squared or exact test.
# 5. We build a dataframe for 4 (description as Output 1 below).
#######################################################
####Output
# Output 1: tablePway columns are :pathway name, pathway IDs, adjusted p-value from the chisq test or exact test,the percentage of genes from our dataset related to the total number of genes in each pathway, the number of genes from our dataset inside the pathway and the total number of genes inside the pathway
#Output 2: NAPways corresponds to the pathway names of the pathway having less than 10 genes for the Chi-Squared or also more than 500 genes for the exact tes.
# Output 3: genesInPway correspond to each pathway with the genes from the datasets belonging to it and in which cluster they were classsify.
# Output 4: refProb is the probability of the reference in each cluster.
# Output 5: refCounts is the genes counts of the reference in each cluster.
# Output 6: pwayCounts is the genes counts of the each pathway in each cluster.
# Output 7: numOfClus is the number of clusters.
# Output 8: varStat is statistics sd, mad, cv or mean chosen for the analysis.
# Output 9: genesInClus is all the genes from the dataset and in which cluster they belong.
# Output 10: vs is the statistics value (sd, mad, cv or mean) for each gene.
########################################################
####Input
#Input 1: dat.mat is a matrix with the genes (gene symbol) on the rows and the samples on the columns.
#Input 2: pways contains the pathways of interest (KEGG, REACTOME, etc...) in the same format that makeDBList
#Input 3: test can be exact or chisq
#Input 4: varstat can be sd, mad, cv or mean.




#sd
pathVarOneSample <- function(dat.mat, pways, test = c("chisq", "exact"), varStat = c("sd", 
    "mean", "mad", "cv")) {
    varStat <- match.arg(varStat)
    test <- match.arg(test)
    # check if any GENES are in the pathway.
    geneIDs <- unique(unlist(pways$GENES))
    if (length(which(rownames(dat.mat) %in% geneIDs)) == 0) {
        stop("None of the genes in the data set are found in the given gene set or pathway")
    }
    # Compute variability for each gene
    if (varStat == "sd") {
        vs <- apply(dat.mat, 1, function(x) sd(x, na.rm = TRUE))
    } else if (varStat == "mean") {
        vs <- rowMeans(dat.mat, na.rm = TRUE)
    } else if (varStat == "mad") {
        vs <- apply(dat.mat, 1, function(x) mad(x, na.rm = TRUE))
    } else if (varStat == "cv") {
        vs <- apply(dat.mat, 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
    }
    # Forcing less than 4 clusters, otherwise the exact test can't be done because too many
    # possibilities to compute
    mm.var.filt <- suppressWarnings(Mclust(vs, G = 1:4))
    # In which cluster each gene is
    mix <- mm.var.filt$classification
    # Number of cluster
    nmix <- mm.var.filt$G
    # olap.pways contains the genes are in each pathway with their cluster number
    olap.pways <- lapply(pways$GENES, function(x) mix[names(mix) %in% x])
    names(olap.pways) <- pways$PATHNAME
    # percentage of all the genes in each cluster
    pexp <- table(mix)/length(mix)
    # list of tables of the number of genes in each cluster per pathway
    pathwayCounts <- lapply(lapply(olap.pways, function(x) table(x, deparse.level = 0)), function(x) if (length(x) != 
        nmix) {
        miss.which <- setdiff(as.character(1:nmix), names(x))
        new.obs.val <- c(rep(0, length(miss.which)), x)
        names(new.obs.val) <- c(miss.which, names(x))
        x <- as.table(new.obs.val[order(names(new.obs.val))])
    } else {
        as.table(c(x))
    })
    # Chi-Square or Exact test to compare the reference and the pathway distribution
    if (test == "chisq") {
        # chisq test and ajustment of the pvalue for each pathway
        pvals.pways <- sapply(pathwayCounts, function(x) if (sum(x) >= 10) {
            exp.val <- sum(x) * pexp  #forgot the.val
            chi <- sum((x - exp.val)^2/exp.val)
            x <- pchisq(chi, df = (length(pexp) - 1), lower.tail = FALSE)
        } else {
            NA
        })
        not_na <- which(!is.na(pvals.pways))
        pval.NA <- names(olap.pways)[-not_na]
        apvals.pways <- p.adjust(pvals.pways[not_na], "BH")
    } else if (test == "exact") {
        # Exact test and ajustment of the pvalue for each pathway
        f <- file()
        sink(file = f)
        # We perform the multinomial test on the pathway containing between 10 and 500 genes because a bigger number will involve too many possibilities to compute.
        pvals.pways <- sapply(pathwayCounts, function(x) if (sum(x) >= 10 & sum(x) < 500) {
            multinomial.test(as.vector(x), as.vector(pexp), useChisq = FALSE)$p.value
        } else {
            NA
        })
        sink()
        close(f)
        not_na <- which(!is.na(pvals.pways))
        pval.NA <- names(olap.pways)[-not_na]
        apvals.pways <- p.adjust(pvals.pways[not_na], "BH")
    }
    # Building the final table with the results
    xtab <- data.table(PwayName = pways$PATHNAME[not_na], PwayID = pways$PATHID[not_na], APval = apvals.pways, 
        PercOfGenesInPway = 100 * (sapply(olap.pways[not_na], length)/pways$SIZE[not_na]), 
        NumOfGenesFromDataSetInPathway = lengths(olap.pways[not_na]), PathwaySize = pways$SIZE[not_na])
    olap.pways <- c(olap.pways[not_na][order(xtab[, APval])], olap.pways[-not_na])
    setorder(xtab, APval, -PercOfGenesInPway, na.last = TRUE)
    rownames(xtab) <- seq(1:length(rownames(xtab)))
    out <- new("geneDistributionSet", tablePway=xtab, NAPways=pval.NA, genesInPway=olap.pways, refProb=pexp, refCounts=pexp * length(mix), pwayCounts=pathwayCounts, numOfClus=nmix, varStat=varStat, genesInClus=mix, var=vs)
    return(out)
}

#######################################################
#pathVarTwoSamplesCont
#######################################################
# 1. Compute the statistics (sd, mad, cv or mean) for each gene. 
# 2. For each pathway, we extract the gene in our dataset.
# 3. For each pathway, we look how its genes are distributed and compare the 2 groups using the bootstrap Kolmogorov-Smirnov test.
# 4. We build a dataframe for 3 named xtab (description as Output 1 below).
#######################################################
####Output
# Output 1: tablePway columns are :pathway name, pathway IDs, adjusted p-value ffrom the boot KS test, the number of genes from our dataset inside the pathway and the total number of genes inside the pathway.
#Output 2: NAPways corresponds to the pathway names of the pathway having no genes inside the dataset.
# Output 3: genesInPway correspond to the genes from the dataset belonging to each pathway
# Output 4: groups in which group each sample belongs to.
# Output 5: groupNames are the names of the two groups.
# Output 6: var1 is the statistics (sd, mad, cv or mean) for each gene for group 1.
# Output 7: var2 is the statistics (sd, mad, cv or mean) for each gene for group 2.
# Output 8: varStat is statistics sd, mad, cv or mean chosen for the analysis.
########################################################
####Input
#Input 1: dat.mat is a matrix with the genes (gene symbol) on the rows and the samples on the columns.
#Input 2: pways contains the pathways of interest (KEGG, REACTOME, etc...) in the same format that makeDBList
#Input 3: groups is a factor containing in which group each sample belongs to.
#Input 4: boot is the number of bootstraps you want to perform when applying the Kolmogorv-smirnov test.
#Input 5: varstat can be sd, mad, cv or mean.


pathVarTwoSamplesCont <- function(dat.mat, pways, groups, boot = 1000, varStat = c("sd", "mean", 
    "mad", "cv")) {
    varStat <- match.arg(varStat)
    # check if any GENES are in the pathway
    geneIDs <- unique(unlist(pways$GENES))
    if (length(which(rownames(dat.mat) %in% geneIDs)) == 0) {
        stop("None of the genes in the data set are found in the given gene set or pathway")
    }
    groupNames <- levels(groups)
    if (length(groupNames) > 2 | groupNames[1] != "1" | groupNames[2] != "2") {
        stop("Error: Only 2 groups may be compared. They must be labeled 1 and 2 in the groups parameter.")
    }
    groups <- as.factor(as.numeric(groups))
    # Compute variability for each gene and each group
    dat.mat_1 <- dat.mat[, which(groups == 1)]
    dat.mat_2 <- dat.mat[, which(groups == 2)]
    if (varStat == "sd") {
        var_1 <- apply(dat.mat_1, 1, function(x) sd(x, na.rm = TRUE))
        var_2 <- apply(dat.mat_2, 1, function(x) sd(x, na.rm = TRUE))
    } else if (varStat == "mean") {
        var_1 <- rowMeans(dat.mat_1, na.rm = TRUE)
        var_2 <- rowMeans(dat.mat_2, na.rm = TRUE)
    } else if (varStat == "mad") {
        var_1 <- apply(dat.mat_1, 1, function(x) mad(x, na.rm = TRUE))
        var_2 <- apply(dat.mat_2, 1, function(x) mad(x, na.rm = TRUE))
    } else if (varStat == "cv") {
        var_1 <- apply(dat.mat_1, 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
        var_2 <- apply(dat.mat_2, 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
    }
    genes <- row.names(dat.mat)
    # olap.pways contains the genes from the dataset in each pathway
    olap.pways <- lapply(pways$GENES, function(x) genes[genes %in% x])
    names(olap.pways) <- pways$PATHNAME
    # We compare the two densities (one for each group) of the genes of each pathway with the Kolmogorov-Smirnow test.       
    pval <- sapply(olap.pways, function(x) if (length(x) > 0) {
        ks.boot(var_1[x], var_2[x], nboots = boot)$ks.boot.pvalue
    } else {
        NA
    })
    # If the pvalue=0 we do not write 0 because it will depends on the number of bootstraps you will perform.
    pval[pval == 0] <- (1/boot) - (1/(boot)^2)
    not_na <- which(!is.na(pval))
    pval.NA <- names(pval)[-not_na]
    apvals <- p.adjust(pval[not_na], "BH")
    # We build the final table containing the results     
    xtab <- data.table(PwayName = pways$PATHNAME[not_na], PwayID = pways$PATHID[not_na], APval = apvals, 
     PercOfGenesInPway = 100 * (sapply(olap.pways[not_na], length)/pways$SIZE[not_na]),
        NumOfGenesFromDataSetInPathway = lengths(olap.pways[not_na]), PathwaySize = pways$SIZE[not_na])
    olap.pways <- c(olap.pways[not_na][order(xtab[, APval])], olap.pways[-not_na])
    setorder(xtab, APval, -PercOfGenesInPway, na.last = TRUE)
    out <- new("geneDistributionSet2", tablePway=xtab, NAPways=pval.NA, genesInPway=olap.pways, groups=groups, groupNames=groupNames, var1=var_1, var2=var_2, varStat=varStat)
    return(out)
}

#######################################################
#pathVarTwoSamplesDisc
#######################################################
# 1. Compute the statistics (sd, mad, cv or mean) for each gene.
# 2. Classify the genes with respect to the statistics in 3 clusters.
# 3. For each pathway, we extract the gene in our dataset and in which cluster they belong.
# 4. For each pathway, we look at the gene counts in each category and compare the 2 samples to each other with all the genes from the dataset with the Chi-Squared or exact test.
#######################################################
####Output
# Output 1: tablePway columns are :pathway name, pathway IDs, adjusted p-value, the percentage of genes in our dataset related to the total number of genes in each pathway, the number of genes from our dataset inside the pathway and the total number of genes inside the pathway.
#Output 2: NAPways corresponds to the pathway names of the pathway having no genes inside the dataset.
# Output 3: genesInPway1 corresponds to the genes from the dataset belonging to each pathway in the first sample
# Output 4: genesInPway2 corresponds to the genes from the dataset belonging to each pathway in the second sample
# Output 5: pwayCounts1 corresponds to a list of tables of the number of genes in each cluster per pathway for group 1
# Output 6: pwayCounts2 corresponds to a list of tables of the number of genes in each cluster per pathway for group 2
# Output 7: groups in which group each sample belongs to.
# Output 8: groupNames are the names of the two groups.
# Output 9: var1 is the statistics (sd, mad, cv or mean) for each gene for group 1.
# Output 10: var2 is the statistics (sd, mad, cv or mean) for each gene for group 2.
# Output 11: varStat is statistics sd, mad, cv or mean chosen for the analysis.
########################################################
####Input
#Input 1: dat.mat is a matrix with the genes (gene symbol) on the rows and the samples on the columns.
#Input 2: pways contains the pathways of interest (KEGG, REACTOME, etc...) in the same format that makeDBList
#Input 3: groups is a factor containing in which group each sample belongs to.
#Input 4: perc is a numeric vector of probabiliities with values between 0 and 1. Used to put genes into clusters
#Input 5: test is a string, either "exact" or "chisq", which are tests to see if clusters in the 2 samples are sig. different from each other
#Input 6: varstat can be sd, mad, cv or mean.


pathVarTwoSamplesDisc <- function(dat.mat, pways, groups, perc = c(1/3, 2/3), test = c("chisq", 
    "exact"), varStat = c("sd", "mean", "mad", "cv")) {
    varStat <- match.arg(varStat)
    test <- match.arg(test)
    # check if any GENES are in the pathway
    geneIDs <- unique(unlist(pways$GENES))
    if (length(which(rownames(dat.mat) %in% geneIDs)) == 0) {
        stop("None of the genes in the data set are found in the given gene set or pathway")
    }
    groupNames <- levels(groups)
    if (length(groupNames) > 2 | groupNames[1] != "1" | groupNames[2] != "2") {
        stop("Error: Only 2 groups may be compared. They must be labeled 1 and 2 in the groups parameter.")
    }
    groups <- as.factor(as.numeric(groups))
    # Compute variability for each gene and each group
    dat.mat_1 <- dat.mat[, which(groups == 1)]
    dat.mat_2 <- dat.mat[, which(groups == 2)]
    if (varStat == "sd") {
        vs <- apply(dat.mat, 1, function(x) sd(x, na.rm = TRUE))
        var_1 <- apply(dat.mat_1, 1, function(x) sd(x, na.rm = TRUE))
        var_2 <- apply(dat.mat_2, 1, function(x) sd(x, na.rm = TRUE))
    } else if (varStat == "mean") {
        vs <- apply(dat.mat, 1, function(x) mean(x, na.rm = TRUE))
        var_1 <- rowMeans(dat.mat_1, na.rm = TRUE)
        var_2 <- rowMeans(dat.mat_2, na.rm = TRUE)
    } else if (varStat == "mad") {
        vs <- apply(dat.mat, 1, function(x) mad(x, na.rm = TRUE))
        var_1 <- apply(dat.mat_1, 1, function(x) mad(x, na.rm = TRUE))
        var_2 <- apply(dat.mat_2, 1, function(x) mad(x, na.rm = TRUE))
    } else if (varStat == "cv") {
        vs <- apply(dat.mat, 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
        var_1 <- apply(dat.mat_1, 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
        var_2 <- apply(dat.mat_2, 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
    }
    cut <- quantile(vs, probs = perc)
    genes <- row.names(dat.mat)
    # In which cluster each gene is
    mix1 <- rep(NA, length(genes))
    mix1[var_1 <= cut[1]] <- 1
    mix1[var_1 <= cut[2] & var_1 > cut[1]] <- 2
    mix1[var_1 > cut[2]] <- 3
    names(mix1) <- genes
    # In which cluster each gene is
    mix2 <- rep(NA, length(genes))
    mix2[var_2 <= cut[1]] <- 1
    mix2[var_2 <= cut[2] & var_2 > cut[1]] <- 2
    mix2[var_2 > cut[2]] <- 3
    names(mix2) <- genes
    # olap.pways contains the genes from the dataset in each pathway
    olap.pways1 <- lapply(pways$GENES, function(x) mix1[names(mix1) %in% x])
    names(olap.pways1) <- pways$PATHNAME
    olap.pways2 <- lapply(pways$GENES, function(x) mix2[names(mix2) %in% x])
    names(olap.pways2) <- pways$PATHNAME
    # list of tables of the number of genes in each cluster per pathway
    pathwayCounts1 <- lapply(lapply(olap.pways1, function(x) table(x, deparse.level = 0)), 
        function(x) if (length(x) != 3) {
            miss.which <- setdiff(as.character(1:3), names(x))
            new.obs.val <- c(rep(0, length(miss.which)), x)
            names(new.obs.val) <- c(miss.which, names(x))
            x <- as.table(new.obs.val[order(names(new.obs.val))])
        } else {
            as.table(c(x))
        })
    # list of tables of the number of genes in each cluster per pathway
    pathwayCounts2 <- lapply(lapply(olap.pways2, function(x) table(x, deparse.level = 0)), 
        function(x) if (length(x) != 3) {
            miss.which <- setdiff(as.character(1:3), names(x))
            new.obs.val <- c(rep(0, length(miss.which)), x)
            names(new.obs.val) <- c(miss.which, names(x))
            x <- as.table(new.obs.val[order(names(new.obs.val))])
        } else {
            as.table(c(x))
        })
    if (test == "chisq") {
        # chisq test and ajustment of the pvalue for each pathway
        pvals.pways <- sapply(pways$PATHNAME, function(x) if (sum(pathwayCounts1[x][[1]]) >= 
            10) {
            exp.val <- pathwayCounts1[x][[1]]  #forgot the.val
            chi <- sum((pathwayCounts2[x][[1]] - exp.val)^2/exp.val)
            x <- pchisq(chi, df = (3 - 1), lower.tail = FALSE)
        } else {
            NA
        })
        not_na <- which(!is.na(pvals.pways))
        pval.NA <- pways$PATHNAME[-not_na]
        apvals.pways <- p.adjust(pvals.pways[not_na], "BH")
    } else if (test == "exact") {
        # Exact test and ajustment of the pvalue for each pathway
        f <- file()
        sink(file = f)
        # We perform the multinomial test on the pathway containing between 10 and 500 genes because a bigger number will involve too many possibilities to compute.
        pvals.pways <- sapply(pways$PATHNAME, function(x) if (sum(pathwayCounts1[x][[1]]) >= 
            10 & sum(pathwayCounts1[x][[1]]) < 500) {
            pexp <- pathwayCounts1[x][[1]]/sum(pathwayCounts1[x][[1]])
            multinomial.test(as.vector(pathwayCounts2[x][[1]]), as.vector(pexp), useChisq = FALSE)$p.value
        } else {
            NA
        })
        sink()
        close(f)
        not_na <- which(!is.na(pvals.pways))
        pval.NA <- pways$PATHNAME[-not_na]
        apvals.pways <- p.adjust(pvals.pways[not_na], "BH")
    }
    # We build the final table containing the results     
    xtab <- data.table(PwayName = pways$PATHNAME[not_na], PwayID = pways$PATHID[not_na], APval = apvals.pways, 
        PercOfGenesInPway = 100 * (sapply(olap.pways1[not_na], length)/pways$SIZE[not_na]), 
        NumOfGenesFromDataSetInPway = lengths(olap.pways1[not_na]), PathwaySize = pways$SIZE[not_na])
    olap.pways1 <- c(olap.pways1[not_na][order(xtab[, APval])], olap.pways1[-not_na])
    olap.pways2 <- c(olap.pways2[not_na][order(xtab[, APval])], olap.pways2[-not_na])
    setorder(xtab, APval, -PercOfGenesInPway, na.last = TRUE)
    rownames(xtab) <- seq(1:length(rownames(xtab)))
    out <- new("geneDistributionSet3", tablePway=xtab, NAPways=pval.NA, genesInPway1=olap.pways1, genesInPway2=olap.pways2, pwayCounts1=pathwayCounts1, pwayCounts2=pathwayCounts2, groups=groups, groupNames=groupNames, var1=var_1, var2=var_2, varStat=varStat)
    return(out)
}

#######################################################
#sigOneSample
#######################################################
#It is a function that returns the significant pathway(s),which category(ies) from this pathway are significant and which gene(s) belongs to this(ese) category(ies).

#######################################################
####Output
# Output 1: genesInSigPways1 contains the genes per significant pathway belonging to the significant category.
#Output 2: sigCatPerPway contains the category(ies) per pathway that are significant.
# Output 3: thresPValue is the chosen p-value for the significance
########################################################
####Input
#Input 1: pvalue_results is result from the pathVarOneSample function
#Input 2: pvalue is the significance we want to test.



sigOneSample <- function(pvalue_results, pvalue) {
    xtab <- pvalue_results@tablePway
    nmix <- pvalue_results@numOfClus
    olap.pways <- pvalue_results@genesInPway
    pexp <- pvalue_results@refProb
    pwayCounts <- pvalue_results@pwayCounts
    if (length(rownames(xtab[APval < pvalue, ])) < 1) {
        warning("There are no significant pathways. Quitting significant_category function and returning empty object")
        sig <- new("significantPathway", genesInSigPways1=list(), sigCatPerPway=list(), thresPValue=numeric())
        return(sig)
    }
    # PathName that were significant in xtab.
    significant <- xtab[APval < pvalue, PwayName]
    # The list of table with the number of genes in each cluster from the significant pathways
    sigPwayCounts <- pwayCounts[significant]
    genes <- vector("list", length(significant))
    category <- vector("list", length(significant))
    names(genes) <- significant
    names(category) <- significant
    # results contain the p-value for each category in each pathway computed with the binomial test.
    f <- file()
    sink(file = f)
    results <- sapply(sigPwayCounts, function(x) apply(rbind(x, pexp), 2, function(y) multinomial.test(c(y[1], 
        sum(x) - y[1]), prob = c(y[2], 1 - y[2]))$p.value))
    sink()
    close(f)
    row.names(results) <- NULL
    # For each significant pathway we look which category(ies) is are significant and the genes
    # belonging to this(ese) category(ies). If no category is significant we will write 0 in
    # the category and the gene names.
    genes <- olap.pways[significant]
    sigCat <- which(apply(results, 2, function(x) sum(x < pvalue)) > 0)
    numOfSigCat <- apply(results[, sigCat], 2, function(x) length(which(x < pvalue)))
    if (length(which(numOfSigCat == nmix)) == length(sigCat)) {
        category[sigCat] <- lapply(as.list(data.frame(apply(results[, sigCat], 2, function(x) which(x < 
            pvalue)))), function(x) unique(x))
    } else {
        category[sigCat] <- apply(results[, sigCat], 2, function(x) which(x < pvalue))
    }
    noSigCat <- which(apply(results, 2, function(x) sum(x < pvalue)) < 1)
    if (length(noSigCat) > 0) {
        category[noSigCat] <- 0
    }
    sig <- new("significantPathway", genesInSigPways1=genes, sigCatPerPway=category, thresPValue=pvalue)
    return(sig)
}
#######################################################
#sigTwoSamplesCont
#######################################################
#It is a function that returns the significant pathways and which genes belongs to these #pathways.

#######################################################
####Output
# Output 1: genesInSigPways1 contains the genes belonging to each significant pathway
# Output 2: thresPValue is the chosen p-value for the significance
########################################################
####Input
#Input 1: pvalue_results is result from the pathVarTwoSamplesCont function
#Input 2: pvalue is the significance we want to test.


sigTwoSamplesCont <- function(pvalue_results, pvalue) {
    xtab <- pvalue_results@tablePway
    olap.pways <- pvalue_results@genesInPway
    if (length(rownames(xtab[APval < pvalue, ])) < 1) {
        warning("There are no significant pathways. Quitting significant_category function and returning empty object")
        sig <- new("significantPathway2", genesInSigPways1=list(), thresPValue=numeric())
        return(sig)
    }
    # Pathways that were significant in xtab.
    significant <- xtab[APval < pvalue, PwayName]
    # Genes from the dataset inside each significant pathway
    genes <- olap.pways[significant]
    sig <- new("significantPathway2", genesInSigPways1=genes, thresPValue=pvalue)
    return(sig)
}
#######################################################
#sigTwoSamplesDisc
#######################################################
#It is a function that returns the significant pathways and which genes belong to these pathways

#######################################################
####Output
# Output 1: genesInSigPways1 contains the genes belonging to each significant pathway in significant categories in the first sample
# Output 2: genesInSigPways2 contains the genes belonging to each significant pathway in significant categories in the second sample
# Output 3: sigCatPerPway contains the significant categories in each pathway
# Output 4: thresPValue is the chosen p-value for the significance
########################################################
####Input
#Input 1: pvalue_results is result from the pathVarTwoSamplesDisc function
#Input 2: pvalue is the significance we want to test.

sigTwoSamplesDisc <- function(pvalue_results, pvalue) {
    xtab <- pvalue_results@tablePway
    olap.pways1 <- pvalue_results@genesInPway1
    olap.pways2 <- pvalue_results@genesInPway2
    pwayCounts1 <- pvalue_results@pwayCounts1
    pwayCounts2 <- pvalue_results@pwayCounts2
    if (length(rownames(xtab[APval < pvalue, ])) < 1) {
        warning("There are no significant pathways. Quitting significant_category function and returning empty object")
        sig <- new("significantPathway3", genesInSigPways1=list(), genesInSigPways2=list(), sigCatPerPway=list(), thresPValue=numeric())
        return(sig)
    }
    # PathName that were significant in xtab.
    significant <- xtab[APval < pvalue, PwayName]
    # The list of table with the number of genes in each cluster from the significant pathways
    sigPwayCounts1 <- pwayCounts1[significant]
    sigPwayCounts2 <- pwayCounts2[significant]
    genes <- vector("list", length(significant))
    category <- vector("list", length(significant))
    names(genes) <- significant
    names(category) <- significant
    # results contain the p-value for each category in each pathway computed with the binomial
    # test
    f <- file()
    sink(file = f)
    results <- sapply(significant, function(x) apply(rbind(sigPwayCounts2[x][[1]], sigPwayCounts1[x][[1]]/sum(sigPwayCounts1[x][[1]])), 
        2, function(y) multinomial.test(c(y[1], sum(sigPwayCounts2[x][[1]]) - y[1]), prob = c(y[2], 
            1 - y[2]))$p.value))
    sink()
    close(f)
    row.names(results) <- NULL
    # For each significant pathway we look which category(ies) is are significant and the genes belonging to this(ese) category(ies). If no category is significant we will write 0 in the category and the gene names.
    genes1 <- olap.pways1[significant]
    genes2 <- olap.pways2[significant]
    sigCat <- which(apply(results, 2, function(x) sum(x < pvalue)) > 0)
    numOfSigCat <- apply(results[, sigCat], 2, function(x) length(which(x < pvalue)))
    if (length(which(numOfSigCat == 3)) == length(sigCat)) {
        category[sigCat] <- lapply(as.list(data.frame(apply(results[, sigCat], 2, function(x) which(x < 
            pvalue)))), function(x) unique(x))
    } else {
        category[sigCat] <- apply(results[, sigCat], 2, function(x) which(x < pvalue))
    }
    noSigCat <- which(apply(results, 2, function(x) sum(x < pvalue)) < 1)
    if (length(noSigCat) > 0) {
        category[noSigCat] <- 0
    }
    sig <- new("significantPathway3", genesInSigPways1=genes1, genesInSigPways2=genes2, sigCatPerPway=category, thresPValue=pvalue)
    return(sig)
}
#######################################################
#sigPway
#######################################################
#It is a function that check if an object is from the one sample or two samples cases and then use sigOneSample, sigTwoSamplesCont, or sigTwoSamplesDisc to find the significant pathways.
#######################################################
####Output
# the result of sigOneSample, sigTwoSamplesCont, or sigTwoSamplesDisc
########################################################
####Input
#Input 1: pvalue_results is result from the pathVarOneSample, pathVarTwoSamplesCont, or pathVarTwoSamplesDisc function.
#Input 2: pvalue is the significance we want to test.
########################################################


sigPway <- function(pvalue_results, pvalue) {
    # Check whether we are in the one sample or two samples case
    if (class(pvalue_results) == "geneDistributionSet") {
        sigOneSample(pvalue_results, pvalue)
    } else if (class(pvalue_results) == "geneDistributionSet2") {
        sigTwoSamplesCont(pvalue_results, pvalue)
    } else {
        sigTwoSamplesDisc(pvalue_results, pvalue)
    }
}


#######################################################
#plotOneSample
#######################################################
#It is a function that returns the plot of the reference counts along with the plot of a chosen #pathway. This function is made for output from pathVarOneSample.
#######################################################
####Output
# plot of the reference and a pathway counts
########################################################
####Input
#Input 1: pvalue_results is result from the pathVarOneSample function
#Input 2: pathway is the chosen pathway you want to plot.
#Input 3: sig should be an output of sigOneSample or NULL.
########################################################
####Remark
#If sig is not NULL, the function will check if the pathway is a significant one and if yes the
#title will be printed in red and the significant category(ies) will be highlighted by a red
#star.



plotOneSample <- function(pvalue_results, pathway, sig) {
    mp <- pathway
    # If the name of the pathway is two long it will cut it into two lines in the plot.
    if (nchar(mp) > 39) {
        mp <- unlist(strsplit(mp, " "))
        if (is.na((which(cumsum(nchar(mp)) > 34)[1]))) {
            l <- length(mp) - 1
        } else {
            l <- which(cumsum(nchar(mp)) > 34)[1] - 1
        }
        mp_1 <- paste(mp[1:l], collapse = " ")
        mp_2 <- paste(mp[(l + 1):length(mp)], collapse = " ")
        mp <- paste(mp_1, "\n", mp_2)
    }
    ref <- pvalue_results@refCounts
    path <- pvalue_results@pwayCounts[[pathway]]
    # if we have the result from the sigOneSample, it will be included in the plot
    if (!is.null(sig)) {
        category <- sig@sigCatPerPway
    } else {
        category <- NULL
    }
    # data frame for the reference counts
    ref <- as.data.frame(ref)
    colnames(ref) <- c("Cluster", "Number_of_genes")
    # data frame for the pathway distribution
    path <- as.data.frame(path)
    colnames(path) <- c("Cluster", "Number_of_genes")
    # plot of the reference counts
    plotRef <- ggplot(ref, aes(x = Cluster, y = Number_of_genes, fill = Cluster)) + geom_bar(stat = "identity", 
        color = "black", fill = "grey") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Number of genes") + 
        theme(legend.position = "none") + ggtitle("Reference") + xlab("")
    # Plot of the pwathway counts
    d <- ggplot(path, aes(x = Cluster, y = Number_of_genes, fill = Cluster)) + geom_bar(stat = "identity", 
        color = "black", fill = "steelblue1") + ylab("Number of genes") + theme(legend.position = "none") + 
        ggtitle(mp) + xlab("")
    # If the pathway is one of the significant ones, the title will be in red. and the categories, if any, we be highlighted with a red star.
    if (pathway %in% names(category)) {
        sigCat <- category[[pathway]]
        plotPathway <- d + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            plot.title = element_text(colour = "red"))
        if (sum(sigCat) > 0) {
            plotPathway <- plotPathway + annotate("text", x = sigCat, y = path[sigCat + 0.1, 
                2], label = rep("*", length(sigCat)), color = "red", size = 15)
        }
    } else {
        plotPathway <- d + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    }
    # plot the reference and pathway counts side by side
    grid.arrange(arrangeGrob(plotRef, plotPathway, nrow = 1))
}
#######################################################
#plotTwoSamplesCont
#######################################################
#It is a function that returns the plot of the two densities (one for each group) of the statistics (sd, mad, cv or mean) of a chosen pathway. This function is made for output from pathVarTwoSamplesCont.
#######################################################
####Output
# plot of the two densities of the statistics.
########################################################
####Input
#Input 1: pvalue_results is result from the pathVarTwoSamplesCont function
#Input 2: pathway is the chosen pathway you want to plot.
#Input 3: sig should be an output of sigTwoSamplesCont or NULL.
########################################################
####Remark
#If sig is not NULL, the function will check if the pathway is a significant one and if yes the title will be printed in red.


plotTwoSamplesCont <- function(pvalue_results, pathway, sig) {
    mp <- pathway
    # If the name of the pathway is two long it will cut it into two lines in the plot.
    if (nchar(mp) > 39) {
        mp <- unlist(strsplit(mp, " "))
        if (is.na((which(cumsum(nchar(mp)) > 34)[1]))) {
            l <- length(mp) - 1
        } else {
            l <- which(cumsum(nchar(mp)) > 34)[1] - 1
        }
        mp_1 <- paste(mp[1:l], collapse = " ")
        mp_2 <- paste(mp[(l + 1):length(mp)], collapse = " ")
        mp <- paste(mp_1, "\n", mp_2)
    }
    xtab <- pvalue_results@tablePway
    # If the number of genes of the pathway is less than 3, it is not possible to draw a density and it will return an empty plot with this message.
    if (xtab[PwayName == pathway, NumOfGenesFromDataSetInPathway] < 3) {
        df <- data.frame()
        plotPway <- ggplot(df) + geom_point() + xlim(0, 3) + ylim(0, 1) + theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
            ylab("Density") + theme(legend.position = "none") + xlab(varStat) + ggtitle(mp) + 
            geom_text(aes(1.5, 0.5, label = "Less than 3 genes, not possible to draw a density.", 
                col = "red"))
    } else {
        var_1 <- pvalue_results@var1
        var_2 <- pvalue_results@var2
        varStat <- pvalue_results@varStat
        genes <- pvalue_results@genesInPway[[pathway]]
        grp_1 <- as.vector(var_1[genes])
        grp_2 <- as.vector(var_2[genes])
        df <- data.frame(cbind(c(grp_1, grp_2), c(rep("Group 1", length(grp_1)), rep("Group 2", 
            length(grp_2)))))
        df.m <- melt(df)
        df.m[, 1] <- as.numeric(as.character(df.m[, 1]))
        colnames(df.m) <- c("value", "group")
        color <- c("steelblue1", "grey")
        # Plot of the two densities (one for each group) of the variability of the genes inside the pathway.
        d <- ggplot(na.omit(df.m)) + geom_density(aes(x = value, y = ..density.., colour = group, 
            fill = group), alpha = 0.3) + scale_colour_manual(values = color) + scale_fill_manual(values = color) + 
            theme(legend.position = "bottom", legend.title = element_blank()) + ylab("Density") + 
            xlab(varStat) + ggtitle(mp)
        if (!is.null(sig)) {
            significant <- names(sig@genesInSigPways1)
        } else {
            significant <- NULL
        }
        # If we included the results of sigTwoSamplesCont, it will verify if the pathway is one of them and if yes the title will be printed in red.
        if (pathway %in% significant) {
            plotPathway <- d + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                plot.title = element_text(colour = "red"))
        } else {
            plotPathway <- d + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
        }
    }
    plot(plotPathway)
}

#######################################################
#plotTwoSamplesDisc
#######################################################
##It is a function that returns 2 plots of the 2 samples for a chosen pathway. This function is made for output from pathVarTwoSamplesDisc.
#######################################################
####Output
# plot of the 2 samples for a significant pathway
########################################################
####Input
#Input 1: pvalue_results is result from the pathVarTwoSamplesDisc function
#Input 2: pathway is the chosen pathway you want to plot.
#Input 3: sig should be an output of sigTwoSamplesDisc or NULL.
########################################################
####Remark
#If sig is not NULL, the function will check if the pathway is a significant one and if yes the title will be printed in red.


plotTwoSamplesDisc <- function(pvalue_results, pathway, sig) {
    mp <- pathway
    # If the name of the pathway is two long it will cut it into two lines in the plot.
    if (nchar(mp) > 39) {
        mp <- unlist(strsplit(mp, " "))
        if (is.na((which(cumsum(nchar(mp)) > 34)[1]))) {
            l <- length(mp) - 1
        } else {
            l <- which(cumsum(nchar(mp)) > 34)[1] - 1
        }
        mp_1 <- paste(mp[1:l], collapse = " ")
        mp_2 <- paste(mp[(l + 1):length(mp)], collapse = " ")
        mp <- paste(mp_1, "\n", mp_2)
    }
    path1 <- pvalue_results@pwayCounts1[[pathway]]
    path2 <- pvalue_results@pwayCounts2[[pathway]]
    # if we have the result from the sigOneSample, it will be included in the plot
    if (!is.null(sig)) {
        category <- sig@sigCatPerPway
    } else {
        category <- NULL
    }
    # data frame for the reference counts
    path1 <- as.data.frame(path1)
    colnames(path1) <- c("Cluster", "Number_of_genes")
    # data frame for the pathway distribution
    path2 <- as.data.frame(path2)
    colnames(path2) <- c("Cluster", "Number_of_genes")
    yLimMax <- max(path1[, 2], path2[, 2])
    # plot of the reference counts
    plotPath1 <- ggplot(path1, aes(x = Cluster, y = Number_of_genes, fill = Cluster)) + geom_bar(stat = "identity", 
        color = "black", fill = "darkorchid1") + ylab("Number of genes") + theme(legend.position = "none") + 
        ggtitle(paste(mp, "\n Group 1")) + xlab("") + ylim(0, yLimMax)
    # Plot of the pwathway counts
    plotPath2 <- ggplot(path2, aes(x = Cluster, y = Number_of_genes, fill = Cluster)) + geom_bar(stat = "identity", 
        color = "black", fill = "steelblue1") + ylab("Number of genes") + theme(legend.position = "none") + 
        ggtitle(paste(mp, "\n Group 2")) + xlab("") + ylim(0, yLimMax)
    # If the pathway is one of the significant ones, the title will be in red. and the categories, if any, we be highlighted with a red star.
    if (pathway %in% names(category)) {
        sigCat <- category[[pathway]]
        plotPathway1 <- plotPath1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            plot.title = element_text(colour = "red"))
        plotPathway2 <- plotPath2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            plot.title = element_text(colour = "red"))
        if (sum(sigCat) > 0) {
            plotPathway1 <- plotPathway1 + annotate("text", x = sigCat, y = path1[sigCat + 
                0.1, 2], label = rep("*", length(sigCat)), color = "red", size = 15)
            plotPathway2 <- plotPathway2 + annotate("text", x = sigCat, y = path2[sigCat + 
                0.1, 2], label = rep("*", length(sigCat)), color = "red", size = 15)
        }
    } else {
        plotPathway1 <- plotPath1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
        plotPathway2 <- plotPath2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    }
    # plot the reference and pathway counts side by side
    grid.arrange(arrangeGrob(plotPathway1, plotPathway2, nrow = 1))
}

#######################################################
#plotPway
#######################################################
#It is a function that check if an object is from the one sample or two samples cases and then use plotOneSample or plotTwoSamples for the chosen pathway.
#######################################################
####Output
# plot of the results of the one or two samples case for a chosen pathway.
########################################################
####Input
#Input 1: pvalue_results is the result from the pathVarOneSample, pathVarTwoSamplesCont, or pathVarTwoSamplesDisc function
#Input 2: pathway is the chosen pathway you want to plot.
#Input 3: sig should be an output of sigOneSample, sigTwoSamplesCont, sigTwoSamplesDisc or NULL.
########################################################
####Remark
#If sig is not NULL, the function will check if the pathway is a significant one. And they will be highlighted in the resulting plot (see plotOneSample or plotTwosamples)


plotPway <- function(pvalue_results, pathway, sig = NULL) {
    # Check whether we are in the one sample or two samples case
    if (class(pvalue_results) == "geneDistributionSet") {
        plotOneSample(pvalue_results, pathway, sig)
    } else if (class(pvalue_results) == "geneDistributionSet2") {
        plotTwoSamplesCont(pvalue_results, pathway, sig)
    } else {
        plotTwoSamplesDisc(pvalue_results, pathway, sig)
    }
}

#######################################################
#saveAsPDF
#######################################################
#Save as a pdf the plots for the one or two samples case of the significant pathway or a chosen list of pathway..
########################################################
####Output
# Save as a pdf the plots of the significant pathway or a chosen list of pathway.
########################################################
####Input
#Input 1: pvalue_results is the result from the pathVarOneSample, pathVarTwoSamplesCont, or pathVarTwoSamplesDisc function
#Input 2: sig should be an output of sigOneSample, sigTwoSamplesCont, sigTwoSamplesDisc, or NULL.
#Input 3: listPath is "significant" if you want to save the plots of the significant pathways or can be a list of names of pathway of interest.
########################################################
####Remark
#If sig is not NULL, the function will check if the pathway is a significant one. And they will be highlighted in the resulting plot (see plotOneSample or plotTwosamples)


saveAsPDF <- function(pvalue_results, sig, listPath = "significant") {
    # If listPath='significant' we will save as pdf all the plots corresponding to the significant pathway from sig. Other wise it will save the pathways given to listPath.
    if (listPath[1] == "significant") {
        listPath <- names(sig@genesInSigPways1)
    }
    # The name of the file will be the pathname where we replace '/' by '_'
    pathname <- sapply(listPath, function(x) if (length(unlist(strsplit(x, "/"))) > 1) {
        x <- paste(unlist(strsplit(x, "/")), collapse = "_")
    } else {
        x
    })
    # save as PDF all the pathways significant or given in listPath
    for (i in 1:length(pathname)) {
        pdf(file = paste(pathname[i], ".pdf", sep = ""), width = 10, height = 7)
        plotPway(pvalue_results, listPath[i], sig)
        dev.off()
    }
}

#######################################################
#getGenes
#######################################################
#It is a function that returns one list of genes for group 1 and one for group 2 of a chosen pathway having their statistics (sd, mad, cv or mean) inside a chosen interval.
#######################################################
####Output
# Output 1: genes1 contains the genes belonging to the pathway in the given window for group 1.
# Output 2: genes2 contains the genes belonging to the pathway in the given window for group 2.
# Output 3: genesAll contains the genes from the dataset belonging to the pathway
########################################################
####Input
#Input 1: pvalue_results is result from the pathVarTwoSamplesCont function
#Input 2: pathway is the chosen pathway.
#Input 3: window is the chosen interval.



getGenes <- function(pvalue_results, pathway, window) {
    gene_window <- NULL
    olap.pways <- pvalue_results@genesInPway
    var_1 <- pvalue_results@var1
    var_2 <- pvalue_results@var2
    genes <- olap.pways[[pathway]]
    grp_1 <- as.vector(var_1[genes])
    grp_2 <- as.vector(var_2[genes])
    # Take the genes from group 1 from the pathway belonging to the window
    genes_1 <- genes[which(grp_1 >= window[1] & grp_1 <= window[2])]
    # Take the genes from group 3 from the pathway belonging to the window
    genes_2 <- genes[which(grp_2 >= window[1] & grp_2 <= window[2])]
    # Take all the genes from the pathway
    genes_all <- genes
    gene_window <- new("geneSet", genes1=genes_1, genes2=genes_2, genesAll=genes_all)
    return(gene_window)
}
#######################################################
#diagnosticsVarPlots
#######################################################
#It is a function that returns a graph with three plots. One with mean vs sd, one with mean vs mad and one with mean vs cv. It also print in the title the correlation between each variability statistics and the mean. It also draws the loess curve.
#######################################################
####Output
# graph with the three plots of variability vs mean.
########################################################
####Input
#Input 1: dat.mat is a matrix with the genes (gene symbol) on the rows and the samples on the columns.


diagnosticsVarPlots <- function(dat.mat) {
    # Compute the sd, mean,cv and mad for each gene
    dat.sd <- apply(dat.mat, 1, sd, na.rm = TRUE)
    dat.avg <- rowMeans(dat.mat, na.rm = TRUE)
    dat.cv <- apply(dat.mat, 1, function(x) {
        sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
    })
    dat.mad <- apply(dat.mat, 1, mad, na.rm = TRUE)
    data.sd <- data.frame(standDev = dat.sd, avg = dat.avg)
    data.mad <- data.frame(medAbsDev = dat.mad, avg = dat.avg)
    data.cv <- data.frame(cv = dat.cv, avg = dat.avg)
    # Plot the mean vs the standard deviation and give the correlation between both in the title.
    plotSD <- ggplot(data.sd, aes(x = avg, y = standDev)) + geom_point(colour = "grey60") + 
        stat_smooth(method = loess) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Average") + 
        ylab("Standard Deviation (SD)") + ggtitle(paste(c("Standard Deviation [ R=", round(cor(dat.sd, 
        dat.avg), 3), "]"), collapse = " "))
    # Plot the mean vs the median absolute deviation and give the correlation between both in the title.
    plotMAD <- ggplot(data.mad, aes(x = avg, y = medAbsDev)) + geom_point(colour = "grey60") + 
        stat_smooth(method = loess) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Average") + 
        ylab("Median Absolute Deviation (MAD)") + ggtitle(paste(c("Median Absolute Deviation [ R=", 
        round(cor(dat.mad, dat.avg), 3), "]"), collapse = " "))
    # Plot the mean vs the coefficient of variation and give the correlation between both in the title.
    plotCV <- ggplot(data.cv, aes(x = avg, y = cv)) + geom_point(colour = "grey60") + stat_smooth(method = loess) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) + xlab("Average") + ylab("Coefficient of Variation (CV)") + 
        ggtitle(paste(c("Coefficient of Variation [ R=", round(cor(dat.cv, dat.avg), 3), "]"), 
            collapse = " "))
    grid.arrange(arrangeGrob(plotSD, plotMAD, plotCV, nrow = 1))
}

