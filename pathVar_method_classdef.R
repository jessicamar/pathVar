setClass("geneDistributionSet",representation(tablePway="data.table",NAPways="character",genesInPway="list",refProb="table",refCounts="table",pwayCounts="list",numOfClus="numeric",varStat="character",genesInClus="numeric",var="numeric"))


setClass("geneDistributionSet2",representation(tablePway="data.table",NAPways="character",genesInPway="list",groups="factor",groupNames="character",var1="numeric",var2="numeric",varStat="character"))

setClass("geneDistributionSet3",representation(tablePway="data.table",NAPways="character",genesInPway1="list",genesInPway2="list",pwayCounts1="list",pwayCounts2="list",groups="factor",groupNames="character",var1="numeric",var2="numeric",varStat="character"))

setClass("significantPathway",representation(genesInSigPways1="list",sigCatPerPway="list",thresPValue="numeric"))

setClass("significantPathway2",representation(genesInSigPways1="list",thresPValue="numeric"))

setClass("significantPathway3",representation(genesInSigPways1="list",genesInSigPways2="list",sigCatPerPway="list",thresPValue="numeric"))



setClass("geneSet",representation(genes1="character",genes2="character",genesAll="character"))