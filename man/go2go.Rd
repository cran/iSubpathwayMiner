\name{go2go}
\alias{getKGeneFromGene}
\alias{getGeneFromKGene}
\alias{getEnzymeFromKGene}
\alias{getKGeneFromEnzyme}
\alias{getKOFromKGene}
\alias{getKGeneFromKO}
\alias{getPathwayFromKGene}
\alias{getKGeneFromPathway}
\alias{getKGeneFromSymbol}
\alias{getSymbolFromKGene}
\alias{getOrgIdFromOrgName}
\alias{getOrgNameFromOrgId}

\alias{getEnzymeFromGene}
\alias{getGeneFromEnzyme}
\alias{getKOFromGene}
\alias{getGeneFromKO}
\alias{getPathwayFromGene}
\alias{getGeneFromPathway}
\alias{getSymbolFromGene}
\alias{getGeneFromSymbol}
\alias{getEnzymeFromSymbol}
\alias{getSymbolFromEnzyme}
\alias{getKOFromSymbol}
\alias{getSymbolFromKO}

\title{Map between identifiers}
\description{
  Map between two types of identifiers. For example, gene identifiers are converted to enzyme identifiers.
}
\usage{
getKGeneFromGene(geneList)
getGeneFromKGene(keggGeneList)
getEnzymeFromKGene(keggGeneList,ignoreAmbiguousEnzyme=TRUE)
getKGeneFromEnzyme(enzymeList,ignoreAmbiguousEnzyme=TRUE)
getKOFromKGene(keggGeneList)
getKGeneFromKO(KOList)
getPathwayFromKGene(keggGeneList)
getKGeneFromPathway(pathwayList)
getKGeneFromSymbol(symbolList)
getSymbolFromKGene(keggGeneList)
getOrgIdFromOrgName(orgName)
getOrgNameFromOrgId(orgId)

getEnzymeFromGene(geneList,ignoreAmbiguousEnzyme=TRUE)
getGeneFromEnzyme(enzymeList,ignoreAmbiguousEnzyme=TRUE)
getKOFromGene(geneList)
getGeneFromKO(KOList)
getPathwayFromGene(geneList)
getGeneFromPathway(pathwayList)
getSymbolFromGene(geneList)
getGeneFromSymbol(symbolList)
getEnzymeFromSymbol(symbolList)
getSymbolFromEnzyme(enzymeList)
getKOFromSymbol(symbolList)
getSymbolFromKO(KOList)
}
\arguments{
  \item{geneList}{A character vector of genes that are the current gene identifiers.}
  \item{keggGeneList}{A character vector of KEGG genes.}
  \item{enzymeList}{A character vector of enzymes.}
  \item{KOList}{A character vector of KOs.}
  
  \item{pathwayList}{A character vector of KEGG pathway identifiers.}
  \item{symbolList}{A character vector of official symbols.}
  \item{orgId}{A character vector of the three- or four-letter KEGG organism codes.}
  \item{orgName}{A character vector of organism names.}
  \item{ignoreAmbiguousEnzyme}{A logical. Should ambiguous enzymes such as 2.6.1.- be ignored?}  
}
\details{
Note that the argument \code{geneList} should be the current gene identifiers. The function \code{\link{getOrgAndIdType}} can get the current the type of gene identifier. The default is ncbi-geneid. The argument \code{keggGeneList} should be KEGG gene identifiers. The argument \code{symbolList} should be gene official symbols.

Note that the result is the union of sets of identifiers.

Note that the argument \code{ignoreAmbiguousEnzyme} is only applied to pathway graphs with enzymes as nodes. If TRUE, genes that ambiguous enzymes (e.g., "-") correspond to are ignored in analyses.

\preformatted{
   ##ignoreAmbiguousEnzyme=TRUE
   getGeneFromEnzyme(c("ec:2.6.1.-"),TRUE)
   ##ignoreAmbiguousEnzyme=FALSE
   getGeneFromEnzyme(c("ec:2.6.1.-"),FALSE)
}

}
\value{
   A character vector of identifiers.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\examples{
\dontrun{

## see the current gene identifier type
getOrgAndIdType()[2] 



## get the list of KEGG gene identifiers of three gene (the current identifier: ncbi-geneid) identifiers.
getKGeneFromGene(c("1","5232","5224"))

## get the list of gene identifiers of two KEGG gene identifiers.
getGeneFromKGene(c("hsa:5232","hsa:5224"))

## get the list of enzymes of two KEGG genes.
getEnzymeFromKGene(c("hsa:5232","hsa:5224"))

## get the list of KEGG genes of two enzymes. 
getKGeneFromEnzyme(c("ec:2.7.2.3","ec:5.4.2.4"))

## get the list of KOs of two KEGG genes.
getKOFromKGene(c("hsa:5232","hsa:5224"))

## get the list of KEGG genes of KOs. 
getKGeneFromKO(c("ko:K00927","ko:K01834"))

## get the list of pathways that two KEGG genes belong to.
getPathwayFromKGene(c("hsa:5232","hsa:5224"))

## get the list of KEGG genes in pathway 00010.
getKGeneFromPathway(c("path:00010"))


## get the list of KEGG genes of two official symbols.
getKGeneFromSymbol(c("PGK2","PGAM2"))

## get the list of symbols of two KEGG genes.
getSymbolFromKGene(c("hsa:5232","hsa:5224"))

## get the three- or four-letter KEGG organism code according to organism name.
getOrgIdFromOrgName("S.cerevisiae")

## get organism name according to the three- or four-letter KEGG organism code.
getOrgNameFromOrgId("sce")




## get the list of enzymes of two genes (ncbi-geneid).
getEnzymeFromGene(c("5232","5224"))

## get the list of gene (ncbi-geneid) of two enzymes. 
getGeneFromEnzyme(c("ec:5.4.2.1","ec:2.7.2.3"))

## get the list of KOs of two genes (ncbi-geneid).
getKOFromGene(c("5232","5224"))

## get the list of gene (ncbi-geneid) of two KOs.
getGeneFromKO(c("ko:K01834","ko:K00927"))

## get the list of pathways that two genes belong to.
getPathwayFromGene(c("5224","5232"))

## get the list of genes (ncbi-geneid) in pathway 00010.
getGeneFromPathway("path:00010")

## get the list of symbols of two genes (ncbi-geneid).
getSymbolFromGene(c("5224","5232"))

## get the list of genes (ncbi-geneid) of two official symbols.
getGeneFromSymbol(c("PGK2","PGAM2"))

## get the list of enzymes of two official symbols.
getEnzymeFromSymbol(c("PGK2","PGAM2"))

## get the list of official symbols of enzymes.
getSymbolFromEnzyme(c("ec:5.4.2.1","ec:2.7.2.3"))

## get the list of KOs of two official symbols.
getKOFromSymbol(c("PGK2","PGAM2"))

## get the list of official symbols of KOs.
getSymbolFromKO(c("ko:K01834","ko:K00927"))

} 
}
\keyword{ file }

