\name{updateOrgAndIdType}
\alias{updateOrgAndIdType}
\title{Update the organism and the type of gene identifier}
\description{
  Update the organism and the type of gene identifier.
}
\usage{
updateOrgAndIdType(org=getOrgAndIdType()[1],idType=getOrgAndIdType()[2],
  path="ftp://ftp.genome.jp/pub/kegg/genes/organisms",symbolData=TRUE,verbose=TRUE)
}
\arguments{
  \item{org}{A character string. The abbreviation of a genome name.}
  \item{idType}{A character string. The type of gene identifier.}
  \item{path}{A character string.}
  \item{symbolData}{A logical. Should the gene symbol data be downloaded?}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{

The function is able to update the variable \code{gene2ec}, \code{gene2ko}, \code{keggGene2gene}, etc. 

Note that if the user don't run the function \code{\link{updateOrgAndIdType}} or \code{\link{loadK2ri}} after starting up R-system and loading the package of the system, then the defalut value of the argument \code{org} and \code{idType} is "hsa" (human) and "ncbi-geneid" (Entrez gene identifiers). The user can get the information from the return value of the function \code{\link{getOrgAndIdType}}.

The argument \code{org} must be the abbreviation of a genome name. For example, the hsa, mms, eco, sce and dme is the abbreviation of human, mouse, E.coli, yeast and fluit fly. Detailed information is provided in \url{http://www.genome.jp/kegg/catalog/org_list.html}. The function \code{\link{getOrgIdFromOrgName}} and \code{\link{getOrgNameFromOrgId}} can also provide conversion between organism name and the three- or four-letter KEGG organism code.

The argument \code{idType} is a character string of the type of identifier. The system supports most KEGG cross-reference identifiers such as Entrez gene IDs (idType="ncbi-geneid"), NCBI gi numbers (idType="ncbi-gi"), UniProt accession numbers (idType="uniprot"), etc. Detailed information is provided in \url{ftp://ftp.genome.jp/pub/kegg/genes/organisms}. For example, because a file name in "hsa" file directory is "hsa ensembl-hsa.list", idType="ensembl-hsa" is available as the input identifier type. Note that the \code{idType} is relative to the genome. Different genomes may support different \code{idType}. For example, "sgd-sce" is supported by yeast. However, it is not supported by human.

The argument \code{path} is the path of file directory of the organism cross-reference identifiers. The default value is "ftp://ftp.genome.jp/pub/kegg/genes/organisms". The setting ensures that the user is able to obtain the up-to-date data from the KEGG FTP site. Of course, the user can also download the organism data of interest from FTP site and change path to the data file for implementing the local update.

Note that the programming may be time consuming. Although the running time of the function is usually less than 2 mins, it depends on your internet connection speed because files need to be downloaded from KEGG ftp site (\url{ftp://ftp.genome.jp/pub/kegg/genes/organisms}).
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{getOrgAndIdType}},\code{\link{updatePathway}}, \code{\link{importPathway}}, \code{\link{saveK2ri}}, \code{\link{loadK2ri}}}

\examples{
\dontrun{
##update organism and the type of gene identifiers
getOrgAndIdType()
updateOrgAndIdType("sce","sgd-sce")
getOrgAndIdType()
updateOrgAndIdType("hsa","ncbi-geneid")
getOrgAndIdType()
}
}
\keyword{ file }

