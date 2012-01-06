\name{updatePathway}
\alias{updatePathway}
\alias{importPathway}
\title{Update or import pathway files}
\description{
  Update pathways from KEGG ftp site or import pathways from the local file systems.
}
\usage{
updatePathway(path="ftp://ftp.genome.jp/pub/kegg/xml/",verbose=TRUE)
importPathway(path,verbose=TRUE)
}
\arguments{
  \item{path}{A character string for setting path of Pathway files.}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{
The function \code{updatePathway} is used as update of the pathway variable \code{metabolicEC}, \code{metabolicKO}, and \code{nonMetabolicKO} in the environment variable \code{\link{k2ri}} through downloading data from KEGG ftp site. After the function is executed, information in pathways such as nodes, relation, reaction, etc. may be changed due to new change for these pathways from KEGG. Note that several pathways may disappear because they may be deleted or merged from new version of KEGG. The detailed information of pathway change is provided in \url{http://www.genome.jp/kegg/docs/upd_map.html}. Currently, we can't use R to know all new added pathway information in ftp site. Therefore, the function can't get new added pathways.
 
The function \code{importPathway} is used as construction of the pathway variable code{metabolicEC}, \code{metabolicKO}, and \code{nonMetabolicKO} in the environment variable \code{\link{k2ri}} from local system. The function can get up-to-date pathways including merged, changed, and even new pathways. Firstly, users need to download KGML pathway files from KEGG ftp site. They can be obtained from \url{ftp://ftp.genome.jp/pub/kegg/xml/kgml/metabolic/ec}, \url{ftp://ftp.genome.jp/pub/kegg/xml/kgml/metabolic/ko}, and \url{ftp://ftp.genome.jp/pub/kegg/xml/kgml/non-metabolic/ko}. Then, store them in fixed-format folder. We provide the example of the folder. The following path can see it: 

\preformatted{
path<-paste(system.file(package="iSubpathwayMiner"),"/localdata/",sep="")
}

Further, The function \code{\link{saveK2ri}} and \code{\link{loadK2ri}} can help save and load the new environment variable k2ri.



}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{updateOrgAndIdType}}, \code{\link{getPathway}}, \code{\link{saveK2ri}}, \code{\link{loadK2ri}}}

\examples{
\dontrun{
##update pathways
updatePathway()

##import pathways
path<-paste(system.file(package="iSubpathwayMiner"),"/localdata/",sep="")
importPathway(path)
}
}
\keyword{ file }

