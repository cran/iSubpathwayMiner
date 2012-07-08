\name{plotAnnGraph}
\alias{plotAnnGraph}
\title{Visualize annotated pathway graphs}
\description{
  Visualize annotated pathway graphs.
}
\usage{
    plotAnnGraph(pathwayId,graphList,ann,gotoKEGG=FALSE,
	  orgSpecific=TRUE,multipleCell=FALSE,displayInR=TRUE,match=TRUE,vertex.frame.color="red",...)
}
\arguments{
  \item{pathwayId}{An character vector of pathway identifier, e.g., c("path:00010","path:00020").}
  \item{graphList}{An graph list.}
  \item{ann}{An list, e.g., the return value of the function \code{\link{identifyGraph}}.}
  \item{gotoKEGG}{A logical. If TRUE, go to KEGG web site. The identified nodes are colored red.}
  \item{orgSpecific}{A logical. If TRUE and gotoKEGG=TRUE, show org specific pathways.} 
  \item{multipleCell}{A logical. If TRUE and displayInR=TRUE, display all pathway graphs on the device.}
  \item{displayInR}{A logical. If TRUE, display pathway graphs on the device.}
  \item{match}{A logical. The default value is TRUE. This means that pathwayId must be equal to Id in graphList. If FALSE, only entire pathway part of pathwayId match Id in graphList.}
  \item{vertex.frame.color}{A character string of color of identified nodes. The default value is "red".}
  \item{...}{ The arguments passed to or from methods. See \code{\link{plotGraph}}, \code{\link{plot.igraph}} and see \code{\link{plot}}.}  
}
\details{
The function can visualize the annotated pathway (or subpathways) graphs. The red nodes in the result graph represent the nodes which include the submitted molecules.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{plotGraph}},\code{\link{identifyGraph}},\code{\link{printGraph}}}

\examples{
\dontrun{
##visualize annotated graphs
##Convert all metabolic pathways to graphs.
metabolicKO<-get("metabolicKO",envir=k2ri)
graphList<-getMetabolicGraph(metabolicKO)

##get a set of genes
geneList<-getExample(geneNumber=1000)

#annotate gene sets to graphs 
#and identify significant graphs
ann<-identifyGraph(geneList,graphList)
#visualize annotated graph
plotAnnGraph("path:00010",graphList,ann,gotoKEGG=FALSE)
#orgSpecific=TRUE
plotAnnGraph("path:00010",graphList,ann,gotoKEGG=TRUE)
#orgSpecific=FALSE
plotAnnGraph("path:00010",graphList,ann,gotoKEGG=TRUE,orgSpecific=FALSE)



nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
graphList<-getNonMetabolicGraph(nonMetabolicKO)

##get a set of genes
geneList<-getExample(geneNumber=1000)

#annotate gene sets to graphs 
#and identify significant graphs
ann<-identifyGraph(geneList,graphList)
plotAnnGraph("path:04010",graphList,ann,vertex.label=getNodeLabel)

plotAnnGraph("path:04010",graphList,ann,vertex.label=getNodeLabel,gotoKEGG=TRUE)

}
}
\keyword{ file }

