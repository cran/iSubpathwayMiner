\name{getNodeLabel}
\alias{getNodeLabel}
\title{Get new node labels according to names attributes of graph}
\description{
  Get new node labels according to names attributes of graph.
}
\usage{
  getNodeLabel(graph,type="symbol",displayNumber=1)
}
\arguments{
  \item{graph}{ The igraph object of a pathway graph.}
  \item{type}{ A character string. "symbol" or "currentId".}
  \item{displayNumber}{ An integer value.}

}
\details{
The function \code{getNodeLabel} is able to get new node labels according to names attributes of graph. The argument \code{type} defines type of new labels. Currently, we support two types of labels: symbol and currentId (see the function \code{\link{getOrgAndIdType}}). If type="symbol", gene node names of graph will be converted into KEGG gene identifiers and then find their symbols. Finally, symbols are displayed on the corresponding nodes. If nodes don't correspond to any symbol, then node labels will not be changed. The argument \code{displayNumber} decides the number of labels displayed on nodes in the graph. When there are many gene identifiers in a node, we can use the argument \code{displayNumber} to limit the number of gene identifiers displayed. The default value of the argument is 1. This means that if there are more than one gene in a node, the symbol "..." is added to the node.

}

\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\examples{
\dontrun{

#For a non-metabolic pathway graph, 
#change colors of edge with high betweenness as red 
path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/non-metabolic/ko/",sep="")
p1<-getPathway(path,c("ko04010.xml"))
g1<-getNonMetabolicGraph(p1)

#change node labels into symbol
plotGraph(g1[[1]],vertex.label=getNodeLabel)

#Get the type of organism and identifier in the current environment variable.
getOrgAndIdType()
#change node labels into the current gene identifiers
plotGraph(g1[[1]],vertex.label=getNodeLabel(g1[[1]],"currentId",1))

}
}
\keyword{ file }

