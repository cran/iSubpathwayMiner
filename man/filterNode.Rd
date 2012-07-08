\name{filterNode}
\alias{filterNode}
\title{filter node according to node type}
\description{
  filter node according to node type
}
\usage{
  filterNode(graphList,nodeType=c("map"))
}
\arguments{

  \item{graphList}{ A list of graph class.}
  \item{nodeType}{ A character vector.}

}
\details{
The function is used to filter "not interesting" nodes. For example, it may be necessary to ignore nodes with type="map" when focusing on molecules such as compounds and gene products. Moreover, it may be necessary to ignore nodes with type="map" and "compound" when only focusing on gene products. 

The function will delete nodes according to the argument \code{nodeType} and thus related edges are also deleted. The function firstly find nodes that meet the argument \code{nodeType}. Secondly, these nodes are deleted by using the function \code{\link{delete.vertices}}. Finally, we will get subgraphs of all graphs that meet conditions. 

Note that the function can't get gene product graph with gene products as nodes and compound as edges for metabolic pathways. The function \code{\link{simplifyGraph}} can do it.

The argument \code{nodeType} may be one or combination of "ortholog", "enzyme", "gene", "compound", and "map".
}
 
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{ \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{mapNode}}, \code{\link{simplifyGraph}}, \code{\link{getUGraph}}, \code{\link{getSimpleGraph}} }
\examples{
\dontrun{

##Convert non-metabolic pathways to graphs, in which ambiguous edges
# is considered as "<-->"
path<-paste(system.file(package="iSubpathwayMiner"),"/localdata/kgml/non-metabolic/ko/",sep="")
p<-getPathway(path,c("ko04010.xml"))
g<-getNonMetabolicGraph(p,ambiguousEdgeDirection="bi-directed")
#type of each node before nodes are filtered
V(g[[1]])$type
#delete nodes with type="map"
g1<-filterNode(g,nodeType=c("map"))
V(g1[[1]])$type
#delete nodes whose types are not genes or gene products
g2<-filterNode(g,nodeType=c("compound","map"))
V(g2[[1]])$type

##Convert metabolic pathways to graphs.
path<-paste(system.file(package="iSubpathwayMiner"),"/localdata/kgml/metabolic/ec/",sep="")
p1<-getPathway(path,c("ec00010.xml","ec00020.xml"))
g3<-getMetabolicGraph(p1)
#type of each node before nodes are filtered
V(g3[[1]])$type
#delete nodes with type="map"
g4<-filterNode(g3,nodeType=c("map"))
V(g4[[1]])$type
}
}
\keyword{ file }

