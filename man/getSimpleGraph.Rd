\name{getSimpleGraph}
\alias{getSimpleGraph}
\title{Convert multiple graphs to simple graphs}
\description{
Convert multiple graphs to simple graphs
}
\usage{
  getSimpleGraph(graphList)
}
\arguments{

  \item{graphList}{ A list of graph class}
}
\details{
A loop edge is an edge where the two endpoints have the same node (vertex). Two edges are multiple edges if they have exactly the same two endpoints. If a graph is simple, it does not contain loop or/and multiple edges. The function \code{\link{is.simple}} can check whether a graph is simple. If not simple, we can use the function getSimpleGraph to get a simple graph. 

Note that the function can keep information of mutiple edges in edge attribute using ";" as separator.

}
 
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{ \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{filterNode}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{mapNode}}, \code{\link{simplifyGraph}}, \code{\link{getUGraph}}}
\examples{
\dontrun{

##Convert maps to graphs, in which ambiguous edges is considered as "<-->"
path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/non-metabolic/ko/",sep="")
p<-getPathway(path,c("ko04330.xml","ko04140.xml"))
g<-getNonMetabolicGraph(p,ambiguousEdgeDirection="bi-directed")

#check whether all graphs in a graph list are simple
all(sapply(g,is.simple))
#get simple graphs
g1<-getSimpleGraph(g)


g2<-getNonMetabolicGraph(get("nonMetabolicKO",envir=k2ri),ambiguousEdgeDirection="bi-directed")
#check whether all non-metabolic pathway graphs are simple
all(sapply(g2,is.simple))
#display non-metabolic pathway graph that belong to mutiple graphs
sapply(g2[!sapply(g2,is.simple)],function(x) x$title)

#check whether all non-metabolic pathway graphs are simple
g3<-getMetabolicGraph(get("metabolicKO",envir=k2ri))
#check whether all graphs in a graph list are simple
all(sapply(g3,is.simple))



}
}
\keyword{ file }

