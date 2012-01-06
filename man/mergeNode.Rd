\name{mergeNode}
\alias{mergeNode}
\title{Merge nodes with the same names}
\description{
  Merge nodes with the same names
}
\usage{
  mergeNode(graphList,simpleGraph=TRUE)
}
\arguments{
  \item{graphList}{A list of graph class.}
  \item{simpleGraph}{A logical. If true, the returned graphs are simple. See the function \code{\link{getSimpleGraph}}}
}
\details{
A pathway usually includes some nodes with the same names. For example, an enzyme may appear repeatedly in a pathway. The function \code{mergeNode} can merge those nodes with the same name, which represent the same compounds, enzymes, or genes. The edges of the merged nodes are all edges of initial nodes.

After nodes are merged, multiple edges or loops may appear. The argument \code{simpleGraph} can delete them by using the function \code{\link{getSimpleGraph}}. Therefore, the returned graphs will be simple.
}
 
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{ \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{filterNode}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{simplifyGraph}}, \code{\link{getUGraph}}, \code{\link{getSimpleGraph}} }
\examples{
\dontrun{
##Convert metabolic pathways to graphs.
path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/metabolic/ec/",sep="")
p<-getPathway(path,c("ec00010.xml","ec00020.xml"))
g<-getMetabolicGraph(p)
vcount(g[[1]])
g1<-mergeNode(g)
vcount(g1[[1]])
op<-par(mfrow = c(1, 2))
plotGraph(g[[1]],main="the initial graph")
plotGraph(g1[[1]],main="the merged graph")
par(op)

}
}
\keyword{ file }

