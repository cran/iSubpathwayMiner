\name{getNonMetabolicGraph}
\alias{getNonMetabolicGraph}
\title{Convert non-metabolic pathways to directed graphs}
\description{
  Convert non-metabolic pathways to directed graphs. Edges are obtained from relations.
}
\usage{
getNonMetabolicGraph(pathwayList,ambiguousEdgeDirection="bi-directed",
  ambiguousEdgeList=c("unknow","compound","hidden compound","state change","binding/association","dissociation"),
  simpleGraph=TRUE,splitCompoundSubType=TRUE,verbose=FALSE)
}
\arguments{
  \item{pathwayList}{A list of pathways. The argument \code{pathwayList} is the returned value of the function \code{\link{getPathway}}.}
  \item{ambiguousEdgeDirection}{Should be one of "single", "bi-directed", and "delete". }
  \item{ambiguousEdgeList}{What types of edges are difined as ambiguousn edges?}
  \item{simpleGraph}{A logical. If TRUE, Convert multiple graphs to simple graphs when the result graphs are multiple graphs.}
  \item{splitCompoundSubType}{A logical. If TRUE, split edge with type=compound as two edges 
                              protein-compound and compound-protein.}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{
The function converts non-metabolic pathways to directed graphs. The returned graphs mainly contain two types of nodes: gene products (enzymes or KOs) and linked pathway maps. Sometimes, there are several compounds in pathways such as ca+. Edges are obtained from relations. Detailed information about nodes, relations is provided in \url{http://www.genome.jp/kegg/xml/docs/}. In particular, two nodes are connected by an edge if they have relationships get from relation element of the KGML file. The relation element specifies relationships between nodes. For example, the attribute PPrel represents protein-protein interaction such as binding and modification. Other information such as node attribute, pathway attribute, etc. is converted to attribute of graphs in igraph.

For non-metabolic pathways, there are many types of edges, which include ECrel, PPrel, GErel. Each types further has many subtypes such as compound, hidden compound, activation, inhibition, expression, repression, indirect effect, state change, binding/association, dissociation, and missing interaction. Detailed information is provided in \url{http://www.genome.jp/kegg/xml/docs/}. However, not all types of edges have definite direction. The argument \code{ambiguousEdgeDirection} provide flexible methods for defining direction according to subtype of nodes. The default ambiguous edge directions include "compound", "hidden compound", "state change", "binding/association", "dissociation". The value of the argument \code{ambiguousEdgeDirection} should be one of "single", "bi-directed", and "delete". They mean to convert ambiguous edges to "-->", "<-->", or to delete them. 
 
The argument \code{pathwayList} is the returned value of the function \code{\link{getPathway}}.
 
Other detailed information is provided in \code{\link{getMetabolicGraph}}.
}
\value{
  A list of directed graphs
}

\references{
Draghici, S., Khatri, P., Tarca, A.L., Amin, K., Done, A., Voichita, C., Georgescu, C. and Romero, R. (2007) A systems biology approach for pathway level analysis. Genome Res, 17, 1537-1545.

Smart, A.G., Amaral, L.A. and Ottino, J.M. (2008) Cascading failure and robustness in metabolic networks. Proc Natl Acad Sci U S A, 105, 13223-13228.

KEGG Markup Language manual \url{http://www.genome.jp/kegg/xml/docs/}
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}} 
\seealso{ \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{filterNode}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{mapNode}}, \code{\link{simplifyGraph}}, \code{\link{getUGraph}}, \code{\link{getSimpleGraph}} }
\examples{
\dontrun{
###################################################
##Convert maps to graphs, in which ambiguous edges is considered as "<-->"
path<-paste(system.file(package="iSubpathwayMiner"),"/localdata/kgml/non-metabolic/ko/",sep="")
p<-getPathway(path,c("ko04010.xml"))
g<-getNonMetabolicGraph(p,ambiguousEdgeDirection="bi-directed")
}
}
\keyword{ file }

