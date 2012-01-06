\name{mapNode}
\alias{mapNode}
\title{Map current organism-specific gene identifiers to nodes in pathway graphs}
\description{
  Map current organism-specific gene identifiers to nodes in pathway graphs
}
\usage{
  mapNode(graphList)
}
\arguments{

  \item{graphList}{ A list of graph class.}

}
\details{
The function can map current organism-specific gene identifiers to nodes of graphs. In particular, the names attribute-value of nodes that can correspond to gene IDs are replaced using gene IDs. Accordingly, the values of type attribute of nodes are revised to "gene". The values of attribute graphics_name, graphics_fgcolor, and graphics_bgcolor are also resived. Finally, pathway graphs will be converted to organism-specific graphs as well as idType-specific ones. 

The method is similar to the method used in KEGG to get organism-specific pathways according to reference pathways. Our method can map different identifiers (e.g., NCBI-geneid, NCBI-gi) to nodes. However, in xml files of KEGG organism-specific pathways, only KEGG gene identifiers are mapped to nodes. 

Since the function can convert reference pathway graphs to organism-specific pathway graphs, the method may be more convenient in obtaining organism-specific pathways than those provided in KEGG. Moreover, NCBI-geneid is also more popular compared with KEGG gene IDs. More important, our system can support most organism-specific pathway graphs through using the method to map gene IDs to reference pathway graphs and the mehod to update organisms \code{\link{updateOrgAndIdType}}.


}
 
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}

\references{
Hung, J.H., et al. (2010) Identification of Functional Modules That Correlate with Phenotypic Difference: The Influence of Network Topology, Genome Biol, 11, R23.
}
\seealso{ \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{filterNode}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{simplifyGraph}}, \code{\link{getUGraph}}, \code{\link{getSimpleGraph}} }
\examples{
\dontrun{

##Convert metabolic pathways to graphs.
path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/metabolic/ec/",sep="")
p<-getPathway(path,c("ec00010.xml","ec00020.xml"))
g<-getMetabolicGraph(p)

g1<-mapNode(g)



##Convert metabolic pathways to graphs.
path2<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/metabolic/ko/",sep="")
p2<-getPathway(path2,c("ko00010.xml","ko00020.xml"))
g4<-getMetabolicGraph(p2)

g5<-mapNode(g4)
}
}
\keyword{ file }

