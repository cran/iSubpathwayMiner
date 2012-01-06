\name{simplifyGraph}
\alias{simplifyGraph}
\title{Simplify pathways to gene product graphs or compound graphs}
\description{
  Simplify pathways as graphs with only gene products (or only compound) as nodes
}
\usage{
  simplifyGraph(graphList,nodeType="geneProduct",directEdge=TRUE,verbose=FALSE)
}
\arguments{
  \item{graphList}{ A list of graphs.}
  \item{nodeType}{ A character string.}
  \item{directEdge}{ A logical. Should the direct edges remain?}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{
When we focus on gene products such as enzymes, genes, KOs, compounds may be not important. Similarly, gene products may be not important when we focus on metabolites (compounds). An useful approach is to get graphs with gene products (or compounds) as nodes and compound (gene products) as edges. We take an example of constructing metabolic pathway graphs with enzymes as nodes and compounds as edges. Firstly, all enzymes in a pathway graph are used as nodes. For undirected, two nodes are then connected by an edge if their corresponding reactions have a common compound. For directed, two nodes are connected by an edge if their corresponding reactions have a common compound and two nodes is reachable through the compound. Finally, compound information is added into edge attribute of new graphs. Similarly, a metabolic pathway graph can be converted to a graph with compounds as nodes. Two nodes are connected by an edge if they belong to the same reaction. Enzyme information is added into edge attribute of new graphs. 

The argument \code{directEdge} is used to determine whether direct edges should be obtained together.

The function can also be applied to non-metabolic pathway graphs. In non-metabolic pathways, most of nodes are gene products and only several nodes are compounds. Many edges belong to relation between gene products. However, some edges also belong to relation between gene products and compounds. For example, in the calcium signaling pathway, most of nodes are gene products, but several nodes are compounds (e.g., ca2+). These compounds usually interact with other gene products. Sometimes, for non-metabolic pathways, we may want to construct graphs with gene products as nodes and maintain these indirect edges that share the same compound. To do it, we can use the function \code{simplifyGraph} to simplify graphs to those with gene products as nodes and maintain two kinds of edges: direct edges and indirect edges. The values of argument \code{nodeType} and \code{directEdge} should be setted as "geneProduct" and TRUE. The function \code{simplifyGraph} can obtain indirected edges between gene products. When the value of argument is TRUE, direct edges between gene products can also be obtained. 

We can use the function \code{\link{filterNode}} to construct non-metabolic graphs with only gene products as nodes through filtering "map" and "compound". But only direct edges can be obtained. 
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}} 
\references{

Xia, J. and Wishart, D.S. (2010) Metpa: A Web-Based Metabolomics Tool for Pathway Analysis and Visualization, Bioinformatics, 26, 2342-2344.

Guimera, R. and Nunes Amaral, L.A. (2005) Functional Cartography of Complex Metabolic Networks, Nature, 433, 895-900.

Jeong, H., et al. (2000) The Large-Scale Organization of Metabolic Networks, Nature, 407, 651-654.

Koyuturk, M., Grama, A. and Szpankowski, W. (2004) An efficient algorithm for detecting frequent subgraphs in biological networks. Bioinformatics, 20 Suppl 1, i200-207.

Draghici, S., Khatri, P., Tarca, A.L., Amin, K., Done, A., Voichita, C., Georgescu, C. and Romero, R. (2007) A systems biology approach for pathway level analysis. Genome Res, 17, 1537-1545.

Schreiber, F. (2002) High quality visualization of biochemical pathways in BioPath. In Silico Biol, 2, 59-73.

Barabasi, A.L. and Oltvai, Z.N. (2004) Network Biology: Understanding the Cell's Functional Organization, Nature reviews, 5, 101-113.

}
\seealso{ \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{filterNode}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{mapNode}}, \code{\link{getUGraph}}, \code{\link{getSimpleGraph}} }
\examples{
\dontrun{
#####################metabolic##########################
##Convert metabolic pathways to graphs.
path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/metabolic/ec/",sep="")
g<-getMetabolicGraph(getPathway(path,c("ec00010.xml","ec00020.xml")))
##Construct graphs with compounds as nodes and enzymes as edges
g1<-simplifyGraph(g,nodeType="geneProduct")

#visualize
plotGraph(g1[[1]])

##Construct graphs with compounds as nodes and enzymes as edges
g2<-simplifyGraph(g,nodeType="compound")

#visualize
plotGraph(g2[[1]])

#####################non-metabolic##########################
#For a non-metabolic pathway graph, 
path1<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/non-metabolic/ko/",sep="")
gn<-getNonMetabolicGraph(getPathway(path1,c("ko04010.xml")))

#convert pathway to graphs with gene products as nodes
gn1<-simplifyGraph(gn,nodeType="geneProduct")

#visualize
plotGraph(gn1[[1]])

##misuse: only indirected edges are obtained
gn2<-simplifyGraph(gn,nodeType="geneProduct",directEdge=FALSE)
#visualize
plotGraph(gn1[[1]])

##use function filterNode: only directed edges are obtained
gn2<-filterNode(gn,nodeType=c("map","compound"))
#compare results: number of edges reduce because indirected edges are ignored 
ecount(gn1[[1]])
ecount(gn2[[1]])

#visualize
plotGraph(gn1[[1]])


}
}
\keyword{ file }

