\name{expandNode}
\alias{expandNode}
\title{Expand nodes according to the type of nodes}
\description{
  Expand nodes with multiple molecules in pathway graphs according to the type of nodes.
}
\usage{
  expandNode(graphList,nodeType=c("ortholog","enzyme","gene","compound","map"))
}
\arguments{
  \item{graphList}{ A list of graph class.}
  \item{nodeType}{ A character vector. Which types of nodes should be expanded?}

}
\details{
In many pathway maps, some nodes may have multiple molecules. For example, node Ubiquinone, which is the compound node in Ubiquinone and other terpenoid-quinone biosynthesis pathway (ec00130), maps to four compounds: Ubiquinone (cpd:C00399), Ubiquinone-6 (cpd:C17568), Ubiquinone-8 (cpd:C17569) and Ubiquinone-9 (cpd:C01967). Node PDE, which is the enzyme node in Purine metabolism (ec00230), maps to two enzymes: PDE (ec:3.1.4.17) and cGMP-PDE (ec:3.1.4.35). Node PDE, which is the human gene node in Purine metabolism (hsa00230), maps to genes: hsa:10846, hsa:27115, hsa:50940, hsa:5136, hsa:5137, hsa:5138, hsa:5139, hsa:5140, hsa:5141, hsa:5142, hsa:5143, hsa:5144, hsa:5145, hsa:5146, hsa:5147, hsa:5148, hsa:5149, hsa:5150, hsa:5151, hsa:5152 hsa:5153, hsa:5158, hsa:8622, and hsa:8654. Node ADH (1.1.1.1), which is the KO node in Glycolysis / Gluconeogenesis metabolic pathway (ko00010), maps to three KOs: adh (ko:K00001), adhC (ko:K00121), and adhE (ko:K04072). Node RASGRF, which is the KO node in MAPK signaling pathway (ko04010), maps to two KOs: RASGRF1 (ko:K04349) and RASGRF1 (ko:K12326). Node RASGRF, which is the human gene node in MAPK signaling pathway (hsa04010), maps to two genes: RASGRF1 (hsa:5923) and RASGRF1 (hsa:5924).

Since one node in KEGG pathway does not necessarily map to only one molecule, to expand them should be useful for some special analyses. The function expandNode is just used to expand those nodes with multiple molecules.

Users can select which types of nodes are expanded. The argument \code{nodeType} is just used to identify the expanded nodes according to the type of nodes. The default values represent that all nodes are expanded. 
}
 
\references{
Draghici, S., Khatri, P., Tarca, A.L., Amin, K., Done, A., Voichita, C., Georgescu, C. and Romero, R. (2007) A systems biology approach for pathway level analysis. Genome Res, 17, 1537-1545.

Zhang, J.D. and Wiemann, S. (2009) KEGGgraph: a graph approach to KEGG PATHWAY in R and bioconductor. Bioinformatics, 25, 1470-1471.

KEGG Markup Language manual \url{http://www.genome.jp/kegg/xml/docs/}
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{ \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{filterNode}}, \code{\link{mergeNode}}, \code{\link{mapNode}}, \code{\link{simplifyGraph}}, \code{\link{getUGraph}}, \code{\link{getSimpleGraph}} }
\examples{
\dontrun{
##expand non-metabolic graphs
path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/non-metabolic/ko/",sep="")
p<-getPathway(path,c("ko04010.xml"))
g<-getNonMetabolicGraph(p,ambiguousEdgeDirection="bi-directed")
#node number before nodes are expanded
vcount(g[[1]])
##expand nodes in Graphs
g1<-expandNode(g)
#node number after nodes are expanded
vcount(g1[[1]])

##display the expanded graph
print(v=TRUE,g=TRUE,e=TRUE,g1[["04010"]])


##expand metabolic graphs
path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/metabolic/ko/",sep="")
p1<-getPathway(path,c("ko00010.xml"))
g2<-getMetabolicGraph(p1)
vcount(g2[[1]])
ecount(g2[[1]])
##expand nodes in Graphs
g3<-expandNode(g2)
#node number
vcount(g3[[1]])
ecount(g3[[1]])
#only expand nodes with type= "enzyme" or "ortholog" and "compound" in Graphs
g4<-expandNode(g2,nodeType=c("ortholog","enzyme","compound"))
#node number
vcount(g4[[1]])
#only expand nodes with type="enzyme" or "ortholog" in graphs
g5<-expandNode(g2,nodeType=c("ortholog","enzyme"))
#node number
vcount(g5[[1]])

}
}
\keyword{ file }

