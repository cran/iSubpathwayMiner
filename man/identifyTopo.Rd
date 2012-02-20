\name{identifyTopo}
\alias{identifyTopo}
\title{Topology-based pathway analysis of gene sets}
\description{
  Topology-based pathway analysis of gene sets.
}
\usage{
identifyTopo(componentList,graphList,type="gene",propertyName="degree",order="pvalue",decreasing=FALSE,
  degree.mode="total",loops =TRUE,betweenness.directed=TRUE,clusteringCoefficient.type="local",closeness.mode="all",
  locateOrg=TRUE,ignoreAmbiguousEnzyme=TRUE,alternative="two.sided",background=getBackground(type))
}
\arguments{
  \item{componentList}{A character vector of components.}
  \item{graphList}{A list. Its elements is a graph with the \code{igraph} type.} 
  \item{type}{A character string. Should be one of "gene", "compound", or "gene_compound".} 
  \item{propertyName}{A character string. Should be one of "degree", "betweenness", "clusteringCoefficient", or "closeness".}  
  \item{order}{A character string. Should be one of "pvalue" or "fdr".}
  \item{decreasing}{A logical. Should the sort order be increasing or decreasing?}
  \item{degree.mode}{A character string. Should be one of "all", "out", "in", or "total".}
  \item{loops}{A logical. When computing degree, are the loop edges considered? See \code{\link[igraph]{degree}}. }
  \item{betweenness.directed}{A logical. whether directed paths should be considered while determining the shortest paths.}
  \item{clusteringCoefficient.type}{A character string. The type of the clusteringCoefficient to calculate. See \code{\link[igraph]{transitivity}}.}
  \item{closeness.mode}{A character string. Should be one of "out", "in", or "all".}  
  \item{locateOrg}{A logical. Should only the current organism pathways be considered?}  
  \item{ignoreAmbiguousEnzyme}{A logical. Should ambiguous enzymes such as 2.6.1.- be ignored?} 
  \item{alternative}{A character string specifying the alternative hypothesis for wilcoxon rank sum test. Should be one of "two.sided", "greater" or "less". See the function }  
  \item{background}{A character vector. The default value is a set of all current components (e.g., all human genes).}   
}
\details{
The function maps interesting component sets onto each pathway. For the annotated components in a pathway, their topological property values can be computed. These values can be compared with property values of all background components in the pathway. Finally, the statistical significance can be computed using wilcoxon rank sum test. The available topological properties contain degree, clustering coefficient, betweenness, and closeness. Degree of a node is the number of its adjacent edges (see the function \code{\link[igraph]{degree}}). Local clustering coefficient quantifies the probability that the neighbours of a node are connected (see the function \code{\link[igraph]{transitivity}}). Node betweenness can be calculated based on the number of shortest path passing through a given node (see the function \code{\link[igraph]{betweenness}}). Closeness measures how many steps is required to access every other nodes from a given node (see the function \code{\link[igraph]{closeness}}). 

Currently, the function can support input of three kinds of cellular component sets: gene sets, compound (metabolite) sets, and gene and compound sets at the same time. Therefore, the system can provide topology-based pathway analysis of gene sets, compound sets, and the combined sets of genes and compounds. 

The argument \code{degree.mode} should be one of "all", "out", "in", and "total". "out" for out-degree, "in" for in-degree or "total" for the sum of the two. For undirected graphs this argument is ignored (see \code{\link[igraph]{degree}}). If the loop edges are considered while computing degree, the argument \code{loops} should be set as TRUE.

Note that a component may be annotated to multiple nodes in a pathway. For these annotated nodes, topological property values can be computed and their average value is considered as the topological property value of the component.

The argument \code{betweenness.directed} determines whether directed paths should be considered while determining the shortest paths.

Note that clusteringCoefficient.type should be set as "local". (see \code{\link[igraph]{transitivity}}).

The argument \code{closeness.mode} defines the types of the paths used for measuring the distance in directed graphs. "in" measures the paths to a vertex, "out" measures paths from a vertex, "all" uses undirected paths. This argument is ignored for undirected graphs (see the function \code{\link{closeness}}).

Note that the argument \code{ignoreAmbiguousEnzyme} is only applied to pathway graphs with enzymes as nodes. If TRUE, genes that ambiguous enzymes (e.g., "-") correspond to are ignored in analyses. See the function \code{\link{getGeneFromEnzyme}}.
}
\value{
A list. Each element of the list is another list. It includes following elements: 'pathwayId', 'pathwayName', 'annComponentList', 'annComponentNumber', 'annBgComponentList', 'annBgNumber', 'ComponentNumber', 'bgNumber', 'propertyName', 'annComponentPropertyValueList', 'propertyValue', 'annBgComponentPropertyValueList', 'bgPropertyValue', 'pvalue', and 'fdr'. They correspond to pathway identifier, pathway name, the submitted components annotated to a pathway, numbers of submitted components annotated to a pathway, the background components annotated to a pathway, numbers of background components annotated to a pathway, numbers of submitted components, numbers of background components, topological property name (e.g., 'degree'), topological property values of submitted components annotated to a pathway, average topological property values of submitted components annotated to a pathway, topological property values of the background components annotated to a pathway, average topological property values of the background components annotated to a pathway, p-value of wilcoxon rank sum test for 'annComponentPropertyValueList' and 'annBgComponentPropertyValueList', and Benjamini-Hochberg fdr values.

The background components annotated to a pathway are equal to all components in the pathway. For example, if the submitted components are human genes, the background components annotated to a pathway are equal to all human genes in the pathway.

The numbers of background components are the number of all components. For example, if the submitted components are human genes, numbers of background components is equal to all human genes.
	
To visualize and save the results, the \code{list} can be converted to the \code{data.frame} by the function \code{\link{printTopo}}. 

Note that \code{componentList} submitted by users must be a 'character' vector.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{printTopo}}}

\examples{
\dontrun{
##############topology-based pathway analysis of gene sets#######################
##metabolic
#get metabolic pathway graphs with enzymes as nodes and compounds as edges. 
#The result graphs are equal to the result of the following commands:
#metabolicEC<-get("metabolicEC",envir=k2ri)
#graphList<-simplifyGraph(filterNode(getMetabolicGraph(metabolicEC),
#   nodeType=c("map")),nodeType="geneProduct")

##We can convert pathways to the directed graphs with enzymes as nodes and compounds as edges.
##We can test topology properties of intersting gene sets in the types of pathways.
##degree represents enzymes that are encoded by interesting genes number of enzymes.

#get pathway graphs with enzymes as nodes.
graphList<-getMetabolicECECGraph()
#get a set of genes
geneList<-getExample(geneNumber=1000,compoundNumber=0)
#topology-based pathway analysis
ann<-identifyTopo(geneList,graphList,type="gene")
result<-printTopo(ann)
#print
result

#output with detail information
result<-printTopo(ann,detail=T)
write.table(result,file="result.txt",row.names=FALSE,sep="\t")

#For above directed pathway graphs, we can calculate "out" degree
ann<-identifyTopo(geneList,graphList,type="gene",degree.mode="out")
printTopo(ann)

#For pathway graphs, we can calculate clustering coefficient
ann<-identifyTopo(geneList,graphList,type="gene",propertyName="clusteringCoefficient")
result<-printTopo(ann)
result[,c(1,3,4,5,6,7,8)]

#For pathway graphs, we can calculate betweenness
ann<-identifyTopo(geneList,graphList,type="gene",propertyName="betweenness")
result<-printTopo(ann)
result[,c(1,3,4,5,6,7,8)]

#For pathway graphs, we can calculate closeness
ann<-identifyTopo(geneList,graphList,type="gene",propertyName="closeness")
result<-printTopo(ann)
result[,c(1,3,4,5,6,7,8)]



##We can convert pathways to the undirected graphs with only genes as nodes.
##In a graph, each node contains only a gene and each gene only appear once.
##We can test topology properties of intersting gene sets in the types of pathways.
##degree of a gene represents number of adjacent genes.

#get undirected pathway graphs with genes as nodes.
graphList<-getMetabolicGEGEUEMGraph(type="EC")
#get a set of genes
geneList<-getExample(geneNumber=1000,compoundNumber=0)
#topology-based pathway analysis
ann<-identifyTopo(geneList,graphList,type="gene",propertyName="degree")
result<-printTopo(ann)
#print
result[,c(1,3,4,5,6,7,8)]







##non-metabolic
#get non-metabolic pathway graphs with KOs as nodes. 
graphList<-getNonMetabolicKOKOGraph()
#get a set of genes
geneList<-getExample(geneNumber=100,compoundNumber=0)
#topology-based pathway analysis
ann<-identifyTopo(geneList,graphList,type="gene")
#print
printTopo(ann)

##############topology-based pathway analysis of compound sets#######################

##metabolic compound
#get metabolic pathway graphs with compounds as nodes. 
graphList<-getMetabolicCOCOGraph()
#get a set of genes
compoundList<-getExample(geneNumber=0,compoundNumber=100)
#topology-based pathway analysis
ann<-identifyTopo(compoundList,graphList,type="compound",propertyName="betweenness")
#data.frame
result<-printTopo(ann)
#print
result[,c(1,3,4,5,6,7,8)]

##############topology-based pathway analysis of gene and compound sets#######################
#get metabolic pathway graphs with enzymes and compounds as nodes. 
graphList<-getMetabolicECCOGraph()
#get a set of genes
componentList<-getExample(geneNumber=1000,compoundNumber=100)
#topology-based pathway analysis
ann<-identifyTopo(componentList,graphList,type="gene_compound",propertyName="betweenness")
#data.frame
result<-printTopo(ann)
#print
result[1:10,c(1,3:8)]


}
}
\keyword{ file }

