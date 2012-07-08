\name{getLocSubGraph}
\alias{getLocSubGraph}
\title{Locate subpathways bu using lenient distance similarity of signature nodes}
\description{
  Locate subpathways by using lenient distance similarity of signature nodes.
}
\usage{
  getLocSubGraph(moleculeList,graphList,type="gene_compound",
                       n=5,s=5,method="shortestPaths",ignoreAmbiguousEnzyme=TRUE)
}
\arguments{
  \item{moleculeList}{A character vector of molecules.}
  \item{graphList}{A list. Its elements is a \code{igraph}.} 
  \item{type}{A character string. Should be one of "gene", "compound", and "gene_compound".}
  \item{n}{An integer. The maximum number of intermediate nodes.}  
  \item{s}{An integer. The subgraphs in which number of nodes are more than s are considered as a subpathway.}  
  \item{method}{A character string. which shorest path algorithms will be used. See the function \code{\link{get.shortest.paths}}.} 
  \item{ignoreAmbiguousEnzyme}{A logical. Should ambiguous enzymes such as 2.6.1.- be ignored?}    
}
\details{
We use lenient distance similarity of signature nodes (interesting genes and metabolites in a pathway) to locate subpathways. Specially, for each pathway that contains signature nodes, we compute the shortest path between any two signature nodes in the given pathway graph. If the shortest path length between two signature nodes is shorter than n+1, then the two signature nodes and other nodes at the shortest path are added into the same node set. We extract the corresponding subgraph in the pathway graph according to each node set. We finally define these subgraphs with node number >= s as the subpathway regions of the pathway. 

The argument \code{n} is maximum permitted intermediate (that is, non-signature) node number at the shortest path between signature nodes.

The argument \code{s} filter subpathways in which the number of nodes are less than the parameter s. The default value is five, which means that subpathways with the number of nodes <5 are deleted.

The argument \code{method} determins which shorest path algorithms will be used. The value is "get.shortest.paths" or "get.all.shortest.paths". If "get.shortest.paths", only a path will be obtained. If "get.shortest.paths", multiple paths will be obtained when there exist multiple shortest paths. See the function \code{\link{get.shortest.paths}}.

Note that the argument \code{ignoreAmbiguousEnzyme} is only applied to pathway graphs with enzymes as nodes. If TRUE, genes that ambiguous enzymes (e.g., "-") correspond to are ignored in analyses. See the function \code{\link{getGeneFromEnzyme}}.
}
\value{
  A list of graphs.
}
\author{Chunquan Li <lcqbio@yahoo.com.cn>}
\seealso{\code{\link{SubpathwayGM}},\code{\link{identifyGraph}},\code{\link{plotAnnGraph}},\code{\link{printGraph}},\code{\link{getMetabolicKOCOGraph}}}

\examples{
\dontrun{

#convert all metabolic pathways to graphs with KOs and compounds as nodes.
graphList<-getMetabolicKOCOGraph()
#get interesting gene and metabolite sets
moleculeList<-getExample(geneNumber=1000,compoundNumber=100)
#locate subpathways
subGraphList<-getLocSubGraph(moleculeList,graphList,type="gene_compound",n=5,s=5)
#annotate gene and metabolite sets to subpathways
#and identify significant subpathways
ann<-identifyGraph(moleculeList,subGraphList,type="gene_compound")
result<-printGraph(ann,detail=TRUE)
#display part of results
result[1:10,c(1:6)]
write.table(result,file="result.txt",row.names=FALSE,sep="\t")

#visualize
plotAnnGraph("path:00010_1",subGraphList,ann,displayInR=TRUE,gotoKEGG=TRUE)

##Case: use the function SubpathwayGM to identify metabolic subpathways associated with colorectal cancer
##use new background
#convert all metabolic pathways to graphs with KOs and compounds as nodes.
graphList<-getMetabolicKOCOGraph()
#read differential genes and metabolites in colorectal cancer from files
path1<-paste(system.file(package="iSubpathwayMiner"),"/localdata/crc_diff_gene.txt",sep="")
geneList<-as.character(read.table(path1,sep="\t")[[1]])
path2<-paste(system.file(package="iSubpathwayMiner"),"/localdata/crc_diff_metabolite.txt",sep="")
metaboliteList<-as.character(read.table(path2,sep="\t")[[1]])
moleculeList<-c(geneList,metaboliteList)
#locate subpathways
subGraphList<-getLocSubGraph(moleculeList,graphList,type="gene_compound",n=5,s=5)
#annotate gene and metabolite sets to subpathways
#and identify significant subpathways
ann<-identifyGraph(moleculeList,subGraphList,type="gene_compound",background=getPrioBackground(type="gene_compound"))
result<-printGraph(ann)
#display part of results
result[which(result[,"pvalue"]<0.01),]

result1<-printGraph(ann,detail=TRUE)
#write the identification result of subpathways to tab delimited file. 
write.table(result1,file="result1.txt",row.names=FALSE,sep="\t")

}
}
\keyword{ file }

