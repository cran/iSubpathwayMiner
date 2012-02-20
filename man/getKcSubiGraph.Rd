\name{getKcSubiGraph}
\alias{getKcSubiGraph}
\title{Mine subpathways}
\description{
  Mine subpathways by using the k-clique concept in social network analysis.
}
\usage{
     getKcSubiGraph(k=4,graphList)
}
\arguments{
  \item{k}{An integer. A distance similarity parameter.}
  \item{graphList}{An undirected graph list.}
}
\details{
The function \code{getKcSubiGraph} uses the k-clique concept in social network analysis to mine subpathways. In social network analysis, a k-clique in a graph is a subgraph where the distance between any two nodes is no greater than k. 

The k-cliques-based method to identify subpathways firstly divide entire pathways to subpathways according to the parameter k. Secondly, the identifier of each subpathway is given with its pathway identifier plus a subpathway number. For each subpathway (k-clique), distance among nodes within it is no greater than the parameter k (a user-defined distance). Component sets can then be annotated to these subpathways. Finally, the significantly enriched subpathways can be identified using hypergeometric test.

Note that the k-clique subpathway identification method is flexible for choice of pathway data. When we choose the pathway data obtained from the function \code{\link{getMetabolicECECUEMGraph}}, the subpathway identification method corresponds to that in the SubpathwayMiner system.

}
\value{
  A list of graphs.
}
\references{
Wasserman,S. and Faust,K. (1994) Social network analysis: methods and applications. Cambridge University Press., New York, America.

Li, C., Li, X., Miao, Y., Wang, Q., Jiang, W., Xu, C., Li, J., Han, J., Zhang, F., Gong, B. et al. (2009) SubpathwayMiner: a software package for flexible identification of pathways. Nucleic Acids Res, 37, e131.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{identifyGraph}}, \code{\link{getMetabolicECECUGraph}}}

\examples{
\dontrun{
#############identify subpathways based on gene sets###########
############enzyme-enzyme###########
##identify metabolic subpathways based on gene sets
#get undirected metabolic pathway graphs with enzymes as nodes and compounds as edges
graphList<-getMetabolicECECUGraph()
#get all 4-clique subgraphs
subGraphList<-getKcSubiGraph(k=4,graphList)
#get a set of genes
geneList<-getExample(geneNumber=1000,compoundNumber=0)
#annotate gene sets to graphs 
#and identify significant graphs
ann<-identifyGraph(geneList,subGraphList,type="gene")
#convert ann to data.frame
result<-printGraph(ann)

#print the results to screen
result
#output of results
write.table(result,file="result.txt",row.names=FALSE,sep="\t")

###metaboli subpathway identification method in the SubpathwayMiner system###
#get undirected metabolic pathway graphs with enzymes as nodes and compounds as edges
#each node contains only one enzyme and each enzyme only appear once.
graphList<-getMetabolicECECUEMGraph()
#get all 4-clique subgraphs
subGraphList<-getKcSubiGraph(k=4,graphList)
##get a set of compounds and genes
componentList<-getExample(geneNumber=1000,compoundNumber=0)
#annotate gene sets and identify significant subpathways
ann<-identifyGraph(componentList,subGraphList,type="gene")
result<-printGraph(ann)

###############KO-KO###################
##identify metabolic subpathways based on gene sets
#get the KO-KO graphs
graphList<-getMetabolicKOKOUGraph()
#get all 4-clique subgraphs
subGraphList<-getKcSubiGraph(k=4,graphList)
#get a set of compounds
geneList<-getExample(geneNumber=1000,compoundNumber=0)
#annotate gene sets to graphs 
#and identify significant graphs
ann<-identifyGraph(geneList,subGraphList,type="gene")
#convert ann to data.frame
result<-printGraph(ann)
result[1:10,c(1,3,4,5)]

##identify non-metabolic subpathways based on gene sets
#get the KO-KO graphs
graphList<-getNonMetabolicKOKOUGraph()
#get all 4-clique subgraphs
subGraphList<-getKcSubiGraph(k=4,graphList)
#get a set of compounds
geneList<-getExample(geneNumber=1000,compoundNumber=0)
#annotate gene sets to graphs 
#and identify significant graphs
ann<-identifyGraph(geneList,subGraphList,type="gene")
#convert ann to data.frame
result<-printGraph(ann)
result[1:10,c(1,3,4,5)]

##identify metabolic and non-metabolic subpathways based on gene sets
#get the KO-KO graphs
graphList1<-getMetabolicKOKOUGraph()
graphList2<-getNonMetabolicKOKOUGraph()
graphList<-c(graphList1,graphList2)
#get all 4-clique subgraphs
subGraphList<-getKcSubiGraph(k=4,graphList)
#get a set of compounds
geneList<-getExample(geneNumber=1000,compoundNumber=0)
#annotate gene sets to graphs 
#and identify significant graphs
ann<-identifyGraph(geneList,subGraphList,type="gene")
#convert ann to data.frame
result<-printGraph(ann)
result[1:10,c(1,3,4,5)]


#############identify metabolic subpathways based on compound sets###########
##identify metabolic subpathways based on compound sets
#get the compound-compound graphs
graphList<-getMetabolicCOCOUGraph()
#get all 4-clique subgraphs
subGraphList<-getKcSubiGraph(k=4,graphList)
#get a set of compounds
compoundList<-getExample(geneNumber=0,compoundNumber=100)
#annotate compound sets and identify significant graphs
ann<-identifyGraph(compoundList,subGraphList,type="compound")
#convert ann to data.frame
result<-printGraph(ann)
result[1:10,c(1,3,4,5)]



}
}
\keyword{ file }

