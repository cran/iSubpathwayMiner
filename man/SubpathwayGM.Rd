\name{SubpathwayGM}
\alias{SubpathwayGM}
\title{Identify metabolic subpathways via joint power of genes and metabolites and their topology within pathways}
\description{
  SubpathwayGM is an attempt to more accurate level of pathway analysis through integrating information of genes and metabolites, their position and number of their cascade regions within the given pathway.
}
\usage{
  SubpathwayGM(moleculeList,n=5,s=5,background=getPrioBackground(type="gene_compound"),graphType="KO",...)
}
\arguments{
  \item{moleculeList}{A character vector of interesting genes and metabolites.}
  \item{n}{An integer. The maximum number of intermediate nodes. The default value is 5.}  
  \item{s}{An integer. The subgraphs in which number of nodes are more than s are considered as a subpathway. The default value is 5.}  
  \item{background}{A character vector of molecules used to identify the statistically significant subpathways.}
  \item{graphType}{A character string. Should be one of "KO" or "EC". The default value is "KO".}
  \item{...}{ The arguments passed to or from methods. See \code{\link{identifyGraph}} and see \code{\link{getLocSubGraph}}.}   
}
\details{

When users input interesting molecules (genes and metabolites), the function can identify metabolic subpathways mainly through: (1) Mapping metabolites and genes of interest to graphs of pathways after graph-based reconstruction of metabolic pathways; (2) Locating subpathways within pathways according to signature nodes. The step depends on lenient distance similarity of signature nodes within pathway structure to locate each potential metabolic subpathway regions through considering intermediate nodes in which each has close distance with at least one signature node in the subpathway. (3) Evaluating the statistical significance of subpathways. The hypergeometric test is used.

The argument \code{n} is maximum permitted intermediate (that is, non-signature) node number at the shortest path between signature nodes.

The argument \code{s} filter subpathways in which the number of nodes are less than the parameter s. The default value is five, which means that subpathways with the number of nodes <5 are deleted.

The argument \code{graphType} is the used pathway graph in the function \code{SubpathwayGM}. See the function \code{\link{getMetabolicKOCOGraph}} and \code{\link{getMetabolicGraph}}.

If users don't set the values of the argument \code{background}, the background distribution will be obtained from the whole genome and metabolome. Detailed information is provided in the function \code{\link{getPrioBackground}}.


}
\value{
A list. It includes elements: \code{ann} and \code{subGraphList}.

\code{ann} is a list. Each element of the list is another list. It includes the following elements: 'pathwayId', 'pathwayName', 'annMoleculeList', 'annMoleculeNumber', 'annBgMoleculeList', 'annBgNumber', 'moleculeNumber', 'bgNumber', 'pvalue', and 'fdr'. They correspond to subpathway identifier, pathway name, the submitted (interesting) molecules annotated to a pathway, numbers of submitted molecules annotated to a pathway, the background molecules annotated to a pathway, numbers of background molecules annotated to a pathway, numbers of submitted molecules, numbers of background molecules, p-value of the hypergeometric test, and Benjamini-Hochberg fdr values.

\code{subGraphList} is a graph list of subpathways.


To visualize and save the results, the \code{list} can be converted to the \code{data.frame} by the function \code{\link{printGraph}}. 

Note that \code{moleculeList} submitted by users must be a 'character' vector.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{plotAnnGraph}},\code{\link{printGraph}},\code{\link{identifyGraph}},\code{\link{getLocSubGraph}}}

\examples{
\dontrun{

##########identify metabolic subpathways by using the function SubpathwayGM###########
#get a set of interesting genes and metabolites
moleculeList<-getExample(geneNumber=1000,compoundNumber=100)
#annotate gene and metabolite sets to metabolic subpathways 
#and identify subpathways
reGM<-SubpathwayGM(moleculeList,n=5,s=5)
#convert ann to data.frame
result<-printGraph(reGM$ann)
#print the results to screen
result[1:10,]

#Each row of the result (data.frame) is a subpathway. Its columns include pathwayId, 
#pathwayName, annMoleculeRatio, annBgRatio, pvalue, and fdr. They correspond to subpathway
#identifier, pathway name, the ratio of the annotated interesting molecules, the ratio 
#of the annotated background, p-value of the hypergeometric test, and Benjamini-Hochberg
#fdr values. For annMoleculeRatio, 29/1100 means that 29 molecules in 1100 interesting molecules
#are annotated to the subpathway. For annBgRatio, 67/25051 means that 67 molecules 
#in 25051 background molecules are annotated to the subpathway. 

result1<-printGraph(reGM$ann,detail=TRUE)
##write the annotation results to tab delimited file. 
write.table(result1,file="result1.txt",row.names=FALSE,sep="\t")


#visualize subpathways in R and KEGG web site
plotAnnGraph("path:00010_1",reGM$subGraphList,reGM$ann,displayInR=TRUE,gotoKEGG=TRUE)


##Case: use the function SubpathwayGM to identify metabolic subpathways associated with colorectal cancer
#read differential genes and metabolites in colorectal cancer from files
path1<-paste(system.file(package="iSubpathwayMiner"),"/localdata/crc_diff_gene.txt",sep="")
geneList<-as.character(read.table(path1,sep="\t")[[1]])
path2<-paste(system.file(package="iSubpathwayMiner"),"/localdata/crc_diff_metabolite.txt",sep="")
metaboliteList<-as.character(read.table(path2,sep="\t")[[1]])
moleculeList<-c(geneList,metaboliteList)
#identify metabolic subpathways 
reGM<-SubpathwayGM(moleculeList,n=5,s=5)
result<-printGraph(reGM$ann)
#print the significant subpathways to screen
result[which(result[,"pvalue"]<0.01),]

result1<-printGraph(reGM$ann,detail=TRUE)
#write the identification result of subpathways to tab delimited file. 
write.table(result1,file="result1.txt",row.names=FALSE,sep="\t")
}
}
\keyword{ file }

