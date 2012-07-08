\name{identifyGraph}
\alias{identifyGraph}
\title{Annotate a set of molecules to pathway graphs and identify statistically significantly enriched pathway graphs}
\description{
  Annotate a set of molecules to pathway graphs and identify statistically significantly enriched pathway graphs.
}
\usage{
  identifyGraph(moleculeList,graphList,type="gene",background=getBackground(type),
   order="pvalue",decreasing=FALSE,locateOrg=TRUE,ignoreAmbiguousEnzyme=TRUE)
}
\arguments{
  \item{moleculeList}{A character vector of molecules.}
  \item{graphList}{A list. Its elements is a graph with the \code{igraph} type.} 
  \item{type}{A character string. Should be one of "gene", "compound", and "gene_compound".}   
  \item{background}{A character vector of molecules used to identify the statistically significantly enriched graphs.}
  \item{order}{A character string. Should be one of "pvalue", "fdr".}
  \item{decreasing}{A logical. Should the sort order be increasing or decreasing?}
  \item{locateOrg}{A logical. Should only the current organism pathways be considered?}
  \item{ignoreAmbiguousEnzyme}{A logical. Should ambiguous enzymes such as 2.6.1.- be ignored?}  
}
\details{

The function can annotate a set of molecules to pathway graphs and identify the statistically significantly enriched pathway graphs. Currently, the function can support input of three kinds of cellular molecule sets: gene sets, compound (metabolite) sets, and gene and compound sets at the same time. Therefore, the function can support the annotation and identification of pathways based on gene sets, compound sets , and gene and compound sets. For example, if users have interesting gene sets and interesting compound sets, the combined set of compounds and genes can be annotated at the same time to pathways. The enrichment significance of pathways can be evaluated using the hypergeometric test.  

The argument \code{type} should be one of "gene", "compound", and "gene_compound", which represent respectively to input a set of genes, compounds, "genes and compounds". Through setting the argument, the function \code{identifyGraph} can annotate gene or/and compound sets to entire pathways or subpathways and identify statistically significantly enriched entire pathways or subpathways.

If users don't set the values of the argument \code{background}, the background distribution will be obtained from the whole-genome molecules. Detailed information is provided in the function \code{\link{getBackground}}.

Note that the argument \code{ignoreAmbiguousEnzyme} is only applied to pathway graphs with enzymes as nodes. If TRUE, genes that ambiguous enzymes (e.g., "-") correspond to are ignored in analyses. See the function \code{\link{getGeneFromEnzyme}}.

When many correlated pathways (entire pathways or subpathways) are considered, a false positive discovery rate is likely to result. We also provide fdr values for reducing the false positive discovery rate.
}
\value{
A list. Each element of the list is another list. It includes the following elements: 'pathwayId', 'pathwayName', 'annMoleculeList', 'annMoleculeNumber', 'annBgMoleculeList', 'annBgNumber', 'MoleculeNumber', 'bgNumber', 'pvalue', and 'fdr'. They correspond to pathway identifier, pathway name, the submitted molecules annotated to a pathway, numbers of submitted molecules annotated to a pathway, the background molecules annotated to a pathway, numbers of background molecules annotated to a pathway, numbers of submitted molecules, numbers of background molecules, p-value of the hypergeometric test, and Benjamini-Hochberg fdr values.

The background molecules annotated to a pathway are equal to all molecules in the pathway. For example, if the submitted molecules are human genes, the background molecules annotated to a pathway are equal to all human genes in the pathway.

The numbers of background molecules are the number of all molecules. For example, if the submitted molecules are human genes, numbers of background molecules is equal to all human genes.

To visualize and save the results, the \code{list} can be converted to the \code{data.frame} by the function \code{\link{printGraph}}. 

note that \code{moleculeList} submitted by users must be a 'character' vector.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{plotAnnGraph}},\code{\link{printGraph}}}

\examples{
\dontrun{

##########get the list of pathway graphs#########
##Convert all metabolic pathways to graphs.
metabolicKO<-get("metabolicKO",envir=k2ri)
g1<-getMetabolicGraph(metabolicKO)
##Convert all non-metabolic pathways to graphs, 
#in which ambiguous edges is considered as "<-->"
nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
g2<-getNonMetabolicGraph(nonMetabolicKO,ambiguousEdgeDirection="bi-directed")
#get all metabolic and non-metabolic graphs
graphList<-c(g1,g2)


##########annotate gene sets and identify enriched metabolic and non-metabolic pathways###########
##get a set of genes
moleculeList<-getExample(geneNumber=1000)
#annotate gene sets to graphs 
#and identify significant graphs
ann<-identifyGraph(moleculeList,graphList)
#convert ann to data.frame
result<-printGraph(ann)
##write the annotation results to tab delimited file. 
write.table(result,file="result.txt",row.names=FALSE,sep="\t")

#print the results to screen
result[1:5]

result1<-printGraph(ann,detail=TRUE)
##write the annotation results to tab delimited file. 
write.table(result1,file="result1.txt",row.names=FALSE,sep="\t")


#visualize annotated graph
plotAnnGraph("path:00010",graphList,ann,gotoKEGG=FALSE)
#goto kegg, and orgSpecific=TRUE
plotAnnGraph("path:00010",graphList,ann,gotoKEGG=TRUE)
#orgSpecific=FALSE
plotAnnGraph("path:00010",graphList,ann,gotoKEGG=TRUE,orgSpecific=FALSE)



##########annotate metabolite sets and identify enriched metiablic pathways###########
#get a set of compounds
moleculeList<-getExample(geneNumber=0,compoundNumber=100)
#annotate compound sets to metabolic graphs 
#and identify significant graphs
ann<-identifyGraph(moleculeList,g1,type="compound")
#convert ann to data.frame
result<-printGraph(ann)
result[1:5]

##########annotate gene and compound sets and identify enriched metabolic pathways###########
#get a set of compounds and genes
moleculeList<-getExample(geneNumber=1000,compoundNumber=100)
#annotate gene and compound sets to metabolic graphs 
#and identify significant graphs
ann<-identifyGraph(moleculeList,g1,type="gene_compound")
#convert ann to data.frame
result<-printGraph(ann)
result[,1:5]




}
}
\keyword{ file }

