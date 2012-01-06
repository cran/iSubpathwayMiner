\name{printGraph}
\alias{printGraph}
\title{Print the results of graph annotation and identification}
\description{
  Print the results of graph annotation and identification.
}
\usage{
     printGraph(ann,detail=FALSE)
}
\arguments{
  \item{ann}{ A list. The results returned from the function \code{\link{identifyGraph}}.}
  \item{detail}{ A logical. If true, gene lists from the function \code{\link{identifyGraph}}
                         are converted into strings, which are used to display and write results with genes. }
}
\value{
A data.frame of the annotation results. Columns include pathwayId, pathwayName, annComponentRatio, annBgRatio, pvalue, 'qvalue', 'lfdr', annComponentList, annBgComponentList. The \code{annComponentRatio} is the ratio of the annotated components. For example, 30/1000 means that 30 components in 1000 components are annotated. The annComponentList and annBgComponentList are the annotated components and the annotated background components. Detailed information is provided in the function \code{\link{identifyGraph}}.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{identifyGraph}},\code{\link{plotAnnGraph}}}

\examples{
\dontrun{
##identify graphs
metabolicKO<-get("metabolicKO",envir=k2ri)
g1<-getMetabolicGraph(metabolicKO)
nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
g2<-getNonMetabolicGraph(nonMetabolicKO,ambiguousEdgeDirection="bi-directed")
graphList<-c(g1,g2)
geneList<-getExample(geneNumber=1000)
ann<-identifyGraph(geneList,graphList)

##print results to screen
printGraph(ann)

##print subset of columns to visilize well
result<-printGraph(ann)
result[,c(1,2,5)]

##print subset of rows to visilize well
result[1:10,]

##
result[1:10,c(1,2,5)]


##write the results to tab delimited file. 
result<-printGraph(ann)
write.table(result,file="result.txt",row.names=FALSE,sep="\t")
##write the results with gene lists.
result<-printGraph(ann,detail=TRUE)
write.table(result,file="result.txt",row.names=FALSE,sep="\t")

}
}
\keyword{ file }

