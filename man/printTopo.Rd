\name{printTopo}
\alias{printTopo}
\title{Print the results of topology-based pathway analysis of gene sets}
\description{
  Print the results of topology-based pathway analysis of gene sets.
}
\usage{
     printTopo(ann,detail=FALSE)
}
\arguments{
  \item{ann}{ A list. The results returned from the function \code{\link{identifyTopo}}.}
  \item{detail}{ A logical. }
}
\value{

A data.frame of the results. Columns include pathwayId, pathwayName, annComponentRatio, annBgRatio, propertyName, propertyValue, bgPropertyValue, pvalue, fdr, annComponentList, annBgComponentList, annComponentPropertyValueList, annBgComponentPropertyValueList. The \code{propertyValue} is average topological property values of submitted components annotated to a pathway. The \code{bgPropertyValue} is average topological property values of the background components annotated to a pathway The \code{annComponentRatio} is the ratio of the annotated components. For example, 30/1000 means that 30 components in 1000 components are annotated. The \code{annComponentList} are the annotated components. The \code{annComponentPropertyValueList} is topological property values of submitted components annotated to a pathway. Detailed information is provided in the function \code{\link{identifyTopo}}.

}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{identifyTopo}}}

\examples{
\dontrun{
#get pathway graphs with enzymes as nodes.
graphList<-getMetabolicECECGraph()
#get a set of genes
geneList<-getExample(geneNumber=1000,compoundNumber=0)
#topology-based pathway analysis
ann<-identifyTopo(geneList,graphList,type="gene",propertyName="degree")
result<-printTopo(ann)
#print a part of the result
result[1:5,]

##write the results to tab delimited file. 
write.table(result,file="result.txt",row.names=FALSE,sep="\t")

#detailed information is also outputed
result1<-printTopo(ann,detail=TRUE)
##write the results to tab delimited file. 
write.table(result1,file="result1.txt",row.names=FALSE,sep="\t")

}
}
\keyword{ file }

