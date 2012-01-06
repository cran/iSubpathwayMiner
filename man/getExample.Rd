\name{getExample}
\alias{getExample}
\title{Get a set of components}
\description{
  Get a set of components.
}
\usage{
 getExample(geneNumber=1000,compoundNumber=0)
}
\arguments{
  \item{geneNumber}{An integer. Set numbers of genes returned from the function.}
  \item{compoundNumber}{An integer. Set numbers of compounds returned from the function.}
}
\details{
The function can get a set of genes and/or compounds. It get the sets starting from the first component in component dataset. If users want to get random components, users can use the function \code{\link{sampleComponent}}.
}
\value{
A vector of character.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{sampleComponent}}}
\examples{
\dontrun{
##get an example of gene sets
geneList<-getExample(10,0)


##get an example of compound sets
compoundList<-getExample(0,10)

##get an example of component sets
componentList<-getExample(5,5)
}
}
\keyword{ file }

