\name{getExample}
\alias{getExample}
\title{Get a set of molecules}
\description{
  Get a set of molecules.
}
\usage{
 getExample(geneNumber=1000,compoundNumber=0)
}
\arguments{
  \item{geneNumber}{An integer. Set numbers of genes returned from the function.}
  \item{compoundNumber}{An integer. Set numbers of compounds returned from the function.}
}
\details{
The function can get a set of genes and/or compounds. It get the sets starting from the first molecule in molecule dataset. If users want to get random molecules, users can use the function \code{\link{sampleMolecule}}.
}
\value{
A vector of character.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{sampleMolecule}}}
\examples{
\dontrun{
##get an example of gene sets
geneList<-getExample(10,0)


##get an example of compound sets
compoundList<-getExample(0,10)

##get an example of molecule sets
moleculeList<-getExample(5,5)
}
}
\keyword{ file }

