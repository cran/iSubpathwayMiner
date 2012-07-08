\name{sampleMolecule}
\alias{sampleMolecule}
\title{Generate a random molecule set}
\description{
  Generate a random molecule set.
}
\usage{
 sampleMolecule(geneNumber=1000,compoundNumber=0)
}
\arguments{
  \item{geneNumber}{An integer. Set numbers of genes returned from the function.}
  \item{compoundNumber}{An integer. Set numbers of compounds returned from the function.}
}

\value{
A vector of character.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{getExample}}}
\examples{
\dontrun{
##get an example of gene sets
geneList<-sampleMolecule(10,0)

##get an example of compound sets
compoundList<-sampleMolecule(0,10)

##get an example of molecule sets
moleculeList<-sampleMolecule(5,5)
}
}
\keyword{ file }

