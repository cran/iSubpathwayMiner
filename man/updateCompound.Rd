\name{updateCompound}
\alias{updateCompound}
\title{Update compound identifiers}
\description{
  update compound identifiers.
}
\usage{
updateCompound(path="ftp://ftp.genome.jp/pub/kegg/ligand",verbose=TRUE)
}
\arguments{
  \item{path}{A character string.}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{

The function is able to update the variable \code{compound} in the environment variable. 

}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{updateOrgAndIdType}}}

\examples{
\dontrun{
##update compounds
updateCompound()
}
}
\keyword{ file }

