\name{updateTaxonomy}
\alias{updateTaxonomy}
\title{Update taxonomy}
\description{
  Update information about organism name and the three- or four-letter KEGG organism code.
}
\usage{
updateTaxonomy(path="ftp://ftp.genome.jp/pub/kegg/genes",verbose=TRUE)
}
\arguments{
  \item{path}{A character string.}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{

The function is able to update the variable \code{taxonomy} in the environment variable. The variable stores information about organism name and the three- or four-letter KEGG organism code.

}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{updateOrgAndIdType}},code{\link{getOrgIdFromOrgName}}}

\examples{
\dontrun{
##update compounds
updateTaxonomy()
}
}
\keyword{ file }

