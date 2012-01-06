\name{getOrgAndIdType}
\alias{getOrgAndIdType}
\title{Get the type of organism and identifier in the current environment variable}
\description{
  Get the type of organism and identifier in the current environment variable.
}
\usage{
     getOrgAndIdType()
}
\details{
  Users should ensure that the type of organism and gene identifiers accord with the return value of the function 
  \code{\link{getOrgAndIdType}}. The function can help you check the type of organism and identifiers in the current
  system.

If the return values are different from the type of your genes, you need to change them with some methods. 
Detailed information is provided in \code{\link{loadK2ri}}, \code{\link{updateOrgAndIdType}}.
}
\value{
   A character vector.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{updateOrgAndIdType}}}

\examples{
getOrgAndIdType()
}
\keyword{ file }

