\name{initializeK2ri}
\alias{initializeK2ri}
\title{initialize environment k2ri}
\description{
  initialize environment k2ri.
}
\usage{
  initializeK2ri()
}
\details{
The environment variable k2ri save many information. We can use the function \code{\link{ls}} to see the variable and use
 \code{ls(k2ri)} to see information in it, which include \code{gene2ec}, \code{gene2ko}, \code{metabolicEC}, \code{metabolicKO}, \code{nonMetabolicKO}, etc. We can use the function \code{\link{get}} to obtain one of them.  
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{saveK2ri}},\code{\link{loadK2ri}}}

\examples{
  
  initializeK2ri()
  
  ##see whether k2ri is exist in R or not.
  ls()
  #if exist, we will see it as follows:
  #[1] "k2ri"
  
  ##see information in environment variable k2ri
  ls(k2ri)
  
  ##obtain one variable in environment variable k2ri
  metabolicEC<-get("metabolicEC",envir=k2ri)
}
\keyword{ file }

