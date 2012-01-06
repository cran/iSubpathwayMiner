\name{getBackground}
\alias{getBackground}
\title{Get the background}
\description{
  Get the background.
}
\usage{
     getBackground(type="gene")
}
\arguments{
  \item{type}{A character string. Should be one of "gene", "compound", and "gene_compound".}
}
\details{
The default background are obtained from the environment variable \code{\link{k2ri}}.
}
\value{
   A character vector of background components.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{identifyGraph}}}
\examples{
\dontrun{

##get all background genes
backgroundGene<-getBackground(type="gene")

##get all background compounds
backgroundCompound<-getBackground(type="compound")

##get all background compounds and genes
backgroundComponent<-getBackground(type="gene_compound")

}
}
\keyword{ file }

