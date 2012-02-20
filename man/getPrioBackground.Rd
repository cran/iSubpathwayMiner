\name{getPrioBackground}
\alias{getPrioBackground}
\title{Get the conservative background}
\description{
  Get the background.
}
\usage{
     getPrioBackground(type="gene")
}
\arguments{
  \item{type}{A character string. Should be one of "gene", "compound", and "gene_compound".}
}
\details{
Compared with result of \code{\link{getBackground}}, size of the background for getPrioBackground is usually small. The analysis method of pathways that use the background will be more conservative. For human, reference background of compounds is union of all metabolites from HMDB and KEGG human pathways. However, \code{\link{getBackground}} use all compounds in KEGG as background.
}
\value{
   A character vector of background components.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{identifyGraph}},\code{\link{getBackground}}}
\examples{
\dontrun{

##get all background genes
backgroundGene<-getPrioBackground(type="gene")

##get all background compounds
backgroundCompound<-getPrioBackground(type="compound")

##get all background compounds and genes
backgroundComponent<-getPrioBackground(type="gene_compound")

}
}
\keyword{ file }

