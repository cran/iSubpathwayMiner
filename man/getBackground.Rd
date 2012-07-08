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
The default background are obtained from the environment variable \code{\link{k2ri}}. For human, reference background of genes is all human genes in KEGG. Reference background of metabolites is all compound in KEGG.
}
\value{
   A character vector of background molecules.
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{identifyGraph}},\code{\link{getPrioBackground}}}
\examples{
\dontrun{

##get all background genes
backgroundGene<-getBackground(type="gene")

##get all background compounds
backgroundCompound<-getBackground(type="compound")

##get all background compounds and genes
backgroundMolecule<-getBackground(type="gene_compound")

}
}
\keyword{ file }

