\name{k2ri}
\alias{k2ri}
\alias{gene2ec}
\alias{gene2ko}
\alias{keggGene2gene}
\alias{metabolicEC}
\alias{metabolicKO}
\alias{nonMetabolicKO}
\docType{data}

\title{The variables in the environment variable k2ri of the system}
\description{
  The variables in the environment variable k2ri of the system.
}
\details{
The environment variable includes the variable \code{gene2ec}, \code{gene2ko}, \code{keggGene2gene}, \code{metabolicEC}, \code{metabolicKO}, \code{nonMetabolicKO}, etc.

The variable \code{metabolicEC}, \code{metabolicKO}, and \code{nonMetabolicKO} are obtained from \url{ftp://ftp.genome.jp/pub/kegg/xml/kgml/metabolic/ec}, \url{ftp://ftp.genome.jp/pub/kegg/xml/kgml/metabolic/ko}, and    \url{ftp://ftp.genome.jp/pub/kegg/xml/kgml/non-metabolic/ko} by using the function \code{\link{getPathway}}, \code{\link{updatePathway}}, or \code{\link{importPathway}}. The variable \code{metabolicEC} is a list of metabolic pathways based on enzyme commission (EC). The variable \code{metabolicKO} is a list of metabolic pathways based on KEGG Orthology (KO). The variable \code{nonMetabolicKO} is a list of non-metabolic pathways based on KEGG Orthology (KO).

The variable \code{gene2ec}, \code{gene2ko}, and \code{keggGene2gene} are obtained from \url{ftp://ftp.genome.jp/pub/kegg/genes/organisms} according to the current organism such as hsa (see the function \code{\link{updateOrgAndIdType}}).

}
\format{An environment variable}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{loadK2ri}},\code{\link{saveK2ri}},\code{\link{initializeK2ri}}}
\keyword{file}