\name{getPathway}
\alias{getPathway}
\title{Convert KGML files to a list variable in R}
\description{
  Convert KGML files to a list data in R.
}
\usage{
  getPathway(path,filelist=list.files(path),verbose=FALSE)
}
\arguments{
  \item{path}{Path of KGML files}
  \item{filelist}{A character vector of names of KGML files}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{
The function convert KGML files to a list variable in R. 

The KEGG Markup Language (KGML) is an exchange format of the KEGG pathway data. In a KGML file, the pathway element is a root element. The entry element stores information about nodes of the pathway, including the attributes (id, name, type, link, and reaction), the "graphics" sub-element, the "component" sub-element. The relation element store information about relationship between gene products (or between gene products and compounds). It includes the attribute information (entry1, entry2, and type), and the "subtype" subelement that specifies more detailed information about the interaction. The reaction element stores chemical reaction between a substrate and a product. It includes the attribute information (id, name, and type), the "substrate" subelement, and the "product" subelement. Detailed information is provided in \url{http://www.genome.jp/kegg/xml/docs/}.

In KEGG, there are two fundamental controlled vocabularies for matching genes to pathways. Enzyme commission (EC) numbers  are traditionally used as an effective vocabulary for annotating genes to metabolic pathways. With the rapid development of KEGG, more and more non-metabolic pathways including genetic information processing, environmental information processing and cellular processes have been added to KEGG PATHWAY database. KEGG Orthology (KO) identifiers, which overcome limitations of enzyme nomenclature and integrate the pathway and genome information, have become a better controlled vocabulary for annotating genes to both metabolic and regulatory pathways. Therefore, KEGG has provided the KGML files of reference metabolic pathways linked to EC identifiers, reference metabolic pathways linked to KO identifiers, and reference non-metabolic pathways linked to KO identifiers. They can be obtained from KEGG ftp site (\url{ftp://ftp.genome.jp/pub/kegg/xml/kgml/metabolic/ec}, \url{ftp://ftp.genome.jp/pub/kegg/xml/kgml/metabolic/ko}, and \url{ftp://ftp.genome.jp/pub/kegg/xml/kgml/non-metabolic/ko}). In addition, the KGML file of organism-specific pathways linked to genes can also be obtained from KEGG ftp site.

The function can convert the above KGML files to a list variable in R, which is used as pathway data in our system (see \code{\link{k2ri}}. The conversion only changes data structure in order to efficiently operate data in R environment. After conversion, most of original information about pathways are not ignored although data structure changed. The list that stores pathway information will be used as the input of other functions such as \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, etc.
}

\value{
  A list. Each element of the list represent a pathway stored in another list. A pathway list includes entry, relation, reaction, etc. The function \code{\link{str}} can display the list structure. 
  
}

\references{
KEGG Markup Language manual \url{http://www.genome.jp/kegg/xml/docs/}
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}

\seealso{\code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{updatePathway}}, \code{\link{importPathway}}}

\examples{
\dontrun{
#Convert two XML file of metabolic pathways
#to a list in R, which is used as storing pathway maps.
path<-paste(system.file(package="iSubpathwayMiner"),"/localdata/kgml/metabolic/ko/",sep="")
p<-getPathway(path,c("ko00010.xml","ko00020.xml"))

#display the list structure of a pathway
str(p[[1]],list.len=9,max.level=5)

}
}
\keyword{ file }

