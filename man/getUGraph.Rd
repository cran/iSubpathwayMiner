\name{getUGraph}
\alias{getUGraph}
\title{Convert pathways to undirected graphs}
\description{
  Convert pathways to undirected graphs.
}
\usage{
   getUGraph(graphList,simpleGraph=TRUE)
}
\arguments{
  \item{graphList}{A list of graph class.}
  \item{simpleGraph}{A logical. Get simple graphs}
}
\details{
The function can convert directed graphs to undirected graphs.

Note that the kind of conversion may cause multiple edges. Setting the argument \code{simpleGraph} value as TRUE can merge these edges by calling the function \code{\link{getSimpleGraph}}, which will get simple graphs. 
}

\value{
  A list of undirected graphs
}

\references{
Ogata, H., Fujibuchi, W., Goto, S. and Kanehisa, M. (2000) A heuristic graph comparison algorithm and its application to detect functionally related enzyme clusters. Nucleic Acids Res, 28, 4021-4028.

Li, C., Li, X., Miao, Y., Wang, Q., Jiang, W., Xu, C., Li, J., Han, J., Zhang, F., Gong, B. et al. (2009) SubpathwayMiner: a software package for flexible identification of pathways. Nucleic Acids Res, 37, e131.

Antonov, A.V., et al. (2008) Kegg Spider: Interpretation of Genomics Data in the Context of the Global Gene Metabolic Network, Genome Biol, 9, R179.

Barabasi, A.L. and Oltvai, Z.N. (2004) Network Biology: Understanding the Cell's Functional Organization, Nature reviews, 5, 101-113.

KEGG Markup Language manual \url{http://www.genome.jp/kegg/xml/docs/}
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{ \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{filterNode}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{mapNode}}, \code{\link{simplifyGraph}}, \code{\link{getSimpleGraph}} }
\examples{
\dontrun{
##Convert metabolic pathways to graphs.
path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/metabolic/ec/",sep="")
g<-getMetabolicGraph(getPathway(path,
c("ec00010.xml","ec00020.xml")),verbose=TRUE)

g1<-getUGraph(g)

##Convert non-metabolic pathways to graphs
path1<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/non-metabolic/ko/",sep="")
g2<-getNonMetabolicGraph(getPathway(path1,c("ko04010.xml")))

g3<-getUGraph(g2)

}
}
\keyword{ file }

