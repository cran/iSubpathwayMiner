\name{iSubpathwayMiner-internal}
\alias{getGraphics}
\alias{getEntry}
\alias{getSubtype}
\alias{getUnknowSubtype}
\alias{getRelation}
\alias{getUnknowRelation}
\alias{getReaction}
\alias{getUnknowReaction}
\alias{getSubstrate}
\alias{getProduct}
\alias{getUnknowSubstrate}
\alias{getUnknowProduct}
\alias{getEdgeLabel}
\alias{getEdgeLty}
\alias{getLayout}
\alias{fdr.est}
\alias{getOneNodePath}
\title{iSubpathwayMiner internal functions}
\description{
  These are functions internally called by other functions in the package \code{iSubpathwayMiner} and not meant to be called by the user.
}

\usage{
getGraphics(graphics)
getEntry(entry)
getSubtype(subtype)
getUnknowSubtype()
getRelation(relation)
getUnknowRelation()
getReaction(reaction)
getUnknowReaction()
getSubstrate(substrate)
getProduct(product)
getUnknowSubstrate()
getUnknowProduct()
getEdgeLabel(graph)
getEdgeLty(graph)
fdr.est(p)
getOneNodePath(current_node,other_nodes,pathway,n,all_shortest_paths_length,directed,method="shortestPaths")

}

\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}

\keyword{internal}
