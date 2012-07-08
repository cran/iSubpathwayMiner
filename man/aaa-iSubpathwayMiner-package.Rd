\docType{package}
\name{iSubpathwayMiner-package}
\alias{iSubpathwayMiner-package}
\alias{iSubpathwayMiner}
\title{The iSubpathwayMiner package}
\description{iSubpathwayMiner is an R platform for graph-based construction and analysis of KEGG pathways.}
\section{Introduction}{
The main goal of iSubpathwayMiner is to provide graph-based reconstruction and analysis of pathways. (1) The system provides many strategies of converting pathways to graph models. Ten functions related to conversion from pathways to graphs are developed (see the functions \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{filterNode}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{mapNode}}, \code{\link{simplifyGraph}}, \code{\link{getUGraph}}, and \code{\link{getSimpleGraph}}). Furthermore, the combinations of these functions can get many (> 20) combined conversion strategies of pathway graphs (see the function \code{\link{getMetabolicECCOGraph}}). The system can also support topology-based pathway analysis of molecule sets (see the function \code{\link{identifyTopo}}). We develop KEGG layout style of pathway graphs in R to simulate the layout of the pathway picture in KEGG website (see the function \code{\link{plotGraph}}). Pathway graphs can also be exported to the GML format supported by Cytoscape.(2) The iSubpathwayMiner can support the annotation and identification of entire pathways based on molecule (gene and/or metabolite) sets (see the function \code{\link{identifyGraph}}). (3) The iSubpathwayMiner provides the k-clique method for identification of metabolic subpathways based on gene sets (see the function \code{\link{getKcSubiGraph}}). (4) The iSubpathwayMiner provides the Subpathway-GM method for identification of metabolic subpathways based on gene and metabolite sets (see the function \code{\link{SubpathwayGM}}). 

}

\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
