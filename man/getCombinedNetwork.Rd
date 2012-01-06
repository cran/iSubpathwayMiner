\name{getCombinedNetwork}
\alias{getMetabolicECCOGraph}
\alias{getMetabolicECCOUGraph}
\alias{getMetabolicECCOEMGraph}
\alias{getMetabolicECCOUEMGraph}

\alias{getNonMetabolicKOCOGraph}
\alias{getNonMetabolicKOCOUGraph}
\alias{getNonMetabolicKOCOEMGraph}
\alias{getNonMetabolicKOCOUEMGraph}

\alias{getMetabolicKOCOGraph}
\alias{getMetabolicKOCOUGraph}
\alias{getMetabolicKOCOEMGraph}
\alias{getMetabolicKOCOUEMGraph}




\alias{getMetabolicECECGraph}
\alias{getMetabolicECECUGraph}
\alias{getMetabolicECECEMGraph}
\alias{getMetabolicECECUEMGraph}

\alias{getNonMetabolicKOKOGraph}
\alias{getNonMetabolicKOKOUGraph}
\alias{getNonMetabolicKOKOEMGraph}
\alias{getNonMetabolicKOKOUEMGraph}

\alias{getMetabolicKOKOGraph}
\alias{getMetabolicKOKOUGraph}
\alias{getMetabolicKOKOEMGraph}
\alias{getMetabolicKOKOUEMGraph}





\alias{getMetabolicCOCOGraph}
\alias{getMetabolicCOCOUGraph}
\alias{getMetabolicCOCOEMGraph}
\alias{getMetabolicCOCOUEMGraph}





\alias{getMetabolicGEGEEMGraph}
\alias{getMetabolicGEGEUEMGraph}

\alias{getNonMetabolicGEGEEMGraph}
\alias{getNonMetabolicGEGEUEMGraph}


\title{Get pathway graphs using the combination of graph conversion functions}

\description{
  Get pathway graphs using the combination of graph conversion functions.
}
\usage{
   getMetabolicECCOGraph()
   getMetabolicECCOUGraph()
   getMetabolicECCOEMGraph()
   getMetabolicECCOUEMGraph()
   
   getNonMetabolicKOCOGraph()
   getNonMetabolicKOCOUGraph()
   getNonMetabolicKOCOEMGraph()
   getNonMetabolicKOCOUEMGraph()
   
   getMetabolicKOCOGraph()
   getMetabolicKOCOUGraph()
   getMetabolicKOCOEMGraph()
   getMetabolicKOCOUEMGraph()  
  


   getMetabolicECECGraph()
   getMetabolicECECUGraph()
   getMetabolicECECEMGraph()
   getMetabolicECECUEMGraph()
   
   getNonMetabolicKOKOGraph()
   getNonMetabolicKOKOUGraph()
   getNonMetabolicKOKOEMGraph()
   getNonMetabolicKOKOUEMGraph()
   
   getMetabolicKOKOGraph()
   getMetabolicKOKOUGraph()
   getMetabolicKOKOEMGraph()
   getMetabolicKOKOUEMGraph() 


   
   getMetabolicCOCOGraph(type="KO")
   getMetabolicCOCOUGraph(type="KO")
   getMetabolicCOCOEMGraph(type="KO")
   getMetabolicCOCOUEMGraph(type="KO")
   
   
   
   getMetabolicGEGEEMGraph(type="KO")
   getMetabolicGEGEUEMGraph(type="KO")
   
   getNonMetabolicGEGEEMGraph(type="KO")
   getNonMetabolicGEGEUEMGraph(type="KO")
   
  
}
\arguments{
  \item{type}{  A character string. "KO" or "EC". If type="EC", pathway graphs are contructed based on enzyme reference pathway data, that is, the variable metabolicEC.}
}
\details{
We have provided some functions for converting pathways to graphs using the combination of graph conversion functions, including \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{getUGraph}}, \code{\link{mapNode}}, \code{\link{filterNode}}, \code{\link{simplifyGraph}}, \code{\link{mergeNode}}, \code{\link{getSimpleGraph}}, and \code{\link{expandNode}}. These functions provide many conversion strategies of pathway graphs, which represent current main applications.

The function \code{getMetabolicECCOGraph} get metabolic pathway graphs with enzymes and compounds as nodes. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-filterNode(getMetabolicGraph(metabolicEC),nodeType=c("map"))
}	
The function \code{getMetabolicECCOUGraph} get the corresponding undirected graphs, that is, the undirected graphs with enzymes and compounds as nodes. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-filterNode(getUGraph(getMetabolicGraph(metabolicEC)),nodeType=c("map"))
}	
The function \code{getMetabolicECCOEMGraph} can get metabolic pathway graphs with enzymes and compounds as nodes, in which each node only contains an enzyme/compound and each enzyme/compound only appears once. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getMetabolicGraph(metabolicEC),nodeType=c("map"))))
}	
The function \code{getMetabolicECCOUEMGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getUGraph(getMetabolicGraph(metabolicEC)),nodeType=c("map"))))
}


The function \code{getNonMetabolicKOCOGraph} can get non-metabolic pathway graphs with KOs and compounds as nodes. The result graphs are equal to the result of the following commands:
\preformatted{
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-filterNode(getNonMetabolicGraph(nonMetabolicKO),nodeType=c("map"))
}	
The function \code{getNonMetabolicKOCOUGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-filterNode(getUGraph(getNonMetabolicGraph(nonMetabolicKO)),nodeType=c("map"))
}	
The function \code{getNonMetabolicKOCOEMGraph} can get non-metabolic pathway graphs with KOs and compounds as nodes. Each node only contains a KO/compound and each KO/compound only appears once in the graph. The result graphs are equal to the result of the following commands:
\preformatted{
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getNonMetabolicGraph(nonMetabolicKO),nodeType=c("map"))))
}	
The function \code{getNonMetabolicKOCOUEMGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getUGraph(getNonMetabolicGraph(nonMetabolicKO)),nodeType=c("map"))))
}	


The function \code{getMetabolicKOCOGraph} get metabolic pathway graphs with KOs and compounds as nodes. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-filterNode(getMetabolicGraph(metabolicKO),nodeType=c("map"))
}	
The function \code{getMetabolicKOCOUGraph} get the corresponding undirected graphs, that is, the undirected graphs with KOs and compounds as nodes. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-filterNode(getUGraph(getMetabolicGraph(metabolicKO)),nodeType=c("map"))
}	
The function \code{getMetabolicKOCOEMGraph} can get metabolic pathway graphs with KOs and compounds as nodes, in which each node only contains a KO/compound and each KO/compound only appears once. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getMetabolicGraph(metabolicKO),nodeType=c("map"))))
}	
The function \code{getMetabolicKOCOUEMGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getUGraph(getMetabolicGraph(metabolicKO)),nodeType=c("map"))))
}


\preformatted{
 
 
 
 
 
 
 
}
	







The function \code{getMetabolicECECGraph} can get metabolic pathway graphs with enzymes as nodes and compounds as edges. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-simplifyGraph(filterNode(getMetabolicGraph(metabolicEC),nodeType=c("map")),nodeType="geneProduct")
}	
The function \code{getMetabolicECECUGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicEC)),
	                     nodeType=c("map")),nodeType="geneProduct")
}	
The function \code{getMetabolicECECEMGraph} can get metabolic pathway graphs with enzymes as nodes and compounds as edges, in which each node only contains an enzyme and each enzyme only appears once. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getMetabolicGraph(metabolicEC),
	                      nodeType=c("map")),nodeType="geneProduct")))
}	
The function \code{getMetabolicECECUEMGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicEC)),
	                     nodeType=c("map")),nodeType="geneProduct")))
}


The function \code{getNonMetabolicKOKOGraph} can get non-metabolic pathway graphs with KOs as nodes and compounds as edges. The result graphs are equal to the result of the following commands:
\preformatted{
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-simplifyGraph(filterNode(getNonMetabolicGraph(nonMetabolicKO),nodeType=c("map")),nodeType="geneProduct")
}	
The function \code{getNonMetabolicKOKOUGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-simplifyGraph(filterNode(getUGraph(getNonMetabolicGraph(nonMetabolicKO)),nodeType=c("map")),nodeType="geneProduct")
}	
The function \code{getNonMetabolicKOKOEMGraph} can get non-metabolic pathway graphs with KOs as nodes and compounds as edges, in which each node only contains a KO and each KO only appears once. The result graphs are equal to the result of the following commands:
\preformatted{
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getNonMetabolicGraph(nonMetabolicKO),
	                           nodeType=c("map")),nodeType="geneProduct")))
}	
The function \code{getNonMetabolicKOKOUEMGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getUGraph(getNonMetabolicGraph(nonMetabolicKO)),
	                        nodeType=c("map")),nodeType="geneProduct")))
}


The function \code{getMetabolicKOKOGraph} can get metabolic pathway graphs with KOs as nodes and compounds as edges. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-simplifyGraph(filterNode(getMetabolicGraph(metabolicKO),nodeType=c("map")),nodeType="geneProduct")
}	
The function \code{getMetabolicKOKOUGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicKO)),
	                          nodeType=c("map")),nodeType="geneProduct")
}	
The function \code{getMetabolicKOKOEMGraph} can get metabolic pathway graphs with KOs as nodes and compounds as edges, in which each node only contains a KO and each KO only appears once. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getMetabolicGraph(metabolicKO),
	                      nodeType=c("map")),nodeType="geneProduct")))
}	
The function \code{getMetabolicKOKOUEMGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicKO)),
	                      nodeType=c("map")),nodeType="geneProduct")))
}


\preformatted{
 
 
 
 
 
 
 
}
	



The function \code{getMetabolicCOCOGraph} can get metabolic pathway graphs with compounds as nodes and enzymes as edges.  The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-simplifyGraph(filterNode(getMetabolicGraph(metabolicKO),
	                          nodeType=c("map")),nodeType="compound")		
or			
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getMetabolicGraph(metabolicEC),
	                          nodeType=c("map")),nodeType="compound")	 
	 
}	
The function \code{getMetabolicCOCOUGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicKO)),
	                          nodeType=c("map")),nodeType="compound")
or							  
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicEC)),
	                          nodeType=c("map")),nodeType="compound")
}	
The function \code{getMetabolicCOCOEMGraph} can get metabolic pathway graphs with compounds as nodes and enzymes as edges, in which each node only contains a compound and each compound only appears once. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getMetabolicGraph(metabolicKO),
	              nodeType=c("map")),nodeType="compound")))
or				  
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getMetabolicGraph(metabolicEC),
	              nodeType=c("map")),nodeType="compound")))					  
}	
The function \code{getMetabolicCOCOUEMGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicKO)),
	                nodeType=c("map")),nodeType="compound")))
or					
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicEC)),
	               nodeType=c("map")),nodeType="compound")))					
}

\preformatted{
 
 
 
 
 
 
 
}
	

The function \code{getMetabolicGEGEEMGraph} can get metabolic pathway graphs with organism-specific genes as nodes and compounds as edges. Moreover, each node contains only a gene and each gene only appear once in the graph. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getMetabolicGraph(metabolicKO)),
	                     nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
or						 
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getMetabolicGraph(metabolicEC)),
	 nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))					 
}	
The function \code{getMetabolicGEGEUEMGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getUGraph(getMetabolicGraph(metabolicKO))),
	                  nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
or					  
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getUGraph(getMetabolicGraph(metabolicEC))),
	                   nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))			  
}
	
The function \code{getNonMetabolicGEGEEMGraph} can get non-metabolic pathway graphs with organism-specific genes as nodes and compounds as edges. Moreover, each node contains only a gene and each gene only appear once in the graph. The result graphs are equal to the result of the following commands:
\preformatted{
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getNonMetabolicGraph(nonMetabolicKO)),
	                 nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
or					 
     nonMetabolicEC<-get("nonMetabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getNonMetabolicGraph(nonMetabolicEC)),
	                 nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))	 
}	
The function \code{getNonMetabolicGEGEUEMGraph} can get the corresponding undirected graphs. The result graphs are equal to the result of the following commands:
\preformatted{
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getUGraph(getNonMetabolicGraph(nonMetabolicKO))),
	                     nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
or						 
     nonMetabolicEC<-get("nonMetabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getUGraph(getNonMetabolicGraph(nonMetabolicEC))),
	                     nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
}
	


}
\value{
  A list of directed graphs
}
\references{
Smart, A.G., Amaral, L.A. and Ottino, J.M. (2008) Cascading failure and robustness in metabolic networks. Proc Natl Acad Sci U S A, 105, 13223-13228.

Schreiber, F. (2002) High quality visualization of biochemical pathways in BioPath. In Silico Biol, 2, 59-73.

Klukas, C. and Schreiber, F. (2007) Dynamic exploration and editing of KEGG pathway diagrams. Bioinformatics, 23, 344-350.

Kanehisa, M., Goto, S., Hattori, M., Aoki-Kinoshita, K.F., Itoh, M., Kawashima, S., Katayama, T., Araki, M. and Hirakawa, M. (2006) From genomics to chemical genomics: new developments in KEGG. Nucleic Acids Res, 34, D354-357.

Goffard, N. and Weiller, G. (2007) PathExpress: a web-based tool to identify relevant pathways in gene expression data. Nucleic Acids Res, 35, W176-181.

Koyuturk, M., Grama, A. and Szpankowski, W. (2004) An efficient algorithm for detecting frequent subgraphs in biological networks. Bioinformatics, 20 Suppl 1, i200-207.

Hung, J.H., et al. (2010) Identification of Functional Modules That Correlate with Phenotypic Difference: The Influence of Network Topology. Genome Biol, 11, R23.

Xia, J. and Wishart, D.S. (2010) Metpa: A Web-Based Metabolomics Tool for Pathway Analysis and Visualization. Bioinformatics, 26, 2342-2344.

Jeong, H., et al. (2000) The Large-Scale Organization of Metabolic Networks. Nature, 407, 651-654.

Antonov, A.V., et al. (2008) Kegg Spider: Interpretation of Genomics Data in the Context of the Global Gene Metabolic Network. Genome Biol, 9, R179.

Guimera, R. and Nunes Amaral, L.A. (2005) Functional Cartography of Complex Metabolic Networks. Nature, 433, 895-900.

Draghici, S., Khatri, P., Tarca, A.L., Amin, K., Done, A., Voichita, C., Georgescu, C. and Romero, R. (2007) A systems biology approach for pathway level analysis. Genome Res, 17, 1537-1545.

Li, C., Li, X., Miao, Y., Wang, Q., Jiang, W., Xu, C., Li, J., Han, J., Zhang, F., Gong, B. et al. (2009) SubpathwayMiner: a software package for flexible identification of pathways. Nucleic Acids Res, 37, e131.

Ogata, H., Fujibuchi, W., Goto, S. and Kanehisa, M. (2000) A heuristic graph comparison algorithm and its application to detect functionally related enzyme clusters. Nucleic Acids Res, 28, 4021-4028.

Hung, J.H., et al. (2010) Identification of Functional Modules That Correlate with Phenotypic Difference: The Influence of Network Topology. Genome Biol, 11, R23.

Barabasi, A.L. and Oltvai, Z.N. (2004) Network Biology: Understanding the Cell's Functional Organization. Nature reviews, 5, 101-113.

KEGG Markup Language manual \url{http://www.genome.jp/kegg/xml/docs/}
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{ \code{\link{getPathway}}, \code{\link{getMetabolicGraph}}, \code{\link{getNonMetabolicGraph}}, \code{\link{filterNode}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{mapNode}}, \code{\link{simplifyGraph}}, \code{\link{getUGraph}}, \code{\link{getSimpleGraph}} }
\examples{
\dontrun{
graphList<-getMetabolicECCOGraph()
}
}
\keyword{ file }

