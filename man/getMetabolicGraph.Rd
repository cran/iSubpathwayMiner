\name{getMetabolicGraph}
\alias{getMetabolicGraph}
\title{Convert metabolic pathways to graphs}
\description{
  Convert metabolic pathways to directed graphs. Nodes are compounds, enzymes (or KOs), and linked pathway maps. Edges are obtained from reactions and relations.
}
\usage{
   getMetabolicGraph(pathwayList,verbose=FALSE)
}
\arguments{
  \item{pathwayList}{A list of pathways. The result of the function \code{\link{getPathway}}.}
  \item{verbose}{A logical. If TRUE, the additional diagnostics are printed.}
}
\details{
The function converts metabolic pathways to directed graphs. The returned graphs mainly contain three types of nodes: compounds, enzymes (or KOs), and linked pathway maps. Edges are obtained from reactions and relations. Detailed information about nodes, reactions, relations is provided in \url{http://www.genome.jp/kegg/xml/docs/}. In particular, if a compound participates in a reaction as a substrate or product, a directed edge connects the compound node to the reaction node (enzymes or KOs). That is, substrates of a reaction are connected to the corresponding enzyme nodes and these enzyme nodes are connected to products. For a substrate, they are directed toward the reaction node. For a product, the reaction node is directed toward them. Reversible reactions have two times more edges than irreversible reactions. In addition, a compound and a linked map may be connected by an edge if they have relationships get from relation element of the KGML file. Other information such as node attribute, pathway attribute (general information about a pathway. e.g., pathway name), etc. is converted to attribute of graphs in igraph.

We try to use the function to maintain all information of metabolic pathway maps, which can get from the KGML files, even coordinate of nodes in the drawn KEGG pathway maps. Therefore, we expect the returned graphs to be displayed to be similar with the original KEGG pathway maps. In fact, this is difficult. For example, the KGML files don't contain all pathway informations in KEGG pathway maps (pictures). KEGG pathway maps is hyper-graph consisting of hyper-edges that may connect more than two nodes. It is difficult for pathways to equally match to graph model in R. For example, if a reaction have two substrates and a product, there are not equally graph used to represent the relation. In addition, the device and functions in R cause some information in pathways not to be simulated or displayed. Although these limitations cause our current methods not to simulate completely pathway maps in R, pathway maps are simulated well by using the function \code{getMetabolicGraph} and \code{\link{plotGraph}}.

The fucntion is fit to be used to get a detailed data of pathway maps under graph model in R, and as the input data of displaying pathway maps through the function \code{\link{plotGraph}}. However, many applications related to analysis of pathways need to convert pathway maps to other forms. For example, some nodes such as linked map nodes may need to be ignored when focusing on compounds or (and) enzymes. When we only focus on enzymes and relationships between them, graphs need to be simplified to graphs with only enzymes as nodes. Therefore, we provide many functions to meed users' potential needs. Currently, these functions include \code{\link{filterNode}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{mapNode}}, \code{\link{simplifyGraph}}, \code{\link{getUGraph}}, and \code{\link{getSimpleGraph}}.

Since the returned results are graphs, users can use all functions in \code{\link{igraph}} package to analyze these graphs. For example, the function \code{\link{V}} is used to get nodes of a graph; A vertex sequence can be indexed with the \sQuote{\code{$}} operator to select (or modify) the attributes of a subset of vertices.

The argument \code{pathwayList} is the returned value of the function \code{\link{getPathway}}.

}
\value{
  A list of directed graphs
}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\references{
Smart, A.G., Amaral, L.A. and Ottino, J.M. (2008) Cascading failure and robustness in metabolic networks. Proc Natl Acad Sci U S A, 105, 13223-13228.

Schreiber, F. (2002) High quality visualization of biochemical pathways in BioPath. In Silico Biol, 2, 59-73.

Klukas, C. and Schreiber, F. (2007) Dynamic exploration and editing of KEGG pathway diagrams. Bioinformatics, 23, 344-350.

Kanehisa, M., Goto, S., Hattori, M., Aoki-Kinoshita, K.F., Itoh, M., Kawashima, S., Katayama, T., Araki, M. and Hirakawa, M. (2006) From genomics to chemical genomics: new developments in KEGG. Nucleic Acids Res, 34, D354-357.

Goffard, N. and Weiller, G. (2007) PathExpress: a web-based tool to identify relevant pathways in gene expression data. Nucleic Acids Res, 35, W176-181.

Koyuturk, M., Grama, A. and Szpankowski, W. (2004) An efficient algorithm for detecting frequent subgraphs in biological networks. Bioinformatics, 20 Suppl 1, i200-207.
}

\seealso{ \code{\link{getPathway}}, \code{\link{getNonMetabolicGraph}}, \code{\link{filterNode}}, \code{\link{expandNode}}, \code{\link{mergeNode}}, \code{\link{mapNode}}, \code{\link{simplifyGraph}}, \code{\link{getUGraph}}, \code{\link{getSimpleGraph}} }
\examples{
\dontrun{
##Convert metabolic pathways to graphs.
path<-paste(system.file(package="iSubpathwayMiner"),"/localdata/kgml/metabolic/ec/",sep="")
p<-getPathway(path,c("ec00010.xml","ec00020.xml"))
g<-getMetabolicGraph(p)

#show graph information including node number, 
#edge number, and list of all attribute names.
summary(g[[1]])
#print the same information, and also list the edges, 
#and graph, node and/or edge attributes.
print(g[["00010"]],v=TRUE,e=TRUE,g=TRUE)

#display a subgraph with 10 nodes.
sg<-subgraph(g[[1]],V(g[[1]])[1:10])
print(sg,g=TRUE,v=TRUE,e=TRUE)


sapply(c("id","names","type"),function(x) get.vertex.attribute(g[[1]],x))


#get a graph attribute
g[[1]]$title

#get names attribute of all nodes
V(g[[1]])$names

#get names attribute of the second node
V(g[[1]])[2]$names

#another method to get names attribute of nodes
get.vertex.attribute(g[[1]],"names",2)

##get node indexes whose atrribute meet the specific conditions
#get node indexes with name attribute="ec:4.1.2.13","ec:1.2.1.59", or "ec:1.2.1.12"
index<-V(g[[1]])[V(g[[1]])$names 
\%in\% c("ec:4.1.2.13","ec:1.2.1.59","ec:1.2.1.12")]
#see index value
as.integer(index)

#note that as.integer(index) return index value.
#however, index return entry ids of KGML
index

##Convert all metabolic pathways to graphs.
metabolicKO<-get("metabolicKO",envir=k2ri)
g<-getMetabolicGraph(metabolicKO["ko00130.xml"])
print(g[[1]],v=TRUE,e=TRUE,g=TRUE)


##Convert all metabolic pathways to graphs.
metabolicKO<-get("metabolicKO",envir=k2ri)
g<-getMetabolicGraph(metabolicKO)
print(g[[1]],v=TRUE,e=TRUE,g=TRUE)

##Convert metabolic pathways to graphs.
path<-paste(system.file(package="iSubpathwayMiner"),"/localdata/kgml/metabolic/ko/",sep="")
p1<-getPathway(path,c("ko00010.xml","ko00020.xml"))
g1<-getMetabolicGraph(p1)

#get node indexes with names attribute including "ko:K01623","ko:K00134", or "path:ko00640"
#A node may include multiple KOs that are splited by " " or ";"
#We can identify nodes as follows:
index2<-V(g1[[1]])[sapply(V(g1[[1]])$names,function(x) any(unlist(strsplit(x,"[ ;]"))
\%in\% c("ko:K01623","ko:K00134","ko:KO00150")))]
#display them
index2


##################export a pathway graph ###########################
#The function write.graph can export a pathway graph to foreign file formats.
#the GML format supported by Cytoscape.
write.graph(g1[[1]], "ec00010.txt", "gml")

}
}
\keyword{ file }

