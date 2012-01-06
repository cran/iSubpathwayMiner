\name{plotGraph}
\alias{plotGraph}
\title{visualize a pathway graph}
\description{
  Visualize a pathway graph.
}
\usage{
plotGraph(graph,margin=0,vertex.label.cex=0.6,vertex.label.font=1,vertex.size=8,
  vertex.size2=6,edge.arrow.size=0.2,edge.arrow.width=3,vertex.label=V(graph)$graphics_name,
  vertex.shape=V(graph)$graphics_type,layout=getLayout(graph),vertex.label.color="black",
  vertex.color=V(graph)$graphics_bgcolor,vertex.frame.color="dimgray",edge.color="dimgray",
  edge.label=getEdgeLabel(graph),edge.label.cex=0.6,edge.label.color="dimgray",edge.lty=getEdgeLty(graph),
  axes=FALSE,xlab="",ylab="",sub=NULL,main=NULL,...)
}
\arguments{
  \item{graph}{ The igraph object of a pathway graph.}
  \item{margin}{ A numeric. The value is usually between -0.5 and 0.5, which is able to zoom in or out a pathway graph. The default is 0.} 
  \item{vertex.label.cex}{ A numeric vector of node label size.}
  \item{vertex.label.font}{ A numeric vector of label font.} 
  \item{vertex.size}{ A numeric vector of Node size. See \code{\link{plot.igraph}}} 
  \item{vertex.size2}{ A numeric vector of Node size.}
  \item{edge.arrow.size}{ Edge arrow size.The default is 0.2.}
  \item{edge.arrow.width}{ Edge arrow width. The default is 3.}
  \item{vertex.label}{ A vector of node label. The default is graphics_name.}  
  \item{vertex.shape}{ A vector of node shape. The default is graphics_type.} 
  \item{layout}{ A matrix of x-y coordinates with two dims. Determine the placement of the nodes for drawing a graph. The default is the KEGG node coordinates that are originally obtained from the KGML file.}  
  \item{vertex.label.color}{ A vector of node label colors. The default is black.}    
  \item{vertex.color}{ A vector of node colors. The default is the KEGG node color.}  
  \item{vertex.frame.color}{ A vector of node frame color. The default is dimgray.}
  \item{edge.color}{ A vector of edge color. The default is dimgray.}
  \item{edge.label}{ A vector of edge label.}  
  \item{edge.label.cex}{ Edge label size.}  
  \item{edge.label.color}{ A vector of edge label color. The default is dimgray.}  
  \item{edge.lty}{ A vector of line type for the edges.}    
  \item{axes}{ A logical. whether to plot axes. The default is FALSE.}
  \item{xlab}{ A character string. The label of the horizontal axis. The default is the empty string.}
  \item{ylab}{ A character string. The label of the vertical axis. The default is the empty string.}
  \item{sub}{ A character string of subtitle. }
  \item{main}{ A character string of main title. }
  \item{...}{ The arguments passed to or from methods. See \code{\link{plot.igraph}} and see \code{\link{plot}}.}  
}
\details{
The function \code{plotGraph} is able to display a pathway graph. 

The argument \code{layout} is used to determine the placement of the nodes for drawing a graph. There are mainly two preprocessed methods to determine the placement of the nodes for drawing a pathway graph: the KEGG pathway layout and \code{\link{layout}} provided in the function \code{\link{plot.igraph}} of the \code{igraph} package. The default layout is the KEGG layout, for which the coordinates of nodes in KEGG is used to determine the placement of the nodes for drawing a graph. Therefore, the returned figure by the function may be very similar to the KEGG pathway graph when information in the pathway graph is complete relatively. The layouts provided in \code{igraph} include \code{layout.reingold.tilford}, \code{layout.random}, \code{layout.circle}, \code{layout.sphere}, \code{layout.sphere}, \code{layout.fruchterman.reingold}, \code{layout.kamada.kawai},\code{layout.spring}, \code{layout.lgl}, \code{layout.graphopt}, \code{layout.mds}, \code{layout.svd}, \code{layout.norm}, \code{layout.drl}. 

The \code{layout.reingold.tilford} generates a tree-like layout, so it is mainly for trees. The \code{layout.random} places the nodes randomly. The \code{layout.circle} places the nodes on a unit circle. The \code{layout.sphere} places the nodes (approximately) uniformly on the surface of a sphere. The \code{layout.fruchterman.reingold} uses a force-based algorithm proposed by Fruchterman and Reingoldl. The \code{layout.kamada.kawai} is another force based algorithm. The \code{layout.spring} is a spring embedder algorithm. The \code{layout.fruchterman.reingold.grid} is similar to \code{layout.fruchterman.reingold} but repelling force is calculated only between nodes that are closer to each other than a limit, so it is faster. The \code{layout.lgl} is for large connected graphs, it is similar to the layout generator of the Large Graph Layout software. The \code{layout.graphopt} is a port of the graphopt layout algorithm by Michael Schmuhl. The \code{layout.mds} uses metric multidimensional scaling for generating the coordinates. The \code{layout.svd} is a currently experimental layout function based on singular value decomposition. The \code{layout.norm} normalizes a layout. It linearly transforms each coordinate separately to fit into the given limits. The \code{layout.drl} is another force-driven layout generator, it is suitable for quite large graphs.


The function \code{\link{tkplot}} is able to edit a pathway graph. One way to visualize a pathway graph is to display the pathway graph with the function \code{tkplot}, handtune the placement of the nodes, query the coordinates by the function \code{\link{tkplot.getcoords}} and use them with the function \code{plotGraph} to visualize the graph.  

}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{plotAnnGraph}},\code{\link{plot.igraph}},\code{\link{plot}}}

\examples{
\dontrun{
#####################metabolic##########################
path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/metabolic/ec/",sep="")
p<-getPathway(path,c("ec00010.xml","ec00020.xml"))
g<-getMetabolicGraph(p)

#visualize
plotGraph(g[[1]])

#change label into symbol
plotGraph(g[[1]],vertex.label=getNodeLabel)

#change label into symbol
plotGraph(g[[1]],vertex.label=getNodeLabel(g[[1]],type="currentId"))

#add red frame to the enzyme "ec:4.1.2.13" 
vertex.frame.color<-ifelse(V(g[[1]])$names=="ec:4.1.2.13","red","dimgray")
plotGraph(g[[1]],vertex.frame.color=vertex.frame.color)
#add green label to the comound "cpd:C00111" 
vertex.label.color<-ifelse(V(g[[1]])$names=="cpd:C00111","green","dimgray")
plotGraph(g[[1]],vertex.label.color=vertex.label.color)
#change node color
vertex.color<-sapply(V(g[[1]])$type,function(x) if(x=="enzyme"){"pink"}
else if(x=="compound"){"yellow"} else{"white"})
plotGraph(g[[1]],vertex.color=vertex.color)
#change node size
size<-ifelse(V(g[[1]])$graphics_name=="Starch and sucrose metabolism",20,8)
plotGraph(g[[1]],vertex.size=size,vertex.size2=size)
#change a compound label
#font size
vertex.label.cex<-ifelse(V(g[[1]])$names=="cpd:C00036",1.0,0.6)
#italic
vertex.label.font<-ifelse(V(g[[1]])$names=="cpd:C00036",3,1)
plotGraph(g[[1]],vertex.label.cex=vertex.label.cex,vertex.label.font=vertex.label.font)
#change y coordinate of an enzyme
layout<-getLayout(g[[1]])
index<-V(g[[1]])[V(g[[1]])$names=="ec:4.1.1.32"]
layout[index+1,2]<-layout[index+1,2]+50
plotGraph(g[[1]],layout=layout)

#####################non-metabolic##########################
#For a non-metabolic pathway graph
path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/non-metabolic/ko/",sep="")
p1<-getPathway(path,c("ko04010.xml"))
g1<-getNonMetabolicGraph(p1)

#change label into symbol
plotGraph(g1[[1]],vertex.label=getNodeLabel)

#change color of edges with top 10% betweenness as red
eb<-edge.betweenness(g1[[1]])
edge.color<-ifelse(eb>quantile(eb,0.9),"red","dimgray")
plotGraph(g1[[1]],edge.color=edge.color)


vertex.frame.color<-ifelse(V(g[[1]])$names 
\%in\% c("ec:4.1.2.13","ec:1.2.1.59","ec:1.2.1.12"),"red","black")

plotGraph(g[[1]],vertex.frame.color=vertex.frame.color)



vertex.frame.color<-sapply(V(g1[[1]])$names, function(x) 
ifelse(any(unlist(strsplit(x,"[ ;]")) \%in\% c("ko:K04861","ko:K06268")),"red","dimgray"))
plotGraph(g1[[1]],vertex.frame.color=vertex.frame.color)


###############zoom in a part of pathway graph############
plotGraph(g1[[1]],vertex.label=getNodeLabel,xlim=c(-1,-0.5),ylim=c(0,1))

###############annotate genes to a graph#############
##annotate genes to a graph and the identified nodes in pathways use "red"
# color to visualize them.
#get a gene list and convert them to enzyme identifiers
path<-paste(system.file(package="iSubpathwayMiner"),"/localdata/kgml/metabolic/ec/",sep="")
g2<-getMetabolicGraph(getPathway(path,c("ec00010.xml")))
#get a gene list
geneList<-getExample(geneNumber=1000)
#get the corresponding enzyme list
ECList<-getEnzymeFromGene(geneList)
#match enzymeList to identify the annotated nodes
vertex.frame.color<-sapply(V(g2[[1]])$names, function(x)
ifelse(length(intersect(unlist(strsplit(x,"[ ;]")),
 ECList))>0,"red","dimgray"))
#visualize the pathway in which red nodes are the identified nodes
plotGraph(g2[[1]],vertex.frame.color=vertex.frame.color)

###############pathway graph layouts provided in the igraph package#############

path<-paste(system.file(package="iSubpathwayMiner"),
"/localdata/kgml/metabolic/ec/",sep="")
p<-getPathway(path,c("ec00010.xml","ec00020.xml"))
g<-getMetabolicGraph(p)


#layout.random places the nodes randomly
plotGraph(g[[1]],layout=layout.random)

#layout.circle places the nodes on a unit circle
plotGraph(g[[1]],layout=layout.circle)

#layout.sphere places the nodes (approximately) uniformly on the surface of a sphere
plotGraph(g[[1]],layout=layout.sphere)

#layout.fruchterman.reingold uses a force-based algorithm proposed by Fruchterman and Reingold
plotGraph(g[[1]],layout=layout.fruchterman.reingold)

#layout.kamada.kawai is another force based algorithm.
plotGraph(g[[1]],layout=layout.kamada.kawai)

#layout.spring is a spring embedder algorithm
plotGraph(g[[1]],layout=layout.spring)

#layout.fruchterman.reingold.grid is similar to layout.fruchterman.reingold but repelling force is calculated only 
#between nodes that are closer to each other than a limit, so it is faster. 
plotGraph(g[[1]],layout=layout.fruchterman.reingold.grid)

#layout.lgl is for large connected graphs, it is similar to the layout generator of the Large Graph Layout software 
#plotGraph(g[[1]],layout=layout.lgl)

#layout.graphopt is a port of the graphopt layout algorithm by Michael Schmuhl
plotGraph(g[[1]],layout=layout.graphopt)

##layout.reingold.tilford generates a tree-like layout, so it is mainly for trees.
plotGraph(g[[1]],layout=layout.reingold.tilford)

#layout.mds uses metric multidimensional scaling for generating the coordinates. 
plotGraph(g[[1]],layout=layout.mds)

#layout.svd is a currently experimental layout function based on singular value decomposition.
plotGraph(g[[1]],layout=layout.svd)

#layout.drl is another force-driven layout generator, it is suitable for quite large graphs
plotGraph(g[[1]],layout=layout.drl,vertex.label="")

#layout.norm normalizes a layout, it linearly transforms each coordinate separately to fit into the given limits
#plotGraph(g[[1]],layout=layout.norm)

##################edit graph ###########################
##The function tkplot is able to edit a pathway graph. 

path<-paste(system.file(package="iSubpathwayMiner"),"/localdata/kgml/metabolic/ec/",sep="")
gm<-getMetabolicGraph(getPathway(path,c("ec00010.xml","ec00020.xml")))
#display a graph in tkplot
tkp.id<-tkplot(gm[[1]],layout=getLayout(gm[[1]]),vertex.label=V(gm[[1]])$graphics_name,
vertex.color=V(gm[[1]])$graphics_bgcolor,vertex.size=8,vertex.label.cex=0.5)
#users can handtune the placement of the nodes and then obtain new 
#coordinates using the function tkplot.getcoords.
layout<-tkplot.getcoords(tkp.id)
layout<-cbind(layout[,1],-layout[,2])
#display the graph with new corrodinates from tkplot
plotGraph(gm[[1]],layout=layout)


##################export a pathway graph ###########################
#The function write.graph can export a pathway graph to foreign file formats.
#the GML format supported by Cytoscape.
write.graph(gm[[1]], "ec00010.txt", "gml")

}
}
\keyword{ file }

