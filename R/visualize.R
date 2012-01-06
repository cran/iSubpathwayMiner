
#plotGraph(graph,vertex.label=getNodeLabel)
#plotGraph(graph,vertex.label=getNodeLabel(graph,"currentId",3))
plotGraph<-function(graph,margin=0,vertex.label.cex=0.6,vertex.label.font=1,vertex.size=8,
  vertex.size2=6,edge.arrow.size=0.2,edge.arrow.width=3,vertex.label=V(graph)$graphics_name,
  vertex.shape=V(graph)$graphics_type,layout=getLayout(graph),vertex.label.color="black",
  vertex.color=V(graph)$graphics_bgcolor,vertex.frame.color="dimgray",edge.color="dimgray",
  edge.label=getEdgeLabel(graph),edge.label.cex=0.6,edge.label.color="dimgray",edge.lty=getEdgeLty(graph),
  axes=FALSE,xlab="",ylab="",sub=NULL,main=NULL,...){
    if(class(graph)!="igraph") stop("the graph should be a igraph graph.")
    if(vcount(graph)==0){
         print("the graph is an empty graph.")
    }else{	 
    vertex.shape<-replace(vertex.shape,which(vertex.shape %in% c("roundrectangle","line")),"crectangle")
    vertex.color<-replace(vertex.color,which(vertex.color %in% c("unknow","none")),"white")
    if(length(vertex.shape)==0) vertex.shape<-NULL
    if(length(vertex.color)==0) vertex.color<-NULL  
    if(length(vertex.label)==0) vertex.label<-NULL 
    if(length(layout)==0) layout<-NULL 
    if(length(edge.label)==0) edge.label<-NULL
    if((axes==FALSE)&&xlab==""&&ylab==""&&is.null(sub)&&is.null(main)){
         old.mai<-par(mai=c(0.01,0.25,0.01,0.3))
         #old.mai<-par(mai=0.01+c(0,0,0,0))
         on.exit(par(mai=old.mai), add=TRUE)
    }
    plot(graph,margin=margin,vertex.label.cex=vertex.label.cex,vertex.label.font=vertex.label.font,
	      vertex.size=vertex.size,vertex.size2=vertex.size2,
         edge.arrow.size=edge.arrow.size,edge.arrow.width=edge.arrow.width,vertex.label=vertex.label,
         vertex.shape=vertex.shape,layout=layout,vertex.label.color=vertex.label.color,
         vertex.color=vertex.color,vertex.frame.color=vertex.frame.color,edge.color=edge.color,
		 edge.label=edge.label,edge.label.cex=edge.label.cex,edge.label.color=edge.label.color,
		 edge.lty=edge.lty,axes=axes,xlab=xlab,ylab=ylab,sub=sub,main=main,...)
    }
	   
}

#node.label1<-getNodeLabel(graph,displayNumber=2)
getNodeLabel<-function(graph,type="symbol",displayNumber=1){
     
     if(displayNumber<1) stop("displayNumber should >=1")
     graphics_name<-V(graph)$graphics_name
     node.names<-V(graph)$names
	 
     current.org<-graph$org
	 if(type=="symbol"){
	     if(current.org=="ec"){
	          node.label<-lapply(node.names, function(x) getSymbolFromEnzyme(unlist(strsplit(x,"[ ;]"))))
	     }else if(current.org=="ko"){
	         node.label<-lapply(node.names, function(x) getSymbolFromKO(unlist(strsplit(x,"[ ;]"))))
	     }else if(current.org==getOrgAndIdType()[1]){
	         node.label<-lapply(node.names, function(x) getSymbolFromKGene(unlist(strsplit(x,"[ ;]"))))
	     }else{stop("It is not ec, ko, or org graph.")}
	 }else if(type=="currentId"){
	     if(current.org=="ec"){
	          node.label<-lapply(node.names, function(x) getGeneFromEnzyme(unlist(strsplit(x,"[ ;]"))))
	     }else if(current.org=="ko"){
	         node.label<-lapply(node.names, function(x) getGeneFromKO(unlist(strsplit(x,"[ ;]"))))
	     }else if(current.org==getOrgAndIdType()[1]){
	         node.label<-lapply(node.names, function(x) getGeneFromKGene(unlist(strsplit(x,"[ ;]"))))
	     }else{stop("It is not ec, ko, or org graph.")} 
	 }else{stop("type should be symbol or currentId.")}
	 
	 node.label.new<-sapply(node.label,function(x) ifelse(length(x)>displayNumber,
	                  paste(paste(x[1:displayNumber],collapse=","),"...",sep=","),paste(x,collapse=",")))
	 for(i in seq(node.label.new)){
	    if(node.label.new[i]==""){
		   node.label.new[i]<-graphics_name[i]
		}
	 }
     return(node.label.new)
}

getEdgeLabel<-function(graph){
     edge.name<-E(graph)$subtype_name
     edge.value<-E(graph)$subtype_value
     #edge.label<-E(graph)$subtype_value
     edge.label<-rep("",len=length(edge.name))
     for(i in seq(edge.name)){
         edge_i<-unlist(strsplit(edge.name[i],";"))
        if("phosphorylation" %in% edge_i){
             edge.label[i]<-paste("+p",edge.label[i],sep=" ")
        }
        if("dephosphorylation" %in% edge_i){
             edge.label[i]<-paste("-p",edge.label[i],sep=" ")
        }
        if("glycosylation"  %in% edge_i){
             edge.label[i]<-paste("+g",edge.label[i],sep=" ")
        }
        if("ubiquitination"  %in% edge_i){
             edge.label[i]<-paste("+u",edge.label[i],sep=" ")
        }
        if("methylation"  %in% edge_i){
             edge.label[i]<-paste("+m",edge.label[i],sep=" ")
        }
        if("missing interaction"  %in% edge_i){
             edge.label[i]<-paste("/",edge.label[i],sep=" ")
        }
        if("dissociation"  %in% edge_i){
             edge.label[i]<-paste("|",edge.label[i],sep=" ")
        }
        if("binding/association"  %in% edge_i){
             edge.label[i]<-paste("---",edge.label[i],sep=" ")
         }
        if("repression"  %in% edge_i){
             edge.label[i]<-paste("-e-|",edge.label[i],sep=" ")
        }
        if("expression"  %in% edge_i){
             edge.label[i]<-paste("-e->",edge.label[i],sep=" ")
        }
        if("inhibition"  %in% edge_i){
             edge.label[i]<-paste("--|",edge.label[i],sep=" ")
        }
        if("activation"  %in% edge_i){
             edge.label[i]<-paste("-->",edge.label[i],sep=" ")
        }
        if("indirect effect"  %in% edge_i){
             edge.label[i]<-paste("..>",edge.label[i],sep=" ")
        }
        if("state change"  %in% edge_i){
             edge.label[i]<-paste("...",edge.label[i],sep=" ")
        }
        if("compound" %in% edge_i){
             compound<-V(graph)[V(graph)$id==edge.value[i]]$graphics_name
	         if(length(compound)==1){
                 edge.label[i]<-paste(compound,edge.label[i],sep=" ")
	         }    
        }           
    }
    return (edge.label)
}
#04350 indirect effect,04620
getEdgeLty<-function(graph){
edge.name<-E(graph)$subtype_name
edge.lty=rep("solid",len=length(edge.name))
for(i in seq(edge.name)){
  if(edge.name[i]=="indirect effect"){
     edge.lty[i]<-"longdash"
  }else if(edge.name[i]=="state change"){
     edge.lty[i]<-"longdash"
  }
}
return(edge.lty)
}
#pathwayId<-c("path:00230","path:00010")
#plotAnnGraph(c("path:00230","path:00010"),g1,anngen,gotoKEGG=TRUE)
#plotAnnGraph(anngen[[1]][1],g1,anngen,gotoKEGG=TRUE)

#########################################################################
plotAnnGraph<-function(pathwayId,graphList,ann,gotoKEGG=FALSE,orgSpecific=TRUE,multipleCell=FALSE,displayInR=TRUE,match=TRUE,vertex.frame.color="red",...){
     url.list<-c()
	 warning_result<-FALSE
     newPathwayId<-sapply(pathwayId,function(x) unlist(strsplit(x,":"))[2])
     entireNewPathwayId<-substring(newPathwayId,0,5)
	 if(displayInR==TRUE&&multipleCell==TRUE){
	    cellnumber<-ceiling(sqrt(length(pathwayId)))
	    op <- par(mfrow=c(cellnumber, cellnumber))
	 }
	 for(i in seq(pathwayId)){
         if(match==FALSE){
             matchGraph<-graphList[[entireNewPathwayId[i]]]
        }else{
             matchGraph<-graphList[[newPathwayId[i]]]
        }
         org_idType<-unlist(strsplit(matchGraph$org,";"))
         org<-org_idType[1]
		 
         #ann[sapply(ann,function(x) ifelse(x$pathwayId==pathwayId,TRUE,FALSE))]
         annComponentList<-ann[sapply(ann,function(x) ifelse(x$pathwayId==pathwayId[i],TRUE,FALSE))][[1]]$annComponentList
         
		 KOList<-""
         if(org=="ec"){
            KOList<-getEnzymeFromGene(annComponentList)
         }else if(org=="ko"){
            KOList<-getKOFromGene(annComponentList)
         }else if(org==getOrgAndIdType()[1]){
			if(length(org_idType)==2){
				 if(org_idType[2]==getOrgAndIdType()[2]){
					 KOList<-annComponentList
				}
		        else{stop(paste("graph ",i,"  error: it is not ec, ko, or current org graph.",sep=""))}
		    }
		    else{
			     KOList<-getKGeneFromGene(annComponentList)
			}
		}
		else{stop(paste("graph ",i,"  error: it is not ec, ko, or current org graph.",sep=""))}	 
		 
		 
		 
		 compound<-unique(get("compound",envir=k2ri))
         annCompoundList<-intersect(compound,annComponentList)
		 
		if(displayInR==TRUE){
             annCompoundList1<-paste("cpd:",annCompoundList,sep="")
             componentList<-c(KOList,annCompoundList1)
             frame.color<-sapply(V(matchGraph)$names, function(x)
                 ifelse(length(intersect(unlist(strsplit(unlist(strsplit(x," ")),";")),
                 componentList))>0,vertex.frame.color,"dimgray"))
             plotGraph(matchGraph,vertex.frame.color=frame.color,...)
		 }
         if(gotoKEGG==TRUE){
		     limit.length<-250
             if(orgSpecific==TRUE){
	             org<-getOrgAndIdType()[1]
	             annGeneList<-sapply(strsplit(getKGeneFromGene(annComponentList),":"), function(x) x[2])
	             temp<- paste(c(paste(org,entireNewPathwayId[i],sep=""),annGeneList,annCompoundList),sep="",collapse="+")
                 url <- paste("http://www.genome.ad.jp/dbget-bin/show_pathway?",temp,sep="")
	             #print(url)
				 if(length(unlist(strsplit(url,"")))>limit.length){
				    warning_result<-TRUE
					url.list<-c(url.list,url)
				    url<-substring(url,0,limit.length)
					print(paste("warning: gene numbers in ", pathwayId[i], " are too large.",sep=""))
				 }
                 browseURL(url)
            }else{	  
	             temp<- paste(c(paste(org,entireNewPathwayId[i],sep=""),KOList,annCompoundList),sep="",collapse="+")
                 url <- paste("http://www.genome.ad.jp/dbget-bin/show_pathway?",temp,sep = "")
				 if(length(unlist(strsplit(url,"")))>limit.length){
				    url.list<-c(url.list,url)
					url<-substring(url,0,limit.length)
					print(paste("warning: gene numbers in ", pathwayId[i], " are too large.",sep=""))
				 }
                 browseURL(url)
	        }
        }
	}
	if(displayInR==TRUE&&multipleCell==TRUE){
	     on.exit(par(mfrow=op),add=TRUE)
	}
	if(warning_result==TRUE){
	return(url.list)
	}
}

getLayout<-function(graph){
   if(length(V(graph)$graphics_x)==0||length(V(graph)$graphics_y)==0) return (NULL)
    x_y<-c()
    graphics_x<-get.vertex.attribute(graph,"graphics_x")
    index<-which(graphics_x=="unknow")
	
    if(length(index)>1){
       temp<-as.numeric(graphics_x[which(graphics_x!="unknow")])
	   if(length(temp)<2){temp<-as.numeric(c(100,600))}
       replace_value<-seq(min(temp),max(temp),by = (max(temp)-min(temp))/(length(index)-1))
       graphics_x<-replace(graphics_x,which(graphics_x=="unknow"),replace_value)
    }else if(length(index)==1){
       temp<-as.numeric(graphics_x[which(graphics_x!="unknow")])
       graphics_x<-replace(graphics_x,which(graphics_x=="unknow"),min(temp))
    } 
    graphics_x <-as.numeric(graphics_x)
	
	graphics_y<-get.vertex.attribute(graph,"graphics_y")
    index<-which(graphics_y=="unknow")
    if(length(index)>0){
       temp<-as.numeric(graphics_y[which(graphics_y!="unknow")])
	   if(length(temp)<2){temp<-as.numeric(c(100,600))}
       graphics_y<-replace(graphics_y,which(graphics_y=="unknow"),max(temp)+100)
    } 
    graphics_y <-as.numeric(graphics_y)
	
    x_y<-as.matrix(data.frame(graphics_x=graphics_x, graphics_y=graphics_y))
    x_y[,2]<--x_y[,2]
	dimnames(x_y)<-NULL
	return (x_y)
}
