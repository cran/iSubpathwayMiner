####################################################################
##get KO sub-pathway annotation
identifyTopo<-function(componentList,graphList,type="gene",propertyName="degree",order="pvalue",decreasing=FALSE,degree.mode="total",loops = TRUE,betweenness.directed=TRUE,clusteringCoefficient.type="local",closeness.mode="all",locateOrg=TRUE,ignoreAmbiguousEnzyme=TRUE,alternative="two.sided",background=getBackground(type)){
      #print(Sys.time())
      if(typeof(componentList)!="character"){
	  print("warning: your componentList must be 'character' vector. Because the type of your current componentList is not correct, it has been conveted arbitrarily using the function as.character().")
	  componentList<-as.character(componentList)
	  }
      if(!exists("k2ri")) initializeK2ri()
	  graphList.length<-length(graphList)
	  if(graphList.length==0){
	     print("warning: there is no any graphs in graphList or these graphs are not available for pathway analyses.")
	  }	  
	  if(locateOrg==TRUE){
	     if(graphList.length>0){
	         gene2path<-get("gene2path",envir=k2ri)
	         org.path<-unique(as.character(gene2path[,2]))
		     org.path<-substring(org.path,nchar(org.path)-4,nchar(org.path))
		     graphList<-graphList[sapply(graphList,function(x) substring(x$number,0,5) %in% org.path)]
		 }
	  }
	  
	  componentList_all<-list()
	  componentList_all[[1]]<-componentList	
      compound_length<-length(componentList[substring(componentList,0,1)=="C"])
	  componentList_length<-length(componentList)
	  randomNumber<-0
	  if(randomNumber>0){
	    if(type=="gene"){
		     for(i in 1:randomNumber){
		     componentList_all[[i+1]]<-sampleComponent(componentList_length,0)
			 }
		}else if(type=="compound"){
		     for(i in 1:randomNumber){
		     componentList_all[[i+1]]<-sampleComponent(0,componentList_length)
			 }
		}else if (type=="gene_compound"){
		     for(i in 1:randomNumber){
		     componentList_all[[i+1]]<-sampleComponent(componentList_length-compound_length,compound_length)
			 }		
		}
	  }
      annList<-list()
      for(i in 1:length(graphList)){
	  #print(paste("number",i))
            ann<-list(pathwayId=character(),pathwayName="not known",annComponentList=character(),annComponentNumber=0,
                      annBgComponentList=character(),annBgNumber=0,componentNumber=0,bgNumber=0,propertyName="not known",
					  annComponentPropertyValueList=character(),propertyValue=0,annBgComponentPropertyValueList=character(),
					  bgPropertyValue=0,pvalue=1,qvalue=1,lfdr=1)

			ann$pathwayId<-paste("path:",graphList[[i]]$number,sep="")
	        node_names<-lapply(V(graphList[[i]])$names, function(x) unlist(strsplit(x,"[ ;]")))	
            #print(i)			
            KOList<-unique(unlist(node_names))	
		if(type=="gene"||type=="gene_compound"){	
            if(graphList[[i]]$org=="ko"){		
              graphGeneList<-getGeneFromKO(KOList) 
			}
			else if(graphList[[i]]$org=="ec"){
			  graphGeneList<-getGeneFromEnzyme(KOList,ignoreAmbiguousEnzyme=ignoreAmbiguousEnzyme)
			}
			else{
			  org_idType<-unlist(strsplit(graphList[[i]]$org,";"))
			  if(org_idType[1]==getOrgAndIdType()[1]){
			     if(length(org_idType)==2){
				    if(org_idType[2]==getOrgAndIdType()[2]){
					    graphGeneList<-KOList
					}
					else{stop(paste("graph ",i,"  error: it is not ec, ko, or current org graph.",sep=""))}
				 }
				 else{
				     graphGeneList<-getGeneFromKGene(KOList)
				 }
			  }
			  else{stop(paste("graph ",i,"  error: it is not ec, ko, or current org graph.",sep=""))}
			}			
            
			
        }	
       if(type=="compound"||type=="gene_compound"){	
            graphCompoundList<-KOList[substring(KOList,0,5)=="cpd:C"] 
			graphCompoundList<-unique(substring(graphCompoundList,5))
	   }
	   if(type=="gene_compound"){
	        graphComponentList<-c(graphGeneList,graphCompoundList)
	   }
	   else if(type=="gene"){
	        graphComponentList<-graphGeneList   
	   }
	   else if(type=="compound"){
	        graphComponentList<-graphCompoundList  	   
	   }
       #annotatedComponentList<-intersect(graphComponentList,componentList)	
	   
	   annotatedComponentList_all<-list()
	   for(k in seq(componentList_all)){
            annotatedComponentList_all[[k]]<-intersect(graphComponentList,componentList_all[[k]])	
       }
	   
	   if(length(annotatedComponentList_all[[1]])>0){
       #计算一个图中每个结点的属性值存入graph_node_propertyValue中。	   
	   graph_node_propertyValue<-0
	    if (propertyName=="degree"){
			graph_node_propertyValue<-igraph::degree(graphList[[i]],mode=degree.mode,loops=loops)
	    }else if(propertyName=="betweenness"){#betweenness
		    graph_node_propertyValue<-betweenness(graphList[[i]],directed=betweenness.directed)				 
	    }else if(propertyName=="clusteringCoefficient"){ #clustering coefficient
		    graph_node_propertyValue<-igraph::transitivity(graphList[[i]],type=clusteringCoefficient.type)
	        graph_node_propertyValue<-replace(graph_node_propertyValue,which(graph_node_propertyValue=="NaN"),0)	
	    }else if(propertyName=="closeness"){#Closeness centrality 
			graph_node_propertyValue<-closeness(graphList[[i]],mode=closeness.mode)	 
	    }	
		#get all components in a pathway graph
		annotatedBackgroundList<-intersect(graphComponentList,background)
		#calculate property value of all components
        #每一个组分可能对应多个结点。计算该组分每个结点的属性值存入node_propertyValue中。
	    #然后取平均值作为该组分的属性值存入component_degree。			 
            component_propertyValue<-0
            for(j in seq(annotatedBackgroundList)){
			     annNodeList<-character()
			     if(substring(annotatedBackgroundList[j],0,1)=="C"){
				    annNodeList<-paste("cpd:",annotatedBackgroundList[j],sep="")
				}else{
                     if(graphList[[i]]$org=="ko"){		
                         annNodeList<-getKOFromGene(annotatedBackgroundList[j]) 
			         }
			         else if(graphList[[i]]$org=="ec"){
			             annNodeList<-getEnzymeFromGene(annotatedBackgroundList[j],ignoreAmbiguousEnzyme=ignoreAmbiguousEnzyme)
			         }
			         else{
			              org_idType<-unlist(strsplit(graphList[[i]]$org,";"))
			               if(org_idType[1]==getOrgAndIdType()[1]){
			                   if(length(org_idType)==2){
				                     if(org_idType[2]==getOrgAndIdType()[2]){
					                      annNodeList<-annotatedBackgroundList[j]
					                 }
					                 else{stop(paste("graph ",i,"  error: it is not ec, ko, or current org graph.",sep=""))}
				                }
				                else{
				                    annNodeList<-getGeneFromKGene(annotatedBackgroundList[j])
				                }
			               }
			              else{stop(paste("graph ",i,"  error: it is not ec, ko, or current org graph.",sep=""))}
			        }	
				}

                 hit<-sapply(node_names, function(x) ifelse(any(x %in% annNodeList),TRUE,FALSE))
                 node_propertyValue<-0
				 node_propertyValue<-graph_node_propertyValue[hit]
				 component_propertyValue[j]<-mean(node_propertyValue)
			}
		#求真实和所有随机的情趣组分集合在该通路的平均属性值
		#然后对注释进通路内的每个情趣组分求属性值，取平均值作为组分集合在该通路中的属性值ann$propertyValue。	
	    #propertyValue_all<-0
	    propertyValue_all_list<-list()
	    for(k in seq(annotatedComponentList_all)){
	        hit_k<-annotatedBackgroundList%in%annotatedComponentList_all[[k]]
			#print(annotatedComponentList_all[[k]])
			#print(annotatedBackgroundList)
			#print(Sys.time())
			component_propertyValue_k<-component_propertyValue[hit_k]
			#propertyValue_all[k]<-mean(component_propertyValue_k)
			#
            propertyValue_all_list[[k]]<-component_propertyValue_k
            
        } 
			
            pathwayName<-graphList[[i]]$title
            if(length(pathwayName)!=0)
                ann$pathwayName<-pathwayName
            ann$annComponentList<-annotatedComponentList_all[[1]] 
            ann$annComponentNumber<-length(annotatedComponentList_all[[1]])
			

			ann$annBgComponentList<-annotatedBackgroundList
            ann$annBgNumber<-length(annotatedBackgroundList)
	        ann$componentNumber<-length(componentList)			
            ann$bgNumber<-length(background)

			
			ann$propertyName<-propertyName		
			ann$annComponentPropertyValueList<-propertyValue_all_list[[1]]
			ann$annBgComponentPropertyValueList<-component_propertyValue			
			ann$propertyValue<-mean(propertyValue_all_list[[1]])
			ann$bgPropertyValue<-mean(component_propertyValue)
			#print(i)
			#print(propertyValue_all_list[[1]])	
			#print(component_propertyValue)
			test<-wilcox.test(propertyValue_all_list[[1]],component_propertyValue,correct = FALSE,exact = FALSE,alternative=alternative)$p.value
			if(!is.na(test)){
			   ann$pvalue<-test
			}
	  }	#end if(length(annotatedComponentList_all[[1]])>0){
            annList[[i]]<-ann	  
	  #print(paste(Sys.time(),"end"))
	  }
      #names(annList)<-sapply(graphList,function(x) x$number)
	  if(length(annList)>0){
	     p_value<-sapply(annList,function(x) return(x$pvalue))
		 #print(p_value)
         fdrtool.List<-fdrtool(p_value,statistic="pvalue",plot=FALSE,verbose=FALSE)	
         #print(fdrtool.List$qval)
         for(i in seq(annList)){
            annList[[i]]$qvalue<-fdrtool.List$qval[i]
			annList[[i]]$lfdr<-fdrtool.List$lfdr[i]
         }
		 
         #names(annList)<-sapply(graphList,function(x) x$number)
         annList<-annList[sapply(annList,function(x) x$annComponentNumber>0)]
         annList<-annList[order(sapply(annList,function(x) x[[order]]),decreasing=decreasing)]   
	  }	  
      return(annList)
}
#####################################################################
printTopo<-function(ann,detail=FALSE){
	  if(detail==FALSE){
	  pathwayId<-sapply(ann,function(x) x$pathwayId)
      pathwayName<-sapply(ann,function(x) x$pathwayName)
      annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
      annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      propertyName<-sapply(ann,function(x) x$propertyName) 
      propertyValue<-sapply(ann,function(x) x$propertyValue) 	 
      bgPropertyValue<-sapply(ann,function(x) x$bgPropertyValue) 	  
      pvalue<-sapply(ann,function(x) x$pvalue)
      qvalue<-sapply(ann,function(x) x$qvalue)
	  lfdr<-sapply(ann,function(x) x$lfdr)
      #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
      #                       annBgRatio,pvalue,qvalue,lfdr))
      ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annComponentRatio=annComponentRatio,
	  annBgRatio=annBgRatio,propertyName=propertyName,propertyValue=propertyValue,bgPropertyValue=bgPropertyValue,
	  pvalue=pvalue,qvalue=qvalue,lfdr=lfdr,stringsAsFactors=FALSE)							 
	  }
	  else{	 
      pathwayId<-sapply(ann,function(x) x$pathwayId)	  
	  pathwayName<-sapply(ann,function(x) x$pathwayName)
	  annComponentList<-sapply(ann, function(x){ paste(x$annComponentList,collapse=";") })
      annBgComponentList<-sapply(ann, function(x){ paste(x$annBgComponentList,collapse=";")})
	  annComponentPropertyValueList<-sapply(ann, function(x){ paste(x$annComponentPropertyValueList,collapse=";") })
      annBgComponentPropertyValueList<-sapply(ann, function(x){ paste(x$annBgComponentPropertyValueList,collapse=";")})	  
	  annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
      annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      propertyName<-sapply(ann,function(x) x$propertyName) 
      propertyValue<-sapply(ann,function(x) x$propertyValue) 
      bgPropertyValue<-sapply(ann,function(x) x$bgPropertyValue)	  
      pvalue<-sapply(ann,function(x) x$pvalue)
      qvalue<-sapply(ann,function(x) x$qvalue)
	  lfdr<-sapply(ann,function(x) x$lfdr)
      #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
      #                       annBgRatio,pvalue,qvalue,lfdr,annComponentList,annBgComponentList))
      ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annComponentRatio=annComponentRatio,
	  annBgRatio=annBgRatio,propertyName=propertyName,propertyValue=propertyValue,bgPropertyValue=bgPropertyValue,pvalue=pvalue,
	  qvalue=qvalue,lfdr=lfdr,annComponentList=annComponentList,
	  annBgComponentList=annBgComponentList,annComponentPropertyValueList,annBgComponentPropertyValueList,stringsAsFactors=FALSE)								 
	  }	  
      return(ann.data.frame)
}