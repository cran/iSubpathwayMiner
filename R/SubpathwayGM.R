SubpathwayGM<-function(moleculeList,n=5,s=5,background=getPrioBackground(type="gene_compound"),graphType="KO",...){
     reGM<-list()
	 if(graphType=="EC"){
     #convert all metabolic pathways to ECCO graphs
     graphList<-getMetabolicECCOGraph()
	 }else{
     #convert all metabolic pathways to KOCO graphs
     graphList<-getMetabolicKOCOGraph()	 
	 }
     #locate subpathways according to the interesting genes and metabolites 
     reGM$subGraphList<-getLocSubGraph(moleculeList,graphList,type="gene_compound",n=n,s=s,...)
     #annotate gene sets and identify significant subpathways
	 reGM$ann<-identifyGraph(moleculeList,reGM$subGraphList,type="gene_compound",background=background,...)
	 #data.frame
     #result<-printGraph(ann,detail=T)
     return(reGM)
}
#obtain all shortest paths between current_node and other_nodes (distance < = n+1)，
getOneNodePath<-function(current_node,other_nodes,pathway,n,all_shortest_paths_length,directed,method="shortestPaths"){
    current_node_length<-length(current_node)
	other_nodes_length<-length(other_nodes)
    if(other_nodes_length<1) warning("should have one nodes in other_nodes")
    if(current_node_length<1) warning("should have one current_node")
    shortest_path_set<-list();
	current_node<-as.numeric(current_node)
	other_nodes<-as.numeric(other_nodes)
	if(other_nodes_length>0){	
    for(i in 1:other_nodes_length){
        #if(all_shortest_paths_length[current_node+1,other_nodes[i]+1]<=n+1){	
        if(all_shortest_paths_length[current_node,other_nodes[i]]<=n+1){
		    if(method=="shortestPaths"){
                 shortest_path<-get.shortest.paths(pathway, current_node,other_nodes[i],mode="out") 
			}else{
                 shortest_path<-get.all.shortest.paths(pathway, current_node,other_nodes[i],mode="out")		
			}
			
            if(length(shortest_path)>0){
                shortest_path_length<-length(shortest_path[[1]])-1
                if(shortest_path_length<=n+1){
                     shortest_path_set<-c(shortest_path_set,shortest_path)
                }
	        }
        }
		# another directed paths
		if(directed==TRUE){
        #if(all_shortest_paths_length[other_nodes[i]+1,current_node+1]<=n+1){		
        if(all_shortest_paths_length[other_nodes[i],current_node]<=n+1){
            shortest_path<-character()		
		    if(method=="shortestPaths"){
                 shortest_path<-get.shortest.paths(pathway, other_nodes[i],current_node,mode="out") 
			}else{
                 shortest_path<-get.all.shortest.paths(pathway, other_nodes[i],current_node,mode="out")		
			}
            if(length(shortest_path)>0){
                shortest_path_length<-length(shortest_path[[1]])-1
                if(shortest_path_length<=n+1){
                     shortest_path_set<-c(shortest_path_set,shortest_path)
                }		 
            }
        }
		}
    }
	}
    return (shortest_path_set)
}

#get the enzyme-enzyme graphs
getLocSubGraph<-function(moleculeList,graphList,type="gene_compound",n=5,s=5,method="shortestPaths",ignoreAmbiguousEnzyme=TRUE){
        if(typeof(moleculeList)!="character"){
	         print("warning: your moleculeList must be 'character' vector. Because the type of your current moleculeList is not correct, it has been conveted arbitrarily using the function as.character().")
	        moleculeList<-as.character(moleculeList)
	    }
		if(method!="shortestPaths"&&method!="allShortestPaths") stop("the argument should be shortestPaths or allShortestPaths.")
        if(!exists("k2ri")) initializeK2ri()
		subgraphList<-list();kk<-0;
		graphGeneNodeList<-character();graphCompoundNodeList<-character();graphGeneCompoundNodeList<-character();
		#map id1 to id2 according to type of pathway graphs 
		beforeOrg<-graphList[[1]]$org
		if(type=="gene"||type=="gene_compound"){	
            if(beforeOrg=="ko"){		
                graphGeneNodeList<-getKOFromGene(moleculeList) 
			}
			else if(beforeOrg=="ec"){
			    graphGeneNodeList<-getEnzymeFromGene(moleculeList,ignoreAmbiguousEnzyme=ignoreAmbiguousEnzyme) 
		    }
			else if(unlist(strsplit(beforeOrg,"_"))[1]==getOrgAndIdType()[1]){

					if(is.na(unlist(strsplit(beforeOrg,"_"))[2])){
					    graphGeneNodeList<-getKGeneFromGene(moleculeList)
					}else{
					    if(unlist(strsplit(beforeOrg,"_"))[2]!=getOrgAndIdType()[2]){
					     stop("organism error")
				        }else{
					         graphGeneNodeList<-moleculeList
				        }						
					}

		    }		
            else{stop(paste("graph ",i,"  error: it is not ec, ko, or org graph.",sep=""))}
        }	
        if(type=="compound"||type=="gene_compound"){	
            graphCompoundNodeList<-paste("cpd:",moleculeList,sep="")			 
	    }		
        if(type=="gene_compound"){
			graphMoleculeNodeList<-c(graphGeneNodeList,graphCompoundNodeList)
        } 
        else if(type=="gene"){
			graphMoleculeNodeList<-graphGeneNodeList      
        }
        else if(type=="compound"){
			 graphMoleculeNodeList<-graphCompoundNodeList      
        }		
        ###################################3	
        if(length(graphList)>0){		
        for(i in 1:length(graphList)){
		    currentOrg<-graphList[[i]]$org
			if(currentOrg!=beforeOrg){
           		if(type=="gene"||type=="gene_compound"){	
                     if(currentOrg=="ko"){		
                         graphGeneNodeList<-getKOFromGene(moleculeList) 
			        }
			        else if(currentOrg=="ec"){
			             graphGeneNodeList<-getEnzymeFromGene(moleculeList,ignoreAmbiguousEnzyme=ignoreAmbiguousEnzyme) 
		            }	
			        else if(unlist(strsplit(currentOrg,"_"))[1]==getOrgAndIdType()[1]){
					if(is.na(unlist(strsplit(beforeOrg,"_"))[2])){
					    graphGeneNodeList<-getKGeneFromGene(moleculeList)
					}else{
					    if(unlist(strsplit(beforeOrg,"_"))[2]!=getOrgAndIdType()[2]){
					     stop("organism error")
				        }else{
					         graphGeneNodeList<-moleculeList
				        }						
					}
			        }					
                    else{stop(paste("graph ",i,"  error: it is not ec, ko, or org graph.",sep=""))}
                }	
                if(type=="compound"||type=="gene_compound"){	
                     graphCompoundNodeList<-paste("cpd:",moleculeList,sep="")			 
	            }		
                if(type=="gene_compound"){
			         graphMoleculeNodeList<-c(graphGeneNodeList,graphCompoundNodeList)
                } 
                else if(type=="gene"){
		           	graphMoleculeNodeList<-graphGeneNodeList      
                }
                else if(type=="compound"){
			         graphMoleculeNodeList<-graphCompoundNodeList      
                }	
                beforeOrg<-currentOrg				
			}
		    nodes<-character()
			directed<-is.directed(graphList[[i]])
			hit<-sapply(V(graphList[[i]])$names, function(x) ifelse(any(unlist(strsplit(x,"[ ;]")) %in% graphMoleculeNodeList),TRUE,FALSE))
			if(any(hit)==TRUE){
	             nodes<-as.character(V(graphList[[i]])[hit])
			}	
            used_nodes<-c() #记录被使用过的当前结点
            subpathways_list<-list()
            subpathways_number<-0
			all_shortest_paths_length<-shortest.paths(graphList[[i]],mode="out")

            while(length(used_nodes)<length(nodes)){

                 shortest_path_set_in_subpathways<-list() 
                 unused_nodes<-setdiff(nodes,used_nodes)
 
                 available_nodes<-unused_nodes[1] 
                 subpathway_nodes<-unused_nodes[1]
                 while(length(available_nodes)>0){

                     current_node<-available_nodes[1] 
                     unused_nodes<-setdiff(nodes,used_nodes)
                     other_nodes<-setdiff(unused_nodes,current_node)
			         if(length(other_nodes)<1||length(current_node)<1){
			             #warning("erro")
						 #print(other_nodes)
						#print("erro")
			        }else{
                         shortest_path_set<-getOneNodePath(current_node,other_nodes,graphList[[i]],n,all_shortest_paths_length,directed=directed,method=method)
				         new_hit_nodes<-setdiff(intersect(unlist(shortest_path_set),nodes),used_nodes)
                         subpathway_nodes<-union(subpathway_nodes,new_hit_nodes)
				         available_nodes<-union(available_nodes,new_hit_nodes)   
				         shortest_path_set_in_subpathways<-c(shortest_path_set_in_subpathways,shortest_path_set)
			        }
                    used_nodes<-union(used_nodes,current_node)

                    available_nodes<-setdiff(available_nodes,current_node)
                }
                subpathways_number<-subpathways_number+1
                subpathways<-list()
                subpathways[[1]]<-subpathway_nodes 
                subpathways[[2]]<-shortest_path_set_in_subpathways
                subpathways_list[[subpathways_number]]<-subpathways
            }
            for(k in seq(subpathways_list)){
                 entry1<-c();entry2<-c();
	             if(length(subpathways_list[[k]][[2]])>0){
					 V<-as.integer(unique(unlist(subpathways_list[[k]][[2]])))
					 if(length(V)>=s){
					 	 kk<-kk+1
					     #subgraphList[[kk]]<-subgraph(graphList[[i]],V)
						 subgraphList[[kk]]<-induced.subgraph(graphList[[i]],V)#new!
                         subgraphList[[kk]]$number<-paste(subgraphList[[kk]]$number,k,sep="_")
				         names(subgraphList)[kk]<-subgraphList[[kk]]$number	
					}
	            }else{
					 V<-as.integer(unique(subpathways_list[[k]][[1]]))
					 if(length(V)>=s){
					 	 kk<-kk+1
						 #subgraphList[[kk]]<-subgraph(graphList[[i]],V)
					     subgraphList[[kk]]<-induced.subgraph(graphList[[i]],V)#new!
                         subgraphList[[kk]]$number<-paste(subgraphList[[kk]]$number,k,sep="_")
				         names(subgraphList)[kk]<-subgraphList[[kk]]$number	
					}		 
	            }
            }
		}
		}
		
		return (subgraphList)
}
