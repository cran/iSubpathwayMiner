
getNonMetabolicGraph<-function(pathwayList,ambiguousEdgeDirection="bi-directed",ambiguousEdgeList=c("unknow","compound","hidden compound","state change","binding/association","dissociation"),simpleGraph=TRUE,splitCompoundSubType=TRUE,verbose=FALSE){
  if(!is.list(pathwayList)) stop("users must input a list type of data which include pathways")
  if(!ambiguousEdgeDirection %in% c("bi-directed","single","delete")) stop("ambiguousEdgeDirection should be one of bi-directed, single, or delete")
  
  pathwayListLength<-length(pathwayList)
  graphList<-list()
  Na<-c()#
  if(pathwayListLength>0){
     for(t in 1:pathwayListLength){
	    if(verbose==TRUE)
           print(paste("deal with the pathway ",t," in ",pathwayListLength, " pathways",sep=""))
        Na[t]<-pathwayList[[t]]$pathwayAttrs$number
        entrylength<-0
        entrylength<-length(pathwayList[[t]]$entry)

        name<-c();ID<-c();id<-c();names<-c();type<-c();reaction<-c();link<-c()
        graphics_name<-c();graphics_fgcolor<-c();graphics_bgcolor<-c()              
        graphics_type<-c();graphics_x<-c();graphics_y<-c()
        graphics_width<-c();graphics_height<-c();graphics_coords<-c()
		group<-list()
		#component_id<-as.character(component_id)

		j<-0
		k<-0
        for(i in 1:entrylength){
		  if(pathwayList[[t]]$entry[[i]]$type!="group"){
		     j<-j+1
             id[j]<-pathwayList[[t]]$entry[[i]]$id 
             names[j]<-pathwayList[[t]]$entry[[i]]$name
             type[j]<-pathwayList[[t]]$entry[[i]]$type
             reaction[j]<-pathwayList[[t]]$entry[[i]]$reaction
             link[j]<-pathwayList[[t]]$entry[[i]]$link
             graphics_name[j]<-pathwayList[[t]]$entry[[i]]$graphics$name 
             graphics_fgcolor[j]<-pathwayList[[t]]$entry[[i]]$graphics$fgcolor
             graphics_bgcolor[j]<-pathwayList[[t]]$entry[[i]]$graphics$bgcolor
             graphics_type[j]<-pathwayList[[t]]$entry[[i]]$graphics$type
             graphics_x[j]<-pathwayList[[t]]$entry[[i]]$graphics$x
             graphics_y[j]<-pathwayList[[t]]$entry[[i]]$graphics$y
             graphics_width[j]<-pathwayList[[t]]$entry[[i]]$graphics$width
             graphics_height[j]<-pathwayList[[t]]$entry[[i]]$graphics$height
             graphics_coords[j]<-pathwayList[[t]]$entry[[i]]$graphics$coords
          }else if(pathwayList[[t]]$entry[[i]]$type=="group"){
		       k<-k+1		       
			   group[[k]]<-pathwayList[[t]]$entry[[i]]$component
			   names(group)[k]<-pathwayList[[t]]$entry[[i]]$id
		  }
        }
		entrylength<-j
########################################################################################################		
#####################################
         relationlength<-0
         relaionlength<-length(pathwayList[[t]]$relation)

         entry1<-c();entry2<-c();type1<-c();subtype_name<-c();subtype_value<-c();ii<-0;
         if(!(pathwayList[[t]]$relation[[1]]$entry1=="unknow")){
             for(i in 1:relaionlength){

				 group_hit<-FALSE
			     if(length(group)>0){	 		 
				    group_left<-group[names(group) %in%  pathwayList[[t]]$relation[[i]]$entry1]
					group_right<-group[names(group) %in%  pathwayList[[t]]$relation[[i]]$entry2]
					if(length(group_left)>0&&length(group_right)>0){#
						group_hit<-TRUE
					    for(uu in seq(group_left[[1]])){
						   for(uuu in seq(group_right[[1]])){

								for(j in seq(pathwayList[[t]]$relation[[i]]$subtype)){
									
							        if(!(pathwayList[[t]]$relation[[i]]$subtype[[j]]$name %in% ambiguousEdgeList)||(ambiguousEdgeDirection %in% c("single","bi-directed"))){
								
										 if(splitCompoundSubType==TRUE&&pathwayList[[t]]$relation[[i]]$subtype[[j]]$name=="compound"){
										 	 ii<-ii+1
						                     entry1[ii]<-group_left[[1]][uu]
											 entry2[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$value
											 type1[ii]<-"unknow"
											 subtype_name[ii]<-"unknow"
											 subtype_value[ii]<-"unknow"
											 
											 ii<-ii+1
											 entry1[ii]<-entry2[ii-1]
						                     entry2[ii]<-group_right[[1]][uuu]				
											 type1[ii]<-type1[ii-1]
											 subtype_name[ii]<-subtype_name[ii-1]
											 subtype_value[ii]<-subtype_value[ii-1]	 
										 }else{
										     ii<-ii+1
						                     entry1[ii]<-group_left[[1]][uu]
						                     entry2[ii]<-group_right[[1]][uuu]
											 type1[ii]<-pathwayList[[t]]$relation[[i]]$type
											 subtype_name[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$name
											 subtype_value[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$value
										 }
									 }
									 
								     if(ambiguousEdgeDirection=="bi-directed"&&(pathwayList[[t]]$relation[[i]]$subtype[[j]]$name %in% ambiguousEdgeList)){
										 if(splitCompoundSubType==TRUE&&pathwayList[[t]]$relation[[i]]$subtype[[j]]$name=="compound"){
				   		                     ii<-ii+1
						                     entry1[ii]<-entry2[ii-2]
					    	                 entry2[ii]<-entry1[ii-2]
                                             type1[ii]<-type1[ii-2] 
                                             subtype_name[ii]<-subtype_name[ii-2]
                                             subtype_value[ii]<-subtype_value[ii-2]
											 
				   		                     ii<-ii+1											 
									         entry1[ii]<-entry2[ii-2]
					    	                 entry2[ii]<-entry1[ii-2]
                                             type1[ii]<-type1[ii-2] 
                                             subtype_name[ii]<-subtype_name[ii-2]
                                             subtype_value[ii]<-subtype_value[ii-2]
										 }else{
				   		                     ii<-ii+1
						                     entry1[ii]<-entry2[ii-1]
					    	                 entry2[ii]<-entry1[ii-1]
                                             type1[ii]<-type1[ii-1] 
                                             subtype_name[ii]<-subtype_name[ii-1]
                                             subtype_value[ii]<-subtype_value[ii-1]
										 }
										 
                                     }  
                                 }									 
                             }##for							  
						}						
				    }else if(length(group_left)>0&&length(group_right)<1){#
						group_hit<-TRUE
					    for(uu in seq(group_left[[1]])){

								for(j in seq(pathwayList[[t]]$relation[[i]]$subtype)){
									
							        if(!(pathwayList[[t]]$relation[[i]]$subtype[[j]]$name %in% ambiguousEdgeList)||(ambiguousEdgeDirection %in% c("single","bi-directed"))){
									    
										 if(splitCompoundSubType==TRUE&&pathwayList[[t]]$relation[[i]]$subtype[[j]]$name=="compound"){
										 	 ii<-ii+1
						                     entry1[ii]<-group_left[[1]][uu]
											 entry2[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$value
											 type1[ii]<-"unknow"
											 subtype_name[ii]<-"unknow"
											 subtype_value[ii]<-"unknow"
											 
											 ii<-ii+1
											 entry1[ii]<-entry2[ii-1]
						                     entry2[ii]<-pathwayList[[t]]$relation[[i]]$entry2				
											 type1[ii]<-type1[ii-1]
											 subtype_name[ii]<-subtype_name[ii-1]
											 subtype_value[ii]<-subtype_value[ii-1]	 
										 }else{
										     ii<-ii+1
						                     entry1[ii]<-group_left[[1]][uu]
						                     entry2[ii]<-pathwayList[[t]]$relation[[i]]$entry2
											 type1[ii]<-pathwayList[[t]]$relation[[i]]$type
											 subtype_name[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$name
											 subtype_value[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$value
										 }
									 }
									
								     if(ambiguousEdgeDirection=="bi-directed"&&(pathwayList[[t]]$relation[[i]]$subtype[[j]]$name %in% ambiguousEdgeList)){
										 if(splitCompoundSubType==TRUE&&pathwayList[[t]]$relation[[i]]$subtype[[j]]$name=="compound"){
				   		                     ii<-ii+1
						                     entry1[ii]<-entry2[ii-2]
					    	                 entry2[ii]<-entry1[ii-2]
                                             type1[ii]<-type1[ii-2] 
                                             subtype_name[ii]<-subtype_name[ii-2]
                                             subtype_value[ii]<-subtype_value[ii-2]
											 
				   		                     ii<-ii+1											 
									         entry1[ii]<-entry2[ii-2]
					    	                 entry2[ii]<-entry1[ii-2]
                                             type1[ii]<-type1[ii-2] 
                                             subtype_name[ii]<-subtype_name[ii-2]
                                             subtype_value[ii]<-subtype_value[ii-2]
										 }else{
				   		                     ii<-ii+1
						                     entry1[ii]<-entry2[ii-1]
					    	                 entry2[ii]<-entry1[ii-1]
                                             type1[ii]<-type1[ii-1] 
                                             subtype_name[ii]<-subtype_name[ii-1]
                                             subtype_value[ii]<-subtype_value[ii-1]
										 }
										 
                                     } 
                                }									 
                        }##	for								  						 
					}else if(length(group_left)<1&&length(group_right)>0){#
					    group_hit<-TRUE
					    for(uuu in seq(group_right[[1]])){

								for(j in seq(pathwayList[[t]]$relation[[i]]$subtype)){
									
							        if(!(pathwayList[[t]]$relation[[i]]$subtype[[j]]$name %in% ambiguousEdgeList)||(ambiguousEdgeDirection %in% c("single","bi-directed"))){
									     
										 if(splitCompoundSubType==TRUE&&pathwayList[[t]]$relation[[i]]$subtype[[j]]$name=="compound"){
										 	 ii<-ii+1
						                     entry1[ii]<-pathwayList[[t]]$relation[[i]]$entry1
											 entry2[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$value
											 type1[ii]<-"unknow"
											 subtype_name[ii]<-"unknow"
											 subtype_value[ii]<-"unknow"
											 
											 ii<-ii+1
											 entry1[ii]<-entry2[ii-1]
						                     entry2[ii]<-group_right[[1]][uuu]				
											 type1[ii]<-type1[ii-1]
											 subtype_name[ii]<-subtype_name[ii-1]
											 subtype_value[ii]<-subtype_value[ii-1]	 
										 }else{
										     ii<-ii+1
						                     entry1[ii]<-pathwayList[[t]]$relation[[i]]$entry1
						                     entry2[ii]<-group_right[[1]][uuu]
											 type1[ii]<-pathwayList[[t]]$relation[[i]]$type
											 subtype_name[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$name
											 subtype_value[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$value
										 }
									 }
									 
								     if(ambiguousEdgeDirection=="bi-directed"&&(pathwayList[[t]]$relation[[i]]$subtype[[j]]$name %in% ambiguousEdgeList)){
										 if(splitCompoundSubType==TRUE&&pathwayList[[t]]$relation[[i]]$subtype[[j]]$name=="compound"){
				   		                     ii<-ii+1
						                     entry1[ii]<-entry2[ii-2]
					    	                 entry2[ii]<-entry1[ii-2]
                                             type1[ii]<-type1[ii-2] 
                                             subtype_name[ii]<-subtype_name[ii-2]
                                             subtype_value[ii]<-subtype_value[ii-2]
											 
				   		                     ii<-ii+1											 
									         entry1[ii]<-entry2[ii-2]
					    	                 entry2[ii]<-entry1[ii-2]
                                             type1[ii]<-type1[ii-2] 
                                             subtype_name[ii]<-subtype_name[ii-2]
                                             subtype_value[ii]<-subtype_value[ii-2]
										 }else{
				   		                     ii<-ii+1
						                     entry1[ii]<-entry2[ii-1]
					    	                 entry2[ii]<-entry1[ii-1]
                                             type1[ii]<-type1[ii-1] 
                                             subtype_name[ii]<-subtype_name[ii-1]
                                             subtype_value[ii]<-subtype_value[ii-1]
										 }
										 
                                     }  
                                 }									 
                        }##	for									  
                    }
				 }#end length(group)
				 if(group_hit==FALSE){#

								for(j in seq(pathwayList[[t]]$relation[[i]]$subtype)){
									
							        if((!(pathwayList[[t]]$relation[[i]]$subtype[[j]]$name %in% ambiguousEdgeList))||(ambiguousEdgeDirection %in% c("single","bi-directed"))){
									     
										 if(splitCompoundSubType==TRUE&&pathwayList[[t]]$relation[[i]]$subtype[[j]]$name=="compound"){
										 	 ii<-ii+1
						                     entry1[ii]<-pathwayList[[t]]$relation[[i]]$entry1
											 entry2[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$value
											 type1[ii]<-"unknow"
											 subtype_name[ii]<-"unknow"
											 subtype_value[ii]<-"unknow"
											 
											 ii<-ii+1
											 entry1[ii]<-entry2[ii-1]
						                     entry2[ii]<-pathwayList[[t]]$relation[[i]]$entry2				
											 type1[ii]<-type1[ii-1]
											 subtype_name[ii]<-subtype_name[ii-1]
											 subtype_value[ii]<-subtype_value[ii-1]	 
										 }else{
										     ii<-ii+1
						                     entry1[ii]<-pathwayList[[t]]$relation[[i]]$entry1
						                     entry2[ii]<-pathwayList[[t]]$relation[[i]]$entry2
											 type1[ii]<-pathwayList[[t]]$relation[[i]]$type
											 subtype_name[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$name
											 subtype_value[ii]<-pathwayList[[t]]$relation[[i]]$subtype[[j]]$value
										 }
									 }
									 
								     if(ambiguousEdgeDirection=="bi-directed"&&(pathwayList[[t]]$relation[[i]]$subtype[[j]]$name %in% ambiguousEdgeList)){
										 if(splitCompoundSubType==TRUE&&pathwayList[[t]]$relation[[i]]$subtype[[j]]$name=="compound"){
				   		                     ii<-ii+1
						                     entry1[ii]<-entry2[ii-2]
					    	                 entry2[ii]<-entry1[ii-2]
                                             type1[ii]<-type1[ii-2] 
                                             subtype_name[ii]<-subtype_name[ii-2]
                                             subtype_value[ii]<-subtype_value[ii-2]
											 
				   		                     ii<-ii+1											 
									         entry1[ii]<-entry2[ii-2]
					    	                 entry2[ii]<-entry1[ii-2]
                                             type1[ii]<-type1[ii-2] 
                                             subtype_name[ii]<-subtype_name[ii-2]
                                             subtype_value[ii]<-subtype_value[ii-2]
										 }else{
				   		                     ii<-ii+1
						                     entry1[ii]<-entry2[ii-1]
					    	                 entry2[ii]<-entry1[ii-1]
                                             type1[ii]<-type1[ii-1] 
                                             subtype_name[ii]<-subtype_name[ii-1]
                                             subtype_value[ii]<-subtype_value[ii-1]
										 }
										 
                                     }         
                                }									 
			     }#end group_hit==FALSE              
            }#end for(i in 1:relaionlength)	
        }#end if(!(pathwayList[[t]]$relation[[1]]$entry1=="unknow"))		
########################################################################################################
###############################################################

             if(length(entry1)==0||(pathwayList[[t]]$relation[[1]]$entry1=="unknow")){
                 graphList[[t]]<-graph.empty(n=0,directed=TRUE)
                 graphList[[t]]<-add.vertices(graphList[[t]],entrylength,name=id,id=id,names=names,
				   type=type,reaction=reaction,link=link,graphics_name=graphics_name, 
                   graphics_fgcolor=graphics_fgcolor,graphics_bgcolor=graphics_bgcolor,graphics_type=graphics_type,
                   graphics_x=graphics_x,graphics_y=graphics_y,graphics_width=graphics_width,
				   graphics_height=graphics_height,graphics_coords=graphics_coords)
             }else{
                 Vertices<-data.frame(ID=id,id=id,names=names,type=type,reaction=reaction,link=link,
				   graphics_name=graphics_name,graphics_fgcolor=graphics_fgcolor,graphics_bgcolor=graphics_bgcolor,
				   graphics_type=graphics_type,graphics_x=graphics_x,graphics_y=graphics_y,
				   graphics_width=graphics_width,graphics_height=graphics_height,graphics_coords=graphics_coords)
                 Edges<-data.frame(entry1=entry1,entry2=entry2,type=type1,subtype_name=subtype_name,
				   subtype_value=subtype_value)
                 graphList[[t]]<-graph.data.frame(Edges,directed=TRUE,Vertices)				 
              } 
              graphList[[t]]<-set.graph.attribute(graphList[[t]],"name",pathwayList[[t]]$pathwayAttrs$name)
              graphList[[t]]<-set.graph.attribute(graphList[[t]],"number",pathwayList[[t]]$pathwayAttrs$number)
              graphList[[t]]<-set.graph.attribute(graphList[[t]],"org",pathwayList[[t]]$pathwayAttrs$org)
              graphList[[t]]<-set.graph.attribute(graphList[[t]],"title",pathwayList[[t]]$pathwayAttrs$title)
              graphList[[t]]<-set.graph.attribute(graphList[[t]],"image",pathwayList[[t]]$pathwayAttrs$image)
              graphList[[t]]<-set.graph.attribute(graphList[[t]],"link",pathwayList[[t]]$pathwayAttrs$link)
            
     } #end for(t in 1:pathwayListLength)
   }#end  if(pathwayListLength>0)
   names(graphList)<-Na
   if(simpleGraph==TRUE)
	  graphList<-getSimpleGraph(graphList)
   return(graphList)
}
