####################################################################
##get KO sub-pathway annotation
identifyGraph<-function(componentList,graphList,type="gene",background=getBackground(type),
   order="pvalue",decreasing=FALSE,locateOrg=TRUE,ignoreAmbiguousEnzyme=TRUE){
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
      annList<-list()
	  if(graphList.length>0){
      for(i in 1:length(graphList)){
            ann<-list(pathwayId=character(),pathwayName="not known",annComponentList=character(),annComponentNumber=0,
                      annBgComponentList=character(),annBgNumber=0,componentNumber=0,bgNumber=0,pvalue=1,qvalue=1,lfdr=1)

			ann$pathwayId<-paste("path:",graphList[[i]]$number,sep="")
            KOList<-unique(unlist(strsplit(unlist(strsplit(V(graphList[[i]])$names," ")),";")))
            #KOList<-unique(unlist(strsplit(V(graphList[[i]])$names,"[ ;]")))			
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
       annotatedComponentList<-intersect(graphComponentList,componentList)
	   annotatedBackgroundList<-intersect(graphComponentList,background)	
            
            pathwayName<-graphList[[i]]$title
            if(length(pathwayName)!=0)
                ann$pathwayName<-pathwayName
            ann$annComponentList<-annotatedComponentList 
         
            ann$annComponentNumber<-length(annotatedComponentList)
			ann$annBgComponentList<-annotatedBackgroundList
            ann$annBgNumber<-length(annotatedBackgroundList)

            ann$componentNumber<-length(componentList)
            ann$bgNumber<-length(background)

            ann$pvalue<-1-phyper(ann$annComponentNumber-1,ann$annBgNumber,
                 ann$bgNumber-ann$annBgNumber,ann$componentNumber)
            
            annList[[i]]<-ann
      } 
	  }
	  if(length(annList)>0){
	     p_value<-sapply(annList,function(x) return(x$pvalue))
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
printGraph<-function(ann,detail=FALSE){
	  if(detail==FALSE){
	  pathwayId<-sapply(ann,function(x) x$pathwayId)
      pathwayName<-sapply(ann,function(x) x$pathwayName)
      annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
      annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      pvalue<-sapply(ann,function(x) x$pvalue)
      qvalue<-sapply(ann,function(x) x$qvalue)
	  lfdr<-sapply(ann,function(x) x$lfdr)
      #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
      #                       annBgRatio,pvalue,qvalue,lfdr))
      ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annComponentRatio=annComponentRatio,
	  annBgRatio=annBgRatio,pvalue=pvalue,qvalue=qvalue,lfdr=lfdr,stringsAsFactors=FALSE)							 
	  }
	  else{	 
      pathwayId<-sapply(ann,function(x) x$pathwayId)	  
	  pathwayName<-sapply(ann,function(x) x$pathwayName)
	  annComponentList<-sapply(ann, function(x){ paste(x$annComponentList,collapse=";") })
      annBgComponentList<-sapply(ann, function(x){ paste(x$annBgComponentList,collapse=";")})
	  annComponentRatio<-sapply(ann,function(x) paste(x$annComponentNumber,x$componentNumber,sep="/"))
      annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      pvalue<-sapply(ann,function(x) x$pvalue)
      qvalue<-sapply(ann,function(x) x$qvalue)
	  lfdr<-sapply(ann,function(x) x$lfdr)
      #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annComponentRatio,
      #                       annBgRatio,pvalue,qvalue,lfdr,annComponentList,annBgComponentList))
      ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annComponentRatio=annComponentRatio,
	  annBgRatio=annBgRatio,pvalue=pvalue,qvalue=qvalue,lfdr=lfdr,annComponentList=annComponentList,
	  annBgComponentList=annBgComponentList,stringsAsFactors=FALSE)								 
	  }
      return(ann.data.frame)
}