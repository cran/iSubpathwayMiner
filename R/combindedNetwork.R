
##################################ECCO or KOCO#####################################################
#############metabolic#############################################
#############################
getMetabolicECCOGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-filterNode(getMetabolicGraph(metabolicEC),nodeType=c("map"))
	 return (graphList)
}
#############################
getMetabolicECCOUGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-filterNode(getUGraph(getMetabolicGraph(metabolicEC)),nodeType=c("map"))
	 return (graphList)
}
#############################
getMetabolicECCOEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getMetabolicGraph(metabolicEC),nodeType=c("map"))))
	 return (graphList)
}
#############################
getMetabolicECCOUEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicEC<-get("metabolicEC",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getUGraph(getMetabolicGraph(metabolicEC)),nodeType=c("map"))))
	 return (graphList)
}



#############non-metabolic################
##############################################################################
getNonMetabolicKOCOGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-filterNode(getNonMetabolicGraph(nonMetabolicKO))
	 return (graphList)
}
#############################
getNonMetabolicKOCOUGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-filterNode(getUGraph(getNonMetabolicGraph(nonMetabolicKO)))
	 return (graphList)
}
#############################
getNonMetabolicKOCOEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getNonMetabolicGraph(nonMetabolicKO),nodeType=c("map"))))
	 return (graphList)
}
#############################
getNonMetabolicKOCOUEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getUGraph(getNonMetabolicGraph(nonMetabolicKO)),nodeType=c("map"))))
	 return (graphList)
}

#############metabolic################
#############################
getMetabolicKOCOGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-filterNode(getMetabolicGraph(metabolicKO),nodeType=c("map"))
	 return (graphList)
}
#############################
getMetabolicKOCOUGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-filterNode(getUGraph(getMetabolicGraph(metabolicKO)),nodeType=c("map"))
	 return (graphList)
}
#############################
getMetabolicKOCOEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getMetabolicGraph(metabolicKO),nodeType=c("map"))))
	 return (graphList)
}
#############################
getMetabolicKOCOUEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicKO<-get("metabolicKO",envir=k2ri)
     graphList<-mergeNode(expandNode(filterNode(getUGraph(getMetabolicGraph(metabolicKO)),nodeType=c("map"))))
	 return (graphList)
}














#############################################KO-KO or EC-EC####################################################
#############metabolic#############################################
getMetabolicECECGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getMetabolicGraph(metabolicEC),nodeType=c("map")),nodeType="geneProduct")
	 return (graphList)
}
#################################################################################################
getMetabolicECECUGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicEC)),nodeType=c("map")),nodeType="geneProduct")
	 return (graphList)
}
#################################################################################################
getMetabolicECECEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getMetabolicGraph(metabolicEC),nodeType=c("map")),nodeType="geneProduct")))
	 return (graphList)
}
#################################################################################################
getMetabolicECECUEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicEC)),nodeType=c("map")),nodeType="geneProduct")))
	 return (graphList)
}

#############non-metabolic#############################################
#################################################################################################
getNonMetabolicKOKOGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getNonMetabolicGraph(nonMetabolicKO),nodeType=c("map")),nodeType="geneProduct")
	 return (graphList)
}
#################################################################################################
getNonMetabolicKOKOUGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getUGraph(getNonMetabolicGraph(nonMetabolicKO)),nodeType=c("map")),nodeType="geneProduct")
	 return (graphList)
}
#################################################################################################
getNonMetabolicKOKOEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getNonMetabolicGraph(nonMetabolicKO),nodeType=c("map")),nodeType="geneProduct")))
	 return (graphList)
}
#################################################################################################
getNonMetabolicKOKOUEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getUGraph(getNonMetabolicGraph(nonMetabolicKO)),nodeType=c("map")),nodeType="geneProduct")))
	 return (graphList)
}

#############metabolic#############################################
getMetabolicKOKOGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicKO<-get("metabolicKO",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getMetabolicGraph(metabolicKO),nodeType=c("map")),nodeType="geneProduct")
	 return (graphList)
}
#################################################################################################
getMetabolicKOKOUGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicKO<-get("metabolicKO",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicKO)),nodeType=c("map")),nodeType="geneProduct")
	 return (graphList)
}
#################################################################################################
getMetabolicKOKOEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicKO<-get("metabolicKO",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getMetabolicGraph(metabolicKO),nodeType=c("map")),nodeType="geneProduct")))
	 return (graphList)
}
#################################################################################################
getMetabolicKOKOUEMGraph<-function(){
     if(!exists("k2ri")) initializeK2ri()
     metabolicKO<-get("metabolicKO",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicKO)),nodeType=c("map")),nodeType="geneProduct")))
	 return (graphList)
}
















#############################################CO-CO####################################################
#############metabolic#############################################
getMetabolicCOCOGraph<-function(type="KO"){
     if(!exists("k2ri")) initializeK2ri()
	 if(type=="KO"){
     metabolicKO<-get("metabolicKO",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getMetabolicGraph(metabolicKO),nodeType=c("map")),nodeType="compound")
	 }
	 else if(type=="EC"){
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getMetabolicGraph(metabolicEC),nodeType=c("map")),nodeType="compound")	 
	 }
	 else{stop("type should be one of KO or EC.")}
	 return (graphList)
}
#################################################################################################
getMetabolicCOCOUGraph<-function(type="KO"){
     if(!exists("k2ri")) initializeK2ri()
	 if(type=="KO"){	 
     metabolicKO<-get("metabolicKO",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicKO)),nodeType=c("map")),nodeType="compound")
	 }
	 else if(type=="EC"){
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicEC)),nodeType=c("map")),nodeType="compound") 
	 }
	 else{stop("type should be one of KO or EC.")}	 
	 return (graphList)
}
#################################################################################################
getMetabolicCOCOEMGraph<-function(type="KO"){
     if(!exists("k2ri")) initializeK2ri()
	 if(type=="KO"){	 	 
     metabolicKO<-get("metabolicKO",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getMetabolicGraph(metabolicKO),nodeType=c("map")),nodeType="compound")))
	 }
	 else if(type=="EC"){
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getMetabolicGraph(metabolicEC),nodeType=c("map")),nodeType="compound")))	
     }	
	 else{stop("type should be one of KO or EC.")}	 
	 return (graphList)
}
#################################################################################################
getMetabolicCOCOUEMGraph<-function(type="KO"){
	 
     if(!exists("k2ri")) initializeK2ri()
	 if(type=="KO"){	 	 
     metabolicKO<-get("metabolicKO",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicKO)),nodeType=c("map")),nodeType="compound")))
	 }
	 else if(type=="EC"){
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(getUGraph(getMetabolicGraph(metabolicEC)),nodeType=c("map")),nodeType="compound")))
     }		 
	 else{stop("type should be one of KO or EC.")}	 
	 return (graphList)
}













#############################################gene-gene###################################################
#############metabolic#############################################
#################################################################################################
getMetabolicGEGEEMGraph<-function(type="KO"){
     if(!exists("k2ri")) initializeK2ri()
	 if(type=="KO"){		 
     metabolicKO<-get("metabolicKO",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getMetabolicGraph(metabolicKO)),nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
	 }	
	 else if(type=="EC"){
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getMetabolicGraph(metabolicEC)),nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
     }
	 else{stop("type should be one of KO or EC.")}	 
	 return (graphList)
}
#################################################################################################
getMetabolicGEGEUEMGraph<-function(type="KO"){
     if(!exists("k2ri")) initializeK2ri()
	 if(type=="KO"){		 
     metabolicKO<-get("metabolicKO",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getUGraph(getMetabolicGraph(metabolicKO))),nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
	 }	
	 else if(type=="EC"){
     metabolicEC<-get("metabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getUGraph(getMetabolicGraph(metabolicEC))),nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
     }	
	 else{stop("type should be one of KO or EC.")}	 
	 return (graphList)
}

#############non-metabolic#############################################
#################################################################################################
getNonMetabolicGEGEEMGraph<-function(type="KO"){
     if(!exists("k2ri")) initializeK2ri()
	 if(type=="KO"){		 
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getNonMetabolicGraph(nonMetabolicKO)),nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
	 }	
	 else if(type=="EC"){
     nonMetabolicEC<-get("nonMetabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getNonMetabolicGraph(nonMetabolicEC)),nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
     }	
	 else{stop("type should be one of KO or EC.")}	 
	 return (graphList)
}
#################################################################################################
getNonMetabolicGEGEUEMGraph<-function(type="KO"){
     if(!exists("k2ri")) initializeK2ri()

	 if(type=="KO"){		 
     nonMetabolicKO<-get("nonMetabolicKO",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getUGraph(getNonMetabolicGraph(nonMetabolicKO))),nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
	 }	
	 else if(type=="EC"){
     nonMetabolicEC<-get("nonMetabolicEC",envir=k2ri)
	 graphList<-mergeNode(expandNode(simplifyGraph(filterNode(mapNode(getUGraph(getNonMetabolicGraph(nonMetabolicEC))),nodeType=c("map","enzyme","ortholog")),nodeType="geneProduct")))
     }		 
	 else{stop("type should be one of KO or EC.")}	 
	 return (graphList)
}





