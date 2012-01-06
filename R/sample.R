###########################################################################
##get an example
getExample<-function(geneNumber=1000,compoundNumber=0){
   if(!exists("k2ri")) initializeK2ri()
     componentList<-character()
   	 gene2path<-get("gene2path",envir=k2ri)
	 allGene1<-getGeneFromKGene(as.character(gene2path[,1]))
     allGene2<-unique(get("compound",envir=k2ri))
	 componentList1<-character()
	 if(geneNumber<=0){}
     else{
     if(geneNumber<=length(allGene1)){
       componentList1<-allGene1[1:geneNumber]
     }
     else{
       componentList1<-allGene1
     }
	 }
	 componentList2<-character()
	 if(compoundNumber<=0){}
	 else{
     if(compoundNumber<=length(allGene2)){
       componentList2<-allGene2[1:compoundNumber]
     }
     else{
       componentList2<-allGene2
     }	
     }	 
	 componentList<-c(componentList1,componentList2)
     return(componentList)
   
}
#componentList<-getAexample(k=1000)

###########################################################################
##sample
sampleComponent<-function(geneNumber=1000,compoundNumber=0){
   if(!exists("k2ri")) initializeK2ri()
     componentList<-character()
   	 gene2path<-get("gene2path",envir=k2ri)
	 allGene1<-getGeneFromKGene(as.character(gene2path[,1]))
     allGene2<-unique(get("compound",envir=k2ri))
	 componentList1<-character()
	 if(geneNumber<=0){}
     else{
     if(geneNumber<=length(allGene1)){
       componentList1<-sample(allGene1,geneNumber)
     }
     else{
       componentList1<-allGene1
     }
	 }
	 componentList2<-character()
	 if(compoundNumber<=0){}
	 else{
     if(compoundNumber<=length(allGene2)){
       componentList2<-sample(allGene2,compoundNumber)
     }
     else{
       componentList2<-allGene2
     }	
     }	 
	 componentList<-c(componentList1,componentList2)
     return(componentList)
   
}
