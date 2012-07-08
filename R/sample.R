###########################################################################
##get an example
getExample<-function(geneNumber=1000,compoundNumber=0){
   if(!exists("k2ri")) initializeK2ri()
     moleculeList<-character()
   	 gene2path<-get("gene2path",envir=k2ri)
	 allGene1<-getGeneFromKGene(as.character(gene2path[,1]))
     allGene2<-unique(get("compound",envir=k2ri))
	 moleculeList1<-character()
	 if(geneNumber<=0){}
     else{
     if(geneNumber<=length(allGene1)){
       moleculeList1<-allGene1[1:geneNumber]
     }
     else{
       moleculeList1<-allGene1
     }
	 }
	 moleculeList2<-character()
	 if(compoundNumber<=0){}
	 else{
     if(compoundNumber<=length(allGene2)){
       moleculeList2<-allGene2[1:compoundNumber]
     }
     else{
       moleculeList2<-allGene2
     }	
     }	 
	 moleculeList<-c(moleculeList1,moleculeList2)
     return(moleculeList)
   
}
#moleculeList<-getAexample(k=1000)

###########################################################################
##sample
sampleMolecule<-function(geneNumber=1000,compoundNumber=0){
   if(!exists("k2ri")) initializeK2ri()
     moleculeList<-character()
   	 gene2path<-get("gene2path",envir=k2ri)
	 allGene1<-getGeneFromKGene(as.character(gene2path[,1]))
     allGene2<-unique(get("compound",envir=k2ri))
	 moleculeList1<-character()
	 if(geneNumber<=0){}
     else{
     if(geneNumber<=length(allGene1)){
       moleculeList1<-sample(allGene1,geneNumber)
     }
     else{
       moleculeList1<-allGene1
     }
	 }
	 moleculeList2<-character()
	 if(compoundNumber<=0){}
	 else{
     if(compoundNumber<=length(allGene2)){
       moleculeList2<-sample(allGene2,compoundNumber)
     }
     else{
       moleculeList2<-allGene2
     }	
     }	 
	 moleculeList<-c(moleculeList1,moleculeList2)
     return(moleculeList)
   
}
