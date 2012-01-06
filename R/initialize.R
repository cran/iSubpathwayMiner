.First.lib<-function(lib, pkgname){
  #library.dynam(pkgname, pkgname, lib)
  initializeK2ri()   
}
########################################################################
##initialize data
initializeK2ri<-function(){
   data("hsa_ncbi-geneid_k2ri")
}

#sampleComponent(10,0)
########################################################################
##update org and idType
#updateOrgAndIdType("hsa","symbol")
updateOrgAndIdType<-function(org=getOrgAndIdType()[1],idType=getOrgAndIdType()[2],
  path="ftp://ftp.genome.jp/pub/kegg/genes/organisms",symbolData=TRUE,verbose=TRUE){
      if(!exists("k2ri")) initializeK2ri()
      if(verbose==TRUE){
        print("Note that the programming may be time consumming!")  
      } 	 
      #new!	  
	   if(symbolData==TRUE){
	  	  if(verbose==TRUE){
             print("download relations between gene and symbol.")
	      }	 
	     orgName<-getOrgNameFromOrgId(org)
         orgName<-tolower(orgName)
	     #org<-getOrgAndIdType()[1]
         gene2symbol1<-read.table(paste(path,"/",org,"/",orgName,".pos",sep=""),header=FALSE,
	         sep = "\t", quote="\"",fill=TRUE,stringsAsFactors=FALSE)
         gene2symbol2<-gene2symbol1[which(gene2symbol1[,2]!=""),1:2]
         gene2symbol<-data.frame(gene=paste(org,as.character(gene2symbol2[,1]),sep=":"),
         symbol=paste("symbol",sapply(gene2symbol2[,2],function(x) unlist(strsplit(x,", "))[1],USE.NAMES=FALSE),sep=":"))	  
	     assign("gene2symbol",gene2symbol,envir=k2ri)
	  }else{
	  	  if(verbose==TRUE){
             print("If symbolData==FALSE, may cause the functions associated with gene symbol to become wrong.")
	      }		     
	  }
	  if(verbose==TRUE){
        print("download relations between gene and keggGene.")
	  }	  
	  file1<-paste(path,"/",org,"/",org,"_",idType,".list",sep="")
	  if(idType!="symbol"){
          #GenekeggGene<-convertFile1List(file1)
	     keggGene2gene<-read.table(file1,header=FALSE,sep = "\t", quote="\"")
         assign("keggGene2gene",keggGene2gene,envir=k2ri)
	  }else{
	  	  orgName<-getOrgNameFromOrgId(org)
          orgName<-tolower(orgName)
	      #org<-getOrgAndIdType()[1]
         gene2symbol1<-read.table(paste(path,"/",org,"/",orgName,".pos",sep=""),header=FALSE,
	         sep = "\t", quote="\"",fill=TRUE,stringsAsFactors=FALSE)
         gene2symbol2<-gene2symbol1[which(gene2symbol1[,2]!=""),1:2]
         gene2symbol<-data.frame(gene=paste(org,as.character(gene2symbol2[,1]),sep=":"),
             symbol=paste("symbol",sapply(gene2symbol2[,2],function(x) unlist(strsplit(x,", "))[1],USE.NAMES=FALSE),sep=":"))
	     assign("keggGene2gene",gene2symbol,envir=k2ri)
	  }
	  
      if(verbose==TRUE)
        print("download relations between gene and enzyme.")
	  file2<-paste(path,"/",org,"/",org,"_enzyme.list",sep="")
      gene2ec<-read.table(file2,header=FALSE,sep = "\t", quote="\"")
      assign("gene2ec",gene2ec,envir=k2ri)
      
      if(verbose==TRUE)
        print("download relations between gene and KO.")
	  file3<-paste(path,"/",org,"/",org,"_ko.list",sep="")
	  gene2ko<-read.table(file3,header=FALSE,sep = "\t", quote="\"")
      assign("gene2ko",gene2ko,envir=k2ri)
	  
	  if(verbose==TRUE)
      print("download relations between KEGG gene and pathway...............")
	  file4<-paste(path,"/",org,"/",org,"_pathway.list",sep="")
	  gene2path<-read.table(file4,header=FALSE,sep = "\t", quote="\"")	  
      assign("gene2path",gene2path,envir=k2ri)
	  
      assign("orgAndIdType",c(org,idType),envir=k2ri)  
}
#updateCompound()
updateCompound<-function(path="ftp://ftp.genome.jp/pub/kegg/ligand",verbose=TRUE){
      if(!exists("k2ri")) initializeK2ri()
      if(verbose==TRUE){
        print("Note that the programming may be time consumming!")  
      }
  
      if(verbose==TRUE)
        print("download compound data.")
	  file1<-paste(path,"/ligand_update.lst",sep="")
	  compound<-as.character(read.table(file1,header=FALSE,sep = "\t")[,1])
	  compound<-compound[substring(compound,0,1)=="C"]
      assign("compound",compound,envir=k2ri)
}

#
updateTaxonomy<-function(path="ftp://ftp.genome.jp/pub/kegg/genes",verbose=TRUE){
      if(verbose==TRUE){
        print("Note that the programming may be time consumming!")  
      }
     taxonomy<-read.table(paste(path,"/","taxonomy",sep=""),header=FALSE,sep = "\t", quote="\"",fill=TRUE,stringsAsFactors=FALSE)
	 assign("taxonomy",taxonomy,envir=k2ri)
}

#updateOrgAndIdType("hsa","ncbi-geneid")
################################################################
#updatPathway
updatePathway<-function(path="ftp://ftp.genome.jp/pub/kegg/xml/",verbose=TRUE){
	  #update metabolicEC
	  if(verbose==TRUE)
	     print("update metabolicEC")
	  detail_path<-paste(path,"kgml/metabolic/ec/",sep="")
	  pathwayList<-get("metabolicEC",envir=k2ri)
	  files_name<-names(pathwayList)
	  metabolicEC<-getPathway(detail_path,files_name,verbose=verbose)
	  assign("metabolicEC",metabolicEC,envir=k2ri)
	  
	  #update metabolicKO
	  if(verbose==TRUE)
	     print("update metabolicKO")
	  detail_path<-paste(path,"kgml/metabolic/ko/",sep="")
	  pathwayList<-get("metabolicKO",envir=k2ri)
	  files_name<-names(pathwayList)
	  metabolicKO<-getPathway(detail_path,files_name,verbose=verbose)
	  assign("metabolicKO",metabolicKO,envir=k2ri)
	  
	  #update nonMetabolicKO
	  if(verbose==TRUE)
	     print("update nonMetabolicKO")
	  detail_path<-paste(path,"kgml/non-metabolic/ko/",sep="")
	  pathwayList<-get("nonMetabolicKO",envir=k2ri)
	  files_name<-names(pathwayList)
	  nonMetabolicKO<-getPathway(detail_path,files_name,verbose=verbose)
	  assign("nonMetabolicKO",nonMetabolicKO,envir=k2ri)
	  
}
#updatPathway()


#########################################################################
#importPathway
importPathway<-function(path,verbose=TRUE){
	  #import kgml/metabolic/ec
	  if(verbose==TRUE)
	     print("import kgml/metabolic/ec")
	  detail_path<-paste(path,"/","kgml/metabolic/ec/",sep="")
	  metabolicEC<-getPathway(detail_path,filelist=list.files(detail_path),verbose=verbose)
	  assign("metabolicEC",metabolicEC,envir=k2ri)
	  
	  #import kgml/metabolic/ko
	  if(verbose==TRUE)
	     print("import kgml/metabolic/ko")
	  detail_path<-paste(path,"/","kgml/metabolic/ko/",sep="")
	  metabolicKO<-getPathway(detail_path,filelist=list.files(detail_path),verbose=verbose)
	  assign("metabolicKO",metabolicKO,envir=k2ri)
	  
	  #import kgml/non-metabolic/ko
	  if(verbose==TRUE)
	     print("import kgml/non-metabolic/ko")
	  detail_path<-paste(path,"/","kgml/non-metabolic/ko/",sep="")
	  nonMetabolicKO<-getPathway(detail_path,filelist=list.files(detail_path),verbose=verbose)
	  assign("nonMetabolicKO",nonMetabolicKO,envir=k2ri)
	  
}

###########################################################################
##save k2ri environment
saveK2ri<-function(file="k2ri.rda"){
      if(!exists("k2ri")) initializeK2ri()
      save(k2ri,file=file,envir=.GlobalEnv)
}
###########################################################################
##load k2ri environment
loadK2ri<-function(file="k2ri.rda"){
      if(!exists("k2ri")) initializeK2ri()
      load(file=file,envir=.GlobalEnv)
}
#updateOrgAndIdType("hsa","ncbi-geneid")
#savek2ri()

###########################################################################
getOrgAndIdType<-function(){
      if(!exists("k2ri")) initializeK2ri()
      orgAndIdType<-get("orgAndIdType",envir=k2ri)
      return(orgAndIdType)
}


