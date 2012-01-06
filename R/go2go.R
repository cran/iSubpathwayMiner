
####################################################################
##
getBackground<-function(type="gene"){
      if(!exists("k2ri")) initializeK2ri()
	  if(type=="gene"){
	  keggGene2gene<-get("keggGene2gene",envir=k2ri) 
	  background<-unique(as.character(keggGene2gene[,2]))
	  newBackground<-sapply(strsplit(background,":"),function(x) x[2])
	  }
	  else if(type=="compound"){
	  newBackground<-unique(get("compound",envir=k2ri))
	  }
	  else if(type=="gene_compound"){
	  keggGene2gene<-get("keggGene2gene",envir=k2ri) 
	  background<-unique(as.character(keggGene2gene[,2]))
	  newBackground<-sapply(strsplit(background,":"),function(x) x[2])	  
	  newBackground<-union(newBackground,get("compound",envir=k2ri))
	  }
      return(newBackground)
}
#############################################################
##get kegg gene list from gene list
getKGeneFromGene<-function(geneList){
	  geneList<-as.character(geneList)
      if(!exists("k2ri")) initializeK2ri()
	  keggGene2gene<-get("keggGene2gene",envir=k2ri)
keggGeneList<-unique(as.character(keggGene2gene[as.character(keggGene2gene[,2]) %in% paste(getOrgAndIdType()[2],geneList,sep=":"),1]))
      return(keggGeneList)
}
#keggGeneList<-getKeggGeneFromGene(geneList)
#############################################################
##get gene list from kegg gene list
getGeneFromKGene<-function(keggGeneList){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("k2ri")) initializeK2ri()
      keggGene2gene<-get("keggGene2gene",envir=k2ri)
      geneList<-unique(as.character(sapply(strsplit(as.character(keggGene2gene[as.character(keggGene2gene[,1]) %in% keggGeneList,2]),":"),function(x) return (x[2]))))
      return(geneList)
}
#getGeneFromKeggGene(keggGeneList)
#############################################################
##get enzyme list from gene List
getEnzymeFromKGene<-function(keggGeneList,ignoreAmbiguousEnzyme=TRUE){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("k2ri")) initializeK2ri()
	  gene2ec<-get("gene2ec",envir=k2ri)
      enzymeList<-unique(as.character(gene2ec[as.character(gene2ec[,1]) %in% keggGeneList,2]))
	  if(ignoreAmbiguousEnzyme==TRUE){
	     enzymeList<-grep("-",enzymeList,value=TRUE,invert=TRUE)
	  }	  
      return(enzymeList)
}
#enzymeList<-getEnzymeFromKeggGene(keggGeneList)
#############################################################
##get gene list from enzyme List
getKGeneFromEnzyme<-function(enzymeList,ignoreAmbiguousEnzyme=TRUE){
	  enzymeList<-as.character(enzymeList)
	  if(ignoreAmbiguousEnzyme==TRUE){
	     enzymeList<-grep("-",enzymeList,value=TRUE,invert=TRUE)
	  }
      if(!exists("k2ri")) initializeK2ri()
	  gene2ec<-get("gene2ec",envir=k2ri)
      keggGeneList<-unique(as.character(gene2ec[as.character(gene2ec[,2]) %in% enzymeList,1]))
      return(keggGeneList)
}
#getKeggGeneFromEnzyme(enzymeList)
#############################################################
##new! get KO list from gene List
getKOFromKGene<-function(keggGeneList){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("k2ri")) initializeK2ri()
	  gene2ko<-get("gene2ko",envir=k2ri)
      KOList<-unique(as.character(gene2ko[as.character(gene2ko[,1]) %in% keggGeneList,2]))
      return(KOList)
}
#KOList<-getKOFromKeggGene(keggGeneList)
#############################################################
##new! get gene list from KO List
getKGeneFromKO<-function(KOList){
	  KOList<-as.character(KOList)
      if(!exists("k2ri")) initializeK2ri()
	  gene2ko<-get("gene2ko",envir=k2ri)
      keggGeneList<-unique(as.character(gene2ko[as.character(gene2ko[,2]) %in% KOList,1]))
      return(keggGeneList)
}
#getKeggGeneFromKO(KOList)
#############################################################
##new! get pathway from gene List
getPathwayFromKGene<-function(keggGeneList){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("k2ri")) initializeK2ri()
	  gene2path<-get("gene2path",envir=k2ri)
      pathwayList<-unique(as.character(gene2path[as.character(gene2path[,1]) %in% keggGeneList,2]))
	  pathwayList<-paste("path:",substring(pathwayList,9),sep="")
      return(pathwayList)
}
#KOList<-getKOFromKeggGene(keggGeneList)
#############################################################
##new! get gene list from pathway List
getKGeneFromPathway<-function(pathwayList){
	  pathwayList<-as.character(pathwayList)
      if(!exists("k2ri")) initializeK2ri()
	  gene2path<-get("gene2path",envir=k2ri)  
      keggGeneList<-unique(as.character(gene2path[as.character(gene2path[,2]) %in% paste("path:",getOrgAndIdType()[1],substring(pathwayList,6),sep=""),1]))
      return(keggGeneList)
}



#############################################################new!
##get kegg gene list from gene list
getKGeneFromSymbol<-function(symbolList){
	  symbolList<-as.character(symbolList)
      if(!exists("k2ri")) initializeK2ri()
	  gene2symbol<-get("gene2symbol",envir=k2ri)
      keggGeneList<-unique(as.character(gene2symbol[as.character(gene2symbol[,2]) %in% paste("symbol",symbolList,sep=":"),1]))
      return(keggGeneList)
}
#keggGeneList<-getKeggGeneFromGene(geneList)
#############################################################new!
##get gene list from kegg gene list
getSymbolFromKGene<-function(keggGeneList){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("k2ri")) initializeK2ri()
      gene2symbol<-get("gene2symbol",envir=k2ri)
      symbolList<-unique(as.character(sapply(strsplit(as.character(gene2symbol[as.character(gene2symbol[,1]) %in% keggGeneList,2]),":"),function(x) return (x[2]))))
      return(symbolList)
}
#####################################################################new!
#get org ID (e.g.,"sce") from org name (e.g.,"s.cerevisiae")
getOrgIdFromOrgName<-function(orgName){
	  orgName<-as.character(orgName)
      if(!exists("k2ri")) initializeK2ri()
	  taxonomy<-get("taxonomy",envir=k2ri)  
      orgId<-unique(as.character(taxonomy[as.character(taxonomy[,3]) %in% orgName,2]))
      return(orgId)
}

#####################################################################new!
#get org ID (e.g.,"sce") from org name (e.g.,"s.cerevisiae")
getOrgNameFromOrgId<-function(orgId){
	  orgId<-as.character(orgId)
	  #orgId<-tolower(orgId)
      if(!exists("k2ri")) initializeK2ri()
	  taxonomy<-get("taxonomy",envir=k2ri)  
      orgName<-unique(as.character(taxonomy[as.character(taxonomy[,2]) %in% orgId,3]))
      return(orgName)
}



#############################################################
##get kegg gene list from gene list
getEnzymeFromGene<-function(geneList,ignoreAmbiguousEnzyme=TRUE){
      return(getEnzymeFromKGene(getKGeneFromGene(geneList),ignoreAmbiguousEnzyme=ignoreAmbiguousEnzyme))
}
getGeneFromEnzyme<-function(enzymeList,ignoreAmbiguousEnzyme=TRUE){
      return(getGeneFromKGene(getKGeneFromEnzyme(enzymeList,ignoreAmbiguousEnzyme=ignoreAmbiguousEnzyme)))
}
getKOFromGene<-function(geneList){
      return(getKOFromKGene(getKGeneFromGene(geneList)))
}
getGeneFromKO<-function(KOList){
      return(getGeneFromKGene(getKGeneFromKO(KOList)))
}

getPathwayFromGene<-function(geneList){
      return(getPathwayFromKGene(getKGeneFromGene(geneList)))
}
getGeneFromPathway<-function(pathwayList){
      return(getGeneFromKGene(getKGeneFromPathway(pathwayList)))
}
#new!
getSymbolFromGene<-function(geneList){
      return(getSymbolFromKGene(getKGeneFromGene(geneList)))
}
getGeneFromSymbol<-function(symbolList){
      return(getGeneFromKGene(getKGeneFromSymbol(symbolList)))
}

getEnzymeFromSymbol<-function(symbolList){
      return(getEnzymeFromKGene(getKGeneFromSymbol(symbolList)))
}
getSymbolFromEnzyme<-function(enzymeList){
      return(getSymbolFromKGene(getKGeneFromEnzyme(enzymeList)))
}
getKOFromSymbol<-function(symbolList){
      return(getKOFromKGene(getKGeneFromSymbol(symbolList)))
}
getSymbolFromKO<-function(KOList){
      return(getSymbolFromKGene(getKGeneFromKO(KOList)))
}