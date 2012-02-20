\name{k2riEnvironment}
\alias{saveK2ri}
\alias{loadK2ri}
\title{Load and save the environment variable k2ri of the system}
\description{
  Load and save the environment variable k2ri of the system.
}
\usage{
    saveK2ri(file="k2ri.rda")
    loadK2ri(file="k2ri.rda")
}
\arguments{
  \item{file}{ A character string.}
}
\details{
The two functions are used to load or save the environment variable \code{\link{k2ri}} of the system. If one has changed the environment variable with some functions (e.g., \code{\link{updateOrgAndIdType}}) and want to use the setting in the future, then the functions can address the problem. For example, one uses the function \code{saveK2ri} to save the environment variable. When one needs to use the setting next time, one can use the function \code{\link{loadK2ri}} to load the last environment variable.

The functions implement the localization of the system. The functions are very important for the user to frequently annotate genes with the different genomes and the type of identifiers. It can increase largely the running speed because one can update the data one time only, and then repeatedly use it.
 
Note that if and only if the function \code{loadK2ri} runs after loading the \code{iSubpathwayMiner} package, it will be in effect.



In addition, the environment variables of organisms with well annotated genomes are provided by the system and users can use the function \code{\link{data}} to load them.

}
\author{ Chunquan Li \email{lcqbio@yahoo.com.cn}}
\seealso{\code{\link{updateOrgAndIdType}}, \code{\link{updatePathway}}, \code{\link{importPathway}}}

\examples{
\dontrun{
##for example 
#save k2ri environment about yeast. 
#Note that the data is saved to the working directory 

path<-paste(system.file(package="iSubpathwayMiner"),"/localdata",sep="")
updateOrgAndIdType("sce","sgd-sce",path)
saveK2ri("sce_ncbi-geneid.rda")

#shut down the R system
#start up the R system
#library(iSubpathwayMiner)

#when the user annotate yeast genes. Note that you need to 
#working directory to the directory of the data file.
loadK2ri("sce_ncbi-geneid.rda") 

}
}
\keyword{ file }

