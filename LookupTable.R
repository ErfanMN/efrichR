##### Preprocess for Sigma and Mean #####


MyRange <- 100
LookupTable <- array(NA , c(nrow(dat),MyRange,2))


for(i in 1:nrow(dat))
{
  for(k in 1:MyRange)
  {
    GeneList <- strsplit2(MainGenList[i , 3] , ",")
    GeneListSize <- length(GeneList)
    Pvalues <- vector()
    for(j in 1:200)
    {
      RandomGeneSet <- sample(AllGenes[,1] , 10*k , replace = FALSE)
      RandomGeneSetSize <- length(RandomGeneSet)
      LCommonGenes <- length(intersect(as.vector(trimws(GeneList)) , as.vector(trimws(RandomGeneSet))))
      CurrentFisher <- fisher.test(matrix(c(length(unique(AllGenes[,1])) - (RandomGeneSetSize+GeneListSize  - LCommonGenes),RandomGeneSetSize - LCommonGenes , GeneListSize - LCommonGenes ,LCommonGenes  ) , nrow =2) , alternative="greater")
      Pvalues[j] <- CurrentFisher$p.value
    }
    
    LookupTable[i,k,1] <- mean(Pvalues)
    LookupTable[i,k,2] <- sd(Pvalues)
  }
  print(i)
}


saveRDS(LookupTable, file="Data/MyLookupTable.rda")