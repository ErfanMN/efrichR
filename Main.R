library(limma)
library(shiny)
library(pracma)
library(ggplot2)
library(rsconnect)
setwd(".//")

################# Load All Genes that exist ######################

AllGenes <- as.vector(read.csv("Data/AllGenes.txt"))

################# Read & Load All Genes as MainGeneList###########
file <- "GenList.gmt"
dat <- read.delim(file, sep = " " , header = F , stringsAsFactors= F )
MainGenList <- data.frame()

for(i in 1:nrow(dat)){
  datTable <- as.character(strsplit2(dat[i,1] , "\t"))
  MainGenList[i , 1] <- gsub("^.*?_","",datTable[1])
  MainGenList[i , 2] <- datTable[2]
  MainGenList[i , 3] <- toString(datTable[3:length(datTable)])
}

################## if you want to create LookupTable you should uncomment this part #############
LookupTable <- readRDS("Data/MyLookupTable.rda")
#source("LookupTable.R")


################# Structure of main Result & it's Components ######################
result <- data.frame(ID = 1:nrow(MainGenList) , Name = MainGenList[,1], Subscribe = 0 , P_Value = NA , Adj_P_Value = NA , Z_Score = NA , Combined_Score = NA , Common_Genes = NA)
SortedResult <- result
TopSortedResult <- SortedResult[1:5,]
UserGeneList <- data.frame()
options(scipen=3)

################ Get User gene List ###############################################
UserGeneList <- read.delim("Data/list.txt" , sep = "\n")


################ Fisher Test P Value ->>>   result[,4] ######################
a <- nrow(UserGeneList)
for(i in 1:nrow(dat))
{
  GeneList <- strsplit2(MainGenList[i , 3] , ",")
  b <- length(GeneList)
  CommonGenes <- intersect(as.vector(trimws(GeneList)) , as.vector(trimws(UserGeneList[1,])))
  result[i,3] <- length(CommonGenes)
  #################### Save CommonGenes Between Each GeneList And UserGeneList ###############
  result[i,8] <- toString(CommonGenes)
  fisherResult <- fisher.test(matrix(c(length(unique(AllGenes[,1])) - (a+b  - result[i,3] ),a - result[i,3] , b - result[i,3] , result[i,3]  ) , nrow =2 ) , alternative = "greater")
  result[i,4] <- fisherResult$p.value
}

################ Adjusted P Value ################################
padjusted <- p.adjust(result[,4] , method = "BH"  )
result[,5] <- padjusted


################ Finding Z-Score Using LookupTable Created Before #################
result[,6] <- (result[,4] - LookupTable[,a/10,1]) / LookupTable[,a/10,2]



############### Calculating Combined Score #######################################
result[,7] <- result[,6]*log(result[,5])


############### Removing infinte Values in z score and Combined score ############
result <- result[is.finite(rowSums(result[,3:7])),]
############### Sorting The Result in order of Combined score
SortedResult <- result[order(result[,7] , decreasing = T),]
############### Create and Save the Result in a text file #########################
write.table(SortedResult[,2:8] , "Result/result.txt" , sep = " ")


################## User interface that user can Enter the data ###################
source("ui.R")
################## Server Side Code that Contain Main Processes ##################
source("server.R")
################## Run a Web based Application ###################################
shinyApp(ui, server)


