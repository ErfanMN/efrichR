library(ggplot2)
library(limma)


################# Structure of main Result & it's Components ######################
result <- data.frame(ID = 1:nrow(MainGenList) , Name = MainGenList[,1], Subscribe = 0 , P_Value = NA , Adj_P_Value = NA , Z_Score = NA , Combined_Score = NA , Common_Genes = NA)
SortedResult <- result
TopSortedResult <- SortedResult[1:5,]
UserGeneList <- data.frame()
options(scipen=3)

Mapper <- read.csv("Data/pathwaysUrl.txt" , header = T , sep = "\t" , colClasses=c('character', 'character'))



server <- function(input, output) {
 
  
ntext <- eventReactive(input$goButton, {

    
    UserGeneList <<- strsplit2(input$TextArea,"\n")
    a <- ncol(UserGeneList)
    for(i in 1:nrow(dat))
    {
      GeneList <- strsplit2(MainGenList[i , 3] , ",")
      b <- length(GeneList)
      CommonGenes <- intersect(as.vector(trimws(GeneList)) , as.vector(trimws(UserGeneList[1,])))
      result[i,3] <<- length(CommonGenes)
      result[i,8] <<- toString(CommonGenes)
      fisherResult <- fisher.test(matrix(c(length(unique(AllGenes[,1])) - (a+b  - result[i,3] ),a - result[i,3] , b - result[i,3] , result[i,3]  ) , nrow =2 ) , alternative = "greater")
      result[i,4] <<- fisherResult$p.value
    }

    padjusted <- p.adjust(result[,4] , method = "BH"  )
    result[,5] <<- padjusted
    result[,6] <<- (result[,4] - LookupTable[,a/10,1]) / LookupTable[,a/10,2]
    result[,7] <<- result[,6]*log(result[,5])

    
    result <<- result[is.finite(rowSums(result[,3:7])),]
    SortedResult <<- result[order(result[,7] , decreasing = T),]
    write.table(SortedResult[,2:8] , "Result/result.txt" , sep = "\t")
    TopSortedResult <<- SortedResult[1:5,]
    src <- paste(c("https://www.genome.jp/kegg/pathway/hsa/hsa", Mapper[TopSortedResult[1,1],2] , ".png"), collapse = "")
    output$myImage = renderText(
      {c('<img src="',src,'" width="100%">')}
    )
    TopSortedResult[5:1,]

  })
  
  output$plot <- renderPlot({
    ggplot(data=ntext(), aes(x=c(1,2,3,4,5) , y=ntext()$Combined_Score , fill = c(1,2,3,4,5) )) +
    geom_bar(stat="identity") +  
      scale_fill_gradient(low = "#ff1a66", high = "#80002a")+
      coord_flip() + 
      geom_text(aes(label=MainGenList[ntext()$ID,1] ),hjust = 1 , color="black",
                                                                      position = position_dodge(0), size=3.5) + 
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) 
  })
  
  output$mytable = DT::renderDataTable(
    data.frame(Pathway = ntext()[5:1,2] , format(ntext()[5:1,4:7], format = "e", digits = 4) , CommonGenes = ntext()[5:1,8])
                                       
                                       , rownames = FALSE )
  

  
}