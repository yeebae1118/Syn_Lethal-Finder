#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(stringi)
library(rsconnect)



## load the data
DEP_Map_effect=read.csv("Achilles_gene_effect.csv")
CCLE_RNA_data=read.delim("CCLE_expression.csv")
sample_info=read.csv("sample_info.csv")

##process the data
#the number of cell line with cripsr dependency
dim(DEP_Map_effect)
CellLine_DEPeffect=as.data.frame(DEP_Map_effect)
colnames(CellLine_DEPeffect)[1]="DepMap_ID"
CellLine_DEPeffect %>% 
    mutate(DepMap_ID=as.character(DepMap_ID))->CellLine_DEPeffect
dim(CellLine_DEPeffect)
genelist_Dep=sub("\\..*","",colnames(CellLine_DEPeffect)[2:18334])
colnames(CellLine_DEPeffect)[2:18334]=genelist_Dep

#get the cell name from sample_info(cell line with crispr dependency value)
sample_info %>% 
    dplyr::mutate(DepMap_ID=as.character(DepMap_ID)) %>% 
    dplyr::inner_join(CellLine_DEPeffect[,1:2],by="DepMap_ID") ->Crispr_dep_cell_info
dim(Crispr_dep_cell_info)
View(Crispr_dep_cell_info)

#choose the RNAseq of cell with cripsr dependency score
gene_id=colnames(CCLE_RNA_data)[2:19145]
gene_string=str_split(gene_id,"\\..")
gene_list=lapply(gene_string, function(i){
    gene_id=i[1]
})
gene_list=unlist(gene_list,use.names = F)
colnames(CCLE_RNA_data)[2:19145]=gene_list
dim(CCLE_RNA_data)

##there are some duplicates gene in the list  
##then I calculate the mean of duplicated column
#non duplicated gene

duplcaited_gene=unique(gene_list[duplicated(gene_list)])
non_duplicated=unique(gene_list)[-which(unique(gene_list)%in%duplcaited_gene)]
length(non_duplicated)
for (i in 1:length(duplcaited_gene)) {
    gene=duplcaited_gene[i]
    data=CCLE_RNA_data[,colnames(CCLE_RNA_data)==gene]
    data %>% 
        as.data.frame() %>% 
        mutate(mean=rowMeans(data[,1:dim(data)[2]])) %>% 
        dplyr::select(mean)->gene_mean
    colnames(gene_mean)=duplcaited_gene[i]
    if(!exists("duplicated_Genetable_1")){
        duplicated_Genetable_1=gene_mean
    }
    else{
        duplicated_Genetable_1=cbind(duplicated_Genetable_1,gene_mean)
    }
}
duplicated_Genetable_1=duplicated_Genetable_1[,2:dim(duplicated_Genetable_1)[2]]

#bind the table
CCLE_RNA_filtered=cbind(CCLE_RNA_data[,c(1,which(colnames(CCLE_RNA_data)%in%non_duplicated))],duplicated_Genetable_1)

#get the RNAseq after filtration
CCLE_RNA_filtered %>% 
    dplyr::rename(DepMap_ID=X) %>% 
    mutate(DepMap_ID=as.character(DepMap_ID)) %>% 
    inner_join(Crispr_dep_cell_info[,1:2],by="DepMap_ID")->RNAseq_of_CRISPRcell
RNAseq_of_CRISPRcell_1=RNAseq_of_CRISPRcell[,2:19038]
dim(RNAseq_of_CRISPRcell_1)


## gene dependency correlation analysis
gene_dependency_test=function(gene){
    interested_gene_expression=RNAseq_of_CRISPRcell_1[,c(19037,which(colnames(RNAseq_of_CRISPRcell_1)==gene))]
    mean_expression=mean(interested_gene_expression[,2])
    interested_gene_expression %>%
        mutate(condition=ifelse(.[,2]>mean_expression,"High","Low"))->expression_table2
    
    #correlation analysis
    DEP_sample=Crispr_dep_cell_info[,c(1,2)]
    
    interested_DEP_gene=lapply(1:18333, function(x){
        Dep_score=CellLine_DEPeffect[,c(1,x+1)]
        colnames(Dep_score)[2]="Dep_score"
        Dep_score %>% 
            dplyr::filter(!is.na(Dep_score))->Dep_score
        Dep_score %>% 
            left_join(DEP_sample,by="DepMap_ID")->Data
        gene=genelist_Dep[x]
        Data %>% 
            left_join(expression_table2, by="stripped_cell_line_name")->test.data
        correlation=cor.test(test.data[,2],test.data[,4])
        p.value=correlation$p.value
        t.score=correlation$statistic
        cor.score=correlation$estimate
        mean=mean(Dep_score$Dep_score)
        table.cor=data_frame(gene,p.value,t.score,cor.score,mean)
        table.cor=as.data.frame(table.cor)
    })
    
    correlation_table=do.call(rbind,interested_DEP_gene)
    correlation_table %>% 
        arrange(cor.score) %>% 
        filter(mean< -0.5) ->p
}


##  make shiny query website struture 
ui=fluidPage(
    # App title ----
    titlePanel("Syn-Lethal Finder"),
    # Input content 
    sidebarLayout(
        sidebarPanel(
            helpText("Input your favorite gene
               and find its synthetic lethality partner"),
            
            textInput("text", h3("Gene"), 
                      value = "Enter Gene name..."),
            
            sliderInput("range", 
                        label = "Cutoff",
                        min = 1, max = 250, value = c(1, 250)),
            radioButtons("radio", h3("Gene Expression Level"),
                         choices = list("Low Expression" = 1, "High Expression" =2),selected = 1),
            actionButton("action","Submit")
            
        ),
        
        mainPanel(
            tableOutput("result")
        )
    ))

server= function(input, output) {
    observeEvent(input$action,{
        Result <-  reactive({
            data=gene_dependency_test(input$text)
            if(input$radio==2){
                data
            }else{
                data %>% 
                    arrange(desc(cor.score))
            }
        })
        
        output$result <- renderTable({
            Result()[input$range[1]:input$range[2],]
        })
    })
}

#load application
shinyApp(ui =ui, server = server)
