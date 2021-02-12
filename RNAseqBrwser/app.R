library(dplyr)
library(tibble)
library(tidyr)
library(limma)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(cluster)
library(shiny)
library(shinythemes)
library(plotly)
library(ggplot2)
library(gridExtra)
library(manhattanly)
library(heatmaply)
library(htmlwidgets)
library(ggpubr)
library("ggthemes")
library(Rtsne)
library(shinyBS)


ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("RNA-seq analysis & Visualization"),
  sidebarPanel(
    selectInput("gct",label = h3("Choose a experiment to browse"),c("KP_Longterm_Exh","Scottbrowne_Exh","gse65850_Tox_ko","GSE58596","GSE83482","GSE89307_schietinger","KP_JG_Flu_PDL1L2KO","EMTAB_6214","GSE84820","GSE83978","AHUA_C1IpiC1Pembro","GSE99531","GSE94581_Vignali_Tregs","AHUA_PBMC_4cycles_Pembro_RNAseq"),selected="GSE94581_Vignali_Tregs")
    , bsTooltip("gct", "Browse through RNAseq datasets by name and choose your dataset of interest",
                "right", options = list(container = "body")),
    uiOutput("groups_input"),
    bsTooltip("groups_input", "Choose groups to be displayed ",
              "right", options = list(container = "body")),
    textInput("gene_input", label=h4("Enter a gene name to check expression"), value = "Lag3"),
    textAreaInput("genes_input", label=h5("Paste your genes of interest to highlight"), "Pdcd1 Tox Batf Ctla4 Lag3"),
    br(),
    uiOutput("c"),
    p(strong("Choose parameters for differential gene expression analysis for groups choosen")),
    selectInput("adjusttype",label=h4("Choose p-value adjustment"),c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none"),selected="none"),
    numericInput(
      "p_value",
      label = h4("Enter p-value cutoff"),
      value = 0.05,max = 1
    ),
    p(strong("Term analysis options")),
    selectInput("species", label = h4("Select Organism for Term analysis"), 
                choices = list("Human" = "Hs", "Mouse" = "Mm"), 
                selected = "Mm"),
    actionButton("runButton", "Run!"),
    br(),
    textOutput("currentTime")
  ),
  mainPanel(
    tabsetPanel(
    type = "tabs",
    tabPanel("Gene expression",plotOutput("gexp",width = "500px",height="500px")
    ),
    tabPanel("Custom Heatmap",plotOutput("custom_heatmap",height=600)),
      tabPanel("PCA plot", plotOutput("pca",width = "800px",height = "400px")),
      tabPanel("t-SNE plot", plotOutput("tsne",width = "600px",height = "400px"))
    ),
  tabsetPanel(
    type = "tabs",
    tabPanel("Volcano plot", plotlyOutput("vplot",width = "500px",height = "500px")),
    tabPanel("MA plot", plotlyOutput("mplot",width = "500px",height = "500px")),
    tabPanel("DE Genes list & Statistics", dataTableOutput("table"))
  ),
  tabsetPanel(
    type = "tabs",
     tabPanel("Term enrichment summary",splitLayout(
      plotOutput("plot1"),
      plotOutput("plot2")
    ),width=1000,height=1000),
     tabPanel("Term enrichment list for UP genes", dataTableOutput("GO_table_UP")),
     tabPanel("Term enrichment list for DOWN genes", dataTableOutput("GO_table_DOWN")),
    tabPanel("Heatmap of Top 100 significant genes",plotOutput("p_heatmap",width = "600px",height = "1200px")),
    tabPanel("Heatmap of Top 100 significant genes in all samples",plotOutput("p_heatmap2",width = "600px",height = "1200px"))
    # tabPanel("Top GSEA hits",plotOutput("gseaplot",width="70%",height="400px"))
  ),
  downloadButton('downloadData', 'Download gene statistics'),
  downloadButton('downloadData_GO', 'Download Term enrichment result')
)
)
server <- function(input,output,session) {
  options(shiny.maxRequestSize = 60 * 1024 ^ 2)
  
  output$c <- renderUI({
    if (is.null(input$gct))
      return(NULL)
    #inFile <- input$ExpDes
    input_data<-input$gct
    arraydatafiles<-read.delim("arraydatafiles.txt",sep="\t",header = TRUE,row.names = 1)
    expdes<-read.delim(as.character(arraydatafiles[input_data,]$expdesfile),stringsAsFactors = FALSE)
    
    
    #Expdes <- read.delim(inFile$datapath,stringsAsFactors = FALSE)
    groups <- unique(expdes[,2])
    checkboxGroupInput("groups", label = h4("Choose any two groups to compare"), groups)
  })
  output$groups_input <- renderUI({
    if (is.null(input$gct))
      return(NULL)
    input_data<-input$gct
    arraydatafiles<-read.delim("arraydatafiles.txt",sep="\t",header = TRUE,row.names = 1)
    Expdes<-read.delim(as.character(arraydatafiles[input_data,]$expdesfile),stringsAsFactors = FALSE)
    
    groups <- unique(Expdes[,2])
    checkboxGroupInput("groups_keep", label = h4("Choose groups to keep to check expression"), groups)
  })
  design.pairs <-function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    design
  }
  
  Plot_counts<-function(origdata,data,gene,groups){
    
    gene_search<-paste("^",gene,"$",sep="")
    
    gene_idx<-which(origdata$geneSymbol == gene)
    
    gene_id<-origdata[gene_idx,]$id
    
    if(length(gene_id)>1){
      p<-list()
      for(i in 1:length(gene_id)){
        exp_df_gene<-data[,c(1,2,grep(gene_id[i],colnames(data),ignore.case = TRUE))]
        
        #geneCounts <- plotCounts(data, as.character(exp_df_gene$Entrez_geneID), intgroup=c("Group"), returnData=TRUE,transform = transform)
        
        #boxplot(count~Group,data=geneCounts,pch=24, boxwex=0.4,col=col,main=gene,
        #    cex.lab=1.5, ylab='log10 Concentration', cex.main=1.5,las=2,font.lab=2,outline=FALSE)
        
        colnames(exp_df_gene)<-c("sample","Group","count")
        
        gene_title<-paste0(gene_id[i],"/",gene)
        p[[i]]<-ggplot(exp_df_gene[exp_df_gene$Group %in% groups ,], aes(x=Group, y=count, color=Group))  +ylab("log2 Normalized Counts") +geom_point(position=position_jitter(width=.1,height=0), size=4) + theme_bw() + ggtitle(gene_title)+theme(axis.text.x = element_text(size=12,angle = 45, hjust = 1,face="bold"),legend.position="none",plot.title = element_text(lineheight=.8,hjust=0.5, face="bold"),axis.text.y = element_text(size=12,face="bold"))
        
      }
      grid.arrange(p[[1]],p[[2]],ncol = 2)
    }
    else if(length(gene_id) == 1){
      exp_df_gene<-data[,c(1,2,grep(gene_id,colnames(data),ignore.case = TRUE))]
      
      #geneCounts <- plotCounts(data, as.character(exp_df_gene$Entrez_geneID), intgroup=c("Group"), returnData=TRUE,transform = transform)
      
      #boxplot(count~Group,data=geneCounts,pch=24, boxwex=0.4,col=col,main=gene,
      #    cex.lab=1.5, ylab='log10 Concentration', cex.main=1.5,las=2,font.lab=2,outline=FALSE)
      
      colnames(exp_df_gene)<-c("sample","Group","count")
      gene_title<-paste0(gene_id,"/",gene)
      ggplot(exp_df_gene[exp_df_gene$Group %in% groups ,], aes(x=Group, y=count, color=Group))  +ylab("log2 Normalized Counts") +geom_point(position=position_jitter(width=.1,height=0), size=4) + theme_bw() + ggtitle(gene_title)+theme(axis.text.x = element_text(size=12,angle = 45, hjust = 1,face="bold"),legend.position="none",plot.title = element_text(lineheight=.8,hjust=0.5, face="bold"),axis.text.y = element_text(size=12,face="bold"))
    }
  }
  
  
  add_random_noise <- function(x) {
    
    r_value_to_add = rnorm(n=1, mean=0, sd=0.4)
    
    #To avoid negative exp values, add the absolute value of the random noise to 0-count data
    if(x == 0) {
      return(x + abs(r_value_to_add))
    } else {
      #Continue to generate random numbers until the sum of X and the random number
      #is >= 0 (this should rarely need to repeat).
      repeat{
        if(x + r_value_to_add >= 0) break
        r_value_to_add = rnorm(n=1, mean=0, sd=0.4)
      }
      return(x + r_value_to_add)
    }
  }
  
  vobj<-reactive({
    if (is.null(input$gct) )
      return(NULL)
    input_data<-input$gct
    arraydatafiles<-read.delim("arraydatafiles.txt",sep="\t",header = TRUE,row.names = 1)
    sampleData<-read.delim(as.character(arraydatafiles[input_data,]$expdesfile),stringsAsFactors = FALSE)
    
    read_count =read.delim(file = as.character(arraydatafiles[input_data,]$gctfile), row.names = NULL) %>% 
      #Strip off "gene:" prefix for gene ID column (if present)
      mutate(id = gsub("gene:", "", id))
    
    l<-ncol(read_count)
    
    group_samples<-colnames(read_count)[2:(l-2)]
    
    sample_indexes =  which(colnames(read_count) %in% unlist(strsplit(group_samples,split = ",",fixed = TRUE)),arr.ind = TRUE)
    
    data_matrix = 
      read_count %>% 
      dplyr::select(id,sample_indexes) %>% 
      column_to_rownames(var = "id") %>% 
      as.matrix()
    
    data.matrix<-data_matrix
     zero_genes = rownames(data.matrix[apply(data.matrix, 1, function(x) (length(unique(x))==1 & x[1] == 0)),])
    # 
    # #Remove NA values from PORT matrix, as they cause the downstream limma-voom analyses to fail
    # #(these arise from filtering out low expressers)
     data.matrix = data.matrix[!apply(data.matrix, 1, anyNA),]
    # 
    # 
    # #Add random noise to each value in the matrix
     data.matrix[] = vapply(data.matrix, add_random_noise, numeric(1))
    # 
    # #Revert genes with no expression across all samples back to 0 expression
     data.matrix[zero_genes,] = 0
    
    lib.sizes = colSums(data.matrix)
    
    group<-factor(sampleData$condition)
    
    design_matrix <- model.matrix(~group)
    
    port.voom_results = voom(data.matrix,design = design_matrix, lib.size = lib.sizes,
                             plot = FALSE, save.plot = FALSE)
    
    port.voom_results
  })
 
  output$gexp<-renderPlot({
    if (is.null(input$gct) )
      return(NULL)
    if (is.null(input$groups_keep) )
      return(NULL)
    voom_obj<-vobj()
    input_data<-input$gct
    mdata<-voom_obj$E
    
    data_in<-data.frame(t(mdata))
    arraydatafiles<-read.delim("arraydatafiles.txt",sep="\t",header = TRUE,row.names = 1)
    sampleData<-read.delim(as.character(arraydatafiles[input_data,]$expdesfile),stringsAsFactors = FALSE)
    
    read_count =read.delim(file = as.character(arraydatafiles[input_data,]$gctfile), row.names = NULL) %>% 
      #Strip off "gene:" prefix for gene ID column (if present)
      mutate(id = gsub("gene:", "", id))
    groups_keep<-input$groups_keep
    mergData<-merge(sampleData,data_in,by.x="sample",by.y=0)
    gene_in<-input$gene_input
    p<-Plot_counts(read_count,mergData,gene_in,groups_keep)
    p
    
  })
  
  
  
  gene_stats <- eventReactive(input$runButton,{
    if (is.null(input$groups))
      return(NULL)
    voom_obj<-vobj()
    input_data<-input$gct
    arraydatafiles<-read.delim("arraydatafiles.txt",sep="\t",header = TRUE,row.names = 1)
    Expdes<-read.delim(as.character(arraydatafiles[input_data,]$expdesfile),stringsAsFactors = FALSE)
   
    group <- as.factor(Expdes[,2])
    
    print(group)
    design_matrix <- model.matrix(~0+group)
    
    lvl<-levels(group)
    
    print(lvl)
    
    fl <- as.factor(lvl)
    
    print(fl)
    
    contrast.matrix <- design.pairs(lvl)
    
    samples <- rownames(contrast.matrix)
    
    print(samples)
    
    c1 <- input$groups
    
    print(c1)
    
    i1 = which(samples == c1[1])
    
    i2 = which(samples == c1[2])
    
    n <- length(levels(fl))
    
    if (i1 > i2) {
      n1 = i2
      n2 = i1
    }else{
      n1 = i1
      n2 = i2
    }
    
    #Show only choosen pair of comparison.
    c0 = 0
    for (i in 1:n1)
    {
      c0 = c0 + n - i
    }
    cn = c0 - n + n2
    
    comparison_names<-colnames(contrast.matrix)
    
    #Use limma to perform DE analysis
    
    fit<-lmFit(voom_obj,design_matrix)
    
    #contrast.matrix<-makeContrasts(cd4_foxp3_pos-cd4_foxp3_neg,cd4_foxp3_pos-cd8,cd4_foxp3_neg-cd8,levels = design_matrix)
    
    fit2<-contrasts.fit(fit, contrast.matrix)
    
    port.deg_results = eBayes(fit2)
    
    adjust_method<-input$adjusttype
    
    deg_results<-topTable(port.deg_results,coef = cn,number =Inf, adjust.method=adjust_method )
    
    deg_results$P<-deg_results$P.Value
    
    deg_results$EFFECTSIZE<-deg_results$logFC
    
    read_count =read.delim(file = as.character(arraydatafiles[input_data,]$gctfile), row.names = NULL) %>% 
      #Strip off "gene:" prefix for gene ID column (if present)
      mutate(id = gsub("gene:", "", id))
    
    mdata<-merge(deg_results,read_count,by.x=0,by.y="id")
    
   # mdata$gene<-paste0(mdata$Row.names,"_",mdata$geneSymbol)
    
    mdata$gene<-mdata$geneSymbol
    
    mdata$comparison<-comparison_names[cn]
    
    mdata
    
  })
  output$custom_heatmap <- renderPlot({ 
    if (is.null(input$gct) )
      return(NULL)
    if (is.null(input$groups_keep) )
      return(NULL)
    inputdata<-input$gct
    voom_obj<-vobj()
    # 
    data_in<-data.frame(voom_obj$E)
    arraydatafiles<-read.delim("arraydatafiles.txt",sep="\t",header = TRUE,row.names = 1)
    read_count =read.delim(file = as.character(arraydatafiles[inputdata,]$gctfile), row.names = NULL) %>% 
      #Strip off "gene:" prefix for gene ID column (if present)
      mutate(id = gsub("gene:", "", id))
    
    # 
    mergData<-merge(data_in,read_count,by.x=0,by.y="id")
    # 
    Expdes<-read.delim(as.character(arraydatafiles[inputdata,]$expdesfile),stringsAsFactors = FALSE)
    
    samples_interest<-Expdes$sample
    samples_interest_fin<-paste0(samples_interest,".x")
    genes_interest<-strsplit(as.character(input$genes_input)," ")
    print(unlist(genes_interest))
    expr_data_in<-mergData[mergData$geneSymbol %in% unlist(genes_interest), ]

    rownames(expr_data_in)<-expr_data_in$geneSymbol
    input_data <-expr_data_in [,samples_interest_fin]
    rownames(Expdes)<-Expdes$sample
    #Expdes$sample<-NULL
    colnames(input_data)<-gsub(".x","",colnames(input_data))
    groups_keep<-input$groups_keep
    expdes_custom<-Expdes[Expdes[,2] %in% groups_keep,]
    idx<-intersect(rownames(expdes_custom),samples_interest)
    expdes_custom$sample<-NULL
    pheatmap(input_data[,as.character(idx)],cellwidth = 8,cellheight = 10,cluster_cols  = TRUE,scale = "row",annotation_col = expdes_custom)
  })
  
  output$vplot<-renderPlotly({
    deg_result<-gene_stats()
    
    vobj<-volcanor(deg_result,snp="gene")
    
    #volcanoly(vobj,genomewideline= 1.301) %>%
     # saveWidget(file="volcanoly_test.html",selfcontained = TRUE)
    comparison_header<-unique(deg_result$comparison)
    print(comparison_header)

    genes_interest<-strsplit(as.character(input$genes_input)," ")
    print(genes_interest)
    volcanoly(vobj,genomewideline= 1.301,title=comparison_header,highlight = genes_interest[[1]],col="grey",highlight_color = "green",effect_size_line_color="red",genomewideline_color="red")
      
  })
  
  
  output$mplot<-renderPlotly({
    deg_result<-gene_stats()
    
    x = list(title = "log2 fold change",showgrid=F)
    y = list(title = "Mean expression")
    m = list(
      l = 50,
      r = 50,
      b = 100,
      t = 100,
      pad = 4
    )
    f <- list(
      family = "Arial",
      size = 16,
      color = "#7f7f7f"
    )
    comparison_hdr<-unique(deg_result$comparison)
    p <-plot_ly(deg_result,x = ~logFC,y = ~AveExpr,type = "scatter",mode = "markers",marker = list(size =7),text = ~gene,opacity=0.6)
    p %>% layout(title=comparison_hdr,margin=m,autosize = F,font=f, width = 500, height = 500, xaxis = x, yaxis = y,hovermode="closest")
    
  })
  
  output$table <- renderDataTable({
    if (is.null(input$groups))
      return(NULL)
    gene_Table<- gene_stats()
    # p_value <- -log10(input$p_value)
    # fc_cutoff<-log2(input$FC)
    # gene_Table<-gene_Table[gene_Table$lp> p_value,]
    # gene_Table<-gene_Table[abs(gene_Table$logFC)>fc_cutoff, ]
    gene_Table<-gene_Table[order(gene_Table$adj.P.Val),]
    gene_Table[,c("Row.names","geneSymbol","logFC","AveExpr","P.Value","adj.P.Val")]
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$groups[1],input$groups[2],"_genes_stats", '.csv', sep = '')
    },
    content = function(file) {
      df<-gene_stats()
      df<-df[order(df$adj.P.Val),]
      write.csv(df[,c("Row.names","geneSymbol","logFC","AveExpr","P.Value","adj.P.Val")], file)
    }
  )
  output$p_heatmap <- renderPlot({
    if (is.null(input$groups))
      return(NULL)
     dtable<-gene_stats()
     voom_obj<-vobj()
    # 
     mdata<-voom_obj$E
    # 
     data_in<-data.frame((mdata))
    # 
     mergData<-merge(data_in,dtable,by.x=0,by.y="Row.names")
    # 
     grps <- input$groups
     input_data<-input$gct
     arraydatafiles<-read.delim("arraydatafiles.txt",sep="\t",header = TRUE,row.names = 1)
     Expdes<-read.delim(as.character(arraydatafiles[input_data,]$expdesfile),stringsAsFactors = FALSE)
     i1 = which(Expdes[,2] == grps[1])
     i2 = which(Expdes[,2] == grps[2])
     
     samples_interest<-Expdes[c(i1,i2),]$sample
    # 
     expr_de_data<-mergData[mergData$P.Value < 0.05,]
    # 
     expr_de_data<-expr_de_data[order(expr_de_data$P.Value),]
     
     expr_de_data_sub<-expr_de_data[expr_de_data$AveExpr >2,]
     
     expr_de_data_sub<-expr_de_data_sub[order(expr_de_data_sub$P.Value),]
    # 
     expr_de_data_sub<-expr_de_data_sub[1:100,]
    # 
     rownames(expr_de_data_sub)<-expr_de_data_sub$geneSymbol
    # 
    # expr_de_data_sub[,1]<-NULL
    # 
     samples_interest_fin<-paste0(samples_interest,".x")
     expr_data_in <- expr_de_data_sub[,samples_interest_fin]
     colnames(expr_data_in)<-gsub(".x","",colnames(expr_data_in))
    # 
    col_df<-data.frame(Expdes[c(i1,i2),])
    # 
     rownames(col_df)<-col_df[,1]
    # 
     col_df[,1]<-NULL
    # 
     col_col<-data.frame(t(col_df))
    
    #heatmaply(expr_data_in,grid_color = "black",scale="row",margins = c(200,300,100,100),col_side_colors=col_col) %>%
      #saveWidget(file="heatmaply_test.html",selfcontained = FALSE)
    
   # heatmaply(expr_data_in,grid_color = "black",scale="row",margins = c(200,300,100,100),col_side_colors=col_col)
    fname<-paste0(input$groups[1],input$groups[2],"top100.pdf")
     expr_data_in_fin<-t(scale(t(expr_data_in)))
     rownames(Expdes)<-Expdes$sample
     Expdes$sample<-NULL
     pdf(file = fname,width = 8,height = 20)
     pheatmap(expr_data_in_fin,scale = "row",border_color = NA,cellwidth = 16,cellheight = 10,annotation_col = col_df)
     dev.off()
     
    require(iheatmapr)
     expr_data_in_fin %>%
       as.matrix %>%
       iheatmap %>%
       add_col_labels(font = list(size = 8)) %>%
       add_row_labels(size = 0.3,font = list(size = 8))%>%add_col_clustering(k=2)%>%
       add_row_clustering(k=2)  -> hm
     
    # hm %>% save_iheatmap("myplot.pdf")
     #hm_final<-hm %>% as_plotly()
     #hm_final
     pheatmap(expr_data_in_fin,scale = "row",border_color = NA,cellwidth = 16,cellheight = 10,annotation_col = col_df)
     
     
  })
  
  output$p_heatmap2 <- renderPlot({
    if (is.null(input$groups))
      return(NULL)
    dtable<-gene_stats()
    voom_obj<-vobj()
    # 
    mdata<-voom_obj$E
    # 
    data_in<-data.frame((mdata))
    # 
    mergData<-merge(data_in,dtable,by.x=0,by.y="Row.names")
    # 
    grps <- input$groups
    input_data<-input$gct
    arraydatafiles<-read.delim("arraydatafiles.txt",sep="\t",header = TRUE,row.names = 1)
    Expdes<-read.delim(as.character(arraydatafiles[input_data,]$expdesfile),stringsAsFactors = FALSE)
    i1 = which(Expdes[,2] == grps[1])
    i2 = which(Expdes[,2] == grps[2])
    
    samples_interest<-Expdes$sample
    # 
    expr_de_data<-mergData[mergData$P.Value < 0.05,]
    # 
    expr_de_data<-expr_de_data[order(expr_de_data$P.Value),]
    
    expr_de_data_sub<-expr_de_data[expr_de_data$AveExpr >2,]
    
    expr_de_data_sub<-expr_de_data_sub[order(expr_de_data_sub$P.Value),]
    # 
    expr_de_data_sub<-expr_de_data_sub[1:100,]
    # 
    rownames(expr_de_data_sub)<-expr_de_data_sub$geneSymbol
    # 
    # expr_de_data_sub[,1]<-NULL
    # 
    samples_interest_fin<-paste0(samples_interest,".x")
    expr_data_in <- expr_de_data_sub[,samples_interest_fin]
    colnames(expr_data_in)<-gsub(".x","",colnames(expr_data_in))
    # 
    col_df<-data.frame(Expdes[c(i1,i2),])
    # 
    rownames(col_df)<-col_df[,1]
    # 
    col_df[,1]<-NULL
    # 
    col_col<-data.frame(t(col_df))
    
    #heatmaply(expr_data_in,grid_color = "black",scale="row",margins = c(200,300,100,100),col_side_colors=col_col) %>%
    #saveWidget(file="heatmaply_test.html",selfcontained = FALSE)
    
    # heatmaply(expr_data_in,grid_color = "black",scale="row",margins = c(200,300,100,100),col_side_colors=col_col)
    
    ngrps<-length(unique(Expdes[,2]))
    fname<-paste0(input$groups[1],input$groups[2],"top100_all_samples.pdf")
    expr_data_in_fin<-t(scale(t(expr_data_in)))
    
    rownames(Expdes)<-Expdes$sample
    Expdes$sample<-NULL
    pdf(file = fname,width = 16,height = 20)
    pheatmap(expr_data_in_fin,scale = "row",border_color = NA,cellwidth = 16,cellheight = 10,annotation_col = Expdes)
    dev.off()
    
    
    require(iheatmapr)
    expr_data_in_fin %>%
      as.matrix %>%
      iheatmap %>%
      add_col_labels(font = list(size = 8)) %>%
      add_row_labels(size = 0.3,font = list(size = 8))%>%add_col_clustering(k=ngrps)%>%
      add_row_clustering(k=2)  -> hm
    
    
   # hm %>% save_iheatmap("myplot2.pdf")
   # hm_final<-hm %>% as_plotly()
   # hm_final
    
    pheatmap(expr_data_in_fin,scale = "row",border_color = NA,cellwidth = 16,cellheight = 10,annotation_col = Expdes)
  })
  
  
  output$tsne<-renderPlot({
    if (is.null(input$gct) )
      return(NULL)
 
    voom_obj<-vobj()
    
    mdata<-voom_obj$E
    
    data_in<-data.frame(t(mdata))
    input_data<-input$gct
    arraydatafiles<-read.delim("arraydatafiles.txt",sep="\t",header = TRUE,row.names = 1)
    sampleData<-read.delim(as.character(arraydatafiles[input_data,]$expdesfile),stringsAsFactors = FALSE)
    
    mergData<-merge(sampleData,data_in,by.x="sample",by.y=0)
  set.seed(1985)
  
  #d <- stats::dist(mergData[,-c(1:2)])
  
  #rtsne_out <- Rtsne(d, is_distance=TRUE, perplexity=1, verbose = TRUE) 
  
  rtsne_out<-Rtsne::Rtsne(as.matrix(mergData[,-c(1:2)]),perplexity=3)
  
  rtsnedims<-rtsne_out$Y
  
  tsne_samples<-cbind(rtsnedims,sampleData)
  tsne_samples$condition<-as.character(tsne_samples$condition)
  
  str(tsne_samples)
  

  
  colnames(tsne_samples)<-c("x","y","Sample","Condition")
  
  rownames(tsne_samples)<-tsne_samples$Sample
 # ggscatter(tsne_samples, x = "x", y = "y",point = TRUE,label = "Sample",
 #           color = "Condition", shape = "Condition",size = 10,mean.point = FALSE,star.plot = FALSE,
 #           rug = FALSE, font.label = 9, repel = TRUE,label.rectangle = TRUE,palette = "aaas")
  
  ggscatter(tsne_samples, x = "x", y = "y",point = TRUE,size=6,
            color = "Condition",mean.point = FALSE,ellipse = FALSE,
            rug = FALSE, font.label = 9, repel = TRUE,label.rectangle = TRUE)+theme_pander() + scale_colour_pander()
  
  
  })
  
  output$pca<-renderPlot({
    if (is.null(input$gct) )
      return(NULL)
    voom_obj<-vobj()
    input_data<-input$gct
    mdata<-voom_obj$E
    
    data_in<-data.frame(t(mdata))
    
    arraydatafiles<-read.delim("arraydatafiles.txt",sep="\t",header = TRUE,row.names = 1)
    sampleData<-read.delim(as.character(arraydatafiles[input_data,]$expdesfile),stringsAsFactors = FALSE)
    mergData<-merge(sampleData,data_in,by.x="sample",by.y=0)
 
    pca_out<-PCA(mergData[,-c(1)],scale.unit=TRUE,ncp=5,graph=FALSE,quali.sup = 1)
    
    theme_custom<-theme(legend.position="top",axis.text.x = element_text(size=12,angle = 45, hjust = 1,face="bold"),plot.title = element_text(lineheight=.8,hjust=0.5, face="bold"),axis.text.y = element_text(size=12,face="bold"))+ theme(panel.border = element_rect(linetype = "solid",colour = "black"))
    
    
    #Amount of variance contributed for each PCA dimension
  #   eigenvalues <- pca_out$eig
     
   #  eigenvalues$dim<-1:5
     
  # p_bar<-ggplot(eigenvalues[,c(2,4)],aes(x=factor(dim),y=`percentage of variance`))+geom_bar(stat = "identity",fill="blue")+theme_bw()+ylim(0,50)+ylab("Dimension")+ggtitle("Variances")+xlab("Percentage of variances")
    # 
    # op <- par(mar = c(4, 4, 2, 2), oma=c(0, 0, 0, 0))
    # barplot(eigenvalues[1:10, 2], names.arg=1:10, 
    #         main = "Variances",
    #         xlab = "Principal Components",
    #         ylab = "Percentage of variances",
    #         col ="steelblue",ylim=c(0,50))
    # # Add connected line segments to the plot
    # lines(x = 1:10, eigenvalues[1:10, 2], 
    #       type="b", pch=19, col = "red")
    
    #You can do the scree plot simply by
    
   p_bar<-fviz_screeplot(pca_out, ncp=5)+theme_bw()
    
    #PCA plot for two dimensions for 
    
   # plot(pca_out,axes=c(1,2),choix='ind',habillage = 1)
    
    
    #Customize PCA plot using ggplot2
    
    pca_coords_df<-data.frame(pca_out$ind$coord)
    
    pca_coords_df$group<- mergData[,2]
    
    pca_coords_df$sample<-mergData[,1]
    
    pca_coords_df$group<-factor(pca_coords_df$group,levels=unique(pca_coords_df$group))
    
    pca_coords_df$sample<-factor(pca_coords_df$sample,levels=unique(pca_coords_df$sample))
    
    p_pca<-ggplot(pca_coords_df, aes(x=Dim.1, y=Dim.2,color=group),height=600, width=600)+theme_linedraw()+geom_point(size=8)+theme_custom
    
    grid.arrange(p_bar,p_pca,ncol = 2)
    
    
  })
  go_bp_enrich_up<-reactive({
    if (is.null(input$groups) )
      return(NULL)
    library(org.Mm.eg.db)
    library(org.Hs.eg.db)
    species<-input$species
    
    if(species=="Mm"){
      mapping_Symbol_eg<-as.list(org.Mm.egSYMBOL2EG)
    }
    else if (species=="Hs"){
      mapping_Symbol_eg<-as.list(org.Hs.egSYMBOL2EG)
    }
    
    
    gene_Table<- gene_stats()

    p_value <- input$p_value
    sig_table<-gene_Table[gene_Table$P.Value< p_value & gene_Table$logFC >0,]
    
    mapped_id_gene<-mapping_Symbol_eg[gene_Table$geneSymbol]
    
    mapped_id_gene_sig<-mapping_Symbol_eg[sig_table$geneSymbol]
    
    
    go.fisher <- goana(unlist(mapped_id_gene_sig),universe = unlist(mapped_id_gene), species=species)
    go_bp<-go.fisher[go.fisher$Ont=="BP",]
    
    go_bp[order(go_bp$P.DE),]
    
  })
  
  go_bp_enrich_down<-reactive({
    if (is.null(input$groups) )
      return(NULL)
    library(org.Mm.eg.db)
    library(org.Hs.eg.db)
    species<-input$species
    
    if(species=="Mm"){
      mapping_Symbol_eg<-as.list(org.Mm.egSYMBOL2EG)
    }
    else if (species=="Hs"){
      mapping_Symbol_eg<-as.list(org.Hs.egSYMBOL2EG)
    }
    
    
    gene_Table<- gene_stats()
    
    p_value <- input$p_value
    sig_table<-gene_Table[gene_Table$P.Value< p_value & gene_Table$logFC <0,]
    
    mapped_id_gene<-mapping_Symbol_eg[gene_Table$geneSymbol]
    
    mapped_id_gene_sig<-mapping_Symbol_eg[sig_table$geneSymbol]

    
    go.fisher <- goana(unlist(mapped_id_gene_sig),universe = unlist(mapped_id_gene), species=species)
    #go_bp<-go.fisher[go.fisher$Ont=="BP",]
    
    go.fisher[order(go.fisher$P.DE),]
    
  })
  
  output$plot1<-renderPlot({
    if (is.null(input$groups) )
      return(NULL)
    go_bp<-go_bp_enrich_up()
    go_bp$LogP<-log(go_bp$P.DE)
    go_bp<-go_bp[go_bp$N>20,]
    data_in <- transform(go_bp, Term = reorder(Term, -LogP))
    
    ggplot(data_in[1:10,], aes(x = Term, y = -LogP,fill="indianred") ) + geom_bar(stat="identity") +theme_bw()+coord_flip()+ xlab("GO Biological process") + ylab("- log P")+ggtitle("Enrichment UP")+theme(legend.position="none",axis.text.x = element_text(size=12,angle = 45, hjust = 1,face="bold"),plot.title = element_text(lineheight=.8,hjust=0.5, face="bold"),axis.text.y = element_text(size=12,face="bold"))
    
  })
  output$plot2<-renderPlot({
    if (is.null(input$groups) )
      return(NULL)
    go_bp<-go_bp_enrich_down()
    go_bp$LogP<-log(go_bp$P.DE)
    go_bp<-go_bp[go_bp$N>20,]
    data_in <- transform(go_bp, Term = reorder(Term, -LogP))
    
    ggplot(data_in[1:10,], aes(x = Term, y = -LogP,fill="blue") ) + geom_bar(stat="identity") +theme_bw()+coord_flip()+ xlab("GO Biological process") + ylab("- log P")+ggtitle("Enrichment DOWN")+theme(legend.position="none",axis.text.x = element_text(size=12,angle = 45, hjust = 1,face="bold"),plot.title = element_text(lineheight=.8,hjust=0.5, face="bold"),axis.text.y = element_text(size=12,face="bold"))
    
  })
  output$GO_table_UP <- renderDataTable({
    if (is.null(input$groups))
      return(NULL)
    go_Table<- go_bp_enrich_up()
    go_Table[go_Table$N>20,]
  })
  
  output$GO_table_DOWN <- renderDataTable({
    if (is.null(input$groups))
      return(NULL)
    go_Table<- go_bp_enrich_down()
    go_Table[go_Table$N>20,]
  })
  output$downloadData_GO <- downloadHandler(
    filename = function() {
      paste(input$groups[1],input$groups[2],input$adjusttype,"_GO_BP", '.csv', sep = '')
    },
    content = function(file) {
      GOBP_up<-go_bp_enrich_up()
      GOBP_down<-go_bp_enrich_down()
      GOBP_up$dir<-"Genes Up"
      GOBP_down$dir<-"Genes down"
      GOBP<-rbind(GOBP_up,GOBP_down)
      write.csv(GOBP, file)
    }
  )
 
}
shinyApp(ui = ui,server = server)