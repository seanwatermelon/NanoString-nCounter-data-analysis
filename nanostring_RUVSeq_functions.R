RUV_total <- function(raw,pData,fData,k,hkgenes = NULL,exclude = NULL){
  
  ### INPUT: raw - p x n raw expressions with p genes and n samples
  ###        pData - phenotype metadata across samples
  ###        fData - feature metadata across genes
  ###        k - number of dimensions of unwanted variation estimated
  ###        exclude - vector of gene names to exclude
  
  library(RUVSeq)
  library(DESeq2)
  library(limma)
  library(matrixStats)
	
  if (!is.null(hkgenes)){
	  
  fData(set)$Class[rownames(set) %in% hkgenes] = 'Housekeeping'
	  
	  }
  
  fData = fData[rownames(raw),]
  int = intersect(rownames(raw),rownames(fData))
  fData = fData[int,]
  raw = raw[int,]

  set <- newSeqExpressionSet(as.matrix(round(raw)),
                             phenoData=pData,
                             featureData=fData)

  cIdx <- rownames(set)[fData(set)$Class == "Housekeeping"]
  cIdx = cIdx[!(cIdx %in% exclude)]
  x <- as.factor(pData$Group)
  set <- betweenLaneNormalization(set, which="upper")
  set <- RUVg(set, cIdx, k=k)
  dds <- DESeqDataSetFromMatrix(counts(set),colData=pData(set),design=~1)
  rowData(dds) <- fData
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  cts <- counts(dds, normalized=TRUE)
  disp <- pmax((rowVars(cts) - rowMeans(cts)),0)/rowMeans(cts)^2
  mcols(dds)$dispGeneEst <- disp
  dds <- estimateDispersionsFit(dds, fitType="mean")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  mat <- assay(vsd)
  covars <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))),drop=FALSE])
  mat <- removeBatchEffect(mat, covariates=covars)
  assay(vsd) <- mat
  return(list(set = set,vsd = vsd))



}


imagingQC <- function(rcc){
  
  
  #### INPUT: rcc - input from rcc
  #### OUTPUT: flag for imaging quality

	fovRatio = as.numeric(rcc$Lane_Attributes[3]) / as.numeric(rcc$Lane_Attributes[2])
	if (!(fovRatio > .75)) {return('Flag')}
	if (fovRatio > .75) {return('No flag')}

}


  
imagingQC_value <- function(rcc){
  
  
  #### INPUT: rcc - input from rcc
  #### OUTPUT: flag for imaging quality
  
  fovRatio = as.numeric(rcc$Lane_Attributes[3]) / as.numeric(rcc$Lane_Attributes[2])
  return(fovRatio)
  
}


bindingDensityQC <- function(rcc,low,high){
  
  
  #### INPUT: rcc - input from rcc
  ####         low, high - the lower and upper limits for binding density
  #### OUTPUT: flag for binding density

	bd = as.numeric(rcc$Lane_Attributes[6])
	if(!(bd < high & bd > low)) {return('Flag')}
	if (bd < high & bd > low) {return('No flag')}
	

}


bindingDensityQC_value <- function(rcc){
  
  
  #### INPUT: rcc - input from rcc
  ####         low, high - the lower and upper limits for binding density
  #### OUTPUT: flag for binding density
  
  bd = as.numeric(rcc$Lane_Attributes[6])
  return(bd)
  
  
}



limitOfDetectionQC <- function(rcc,numSD = 1){

  #### INPUT: rcc - input from rcc
  ####         numSD - number of standard deviations to calibrate the LOD
  #### OUTPUT: flag for limit of detection
  
	counts = rcc$Code_Summary
	posE = as.numeric(counts$Count[counts$Name == 'POS_E'])
	negControls = as.numeric(counts$Count[grepl('NEG',counts$Name)])
	if(!(posE > mean(negControls) + numSD*sd(negControls))) {return('Flag')}
	if (posE > mean(negControls) + numSD*sd(negControls)) {return('No flag')}

}


  
limitOfDetectionQC_value <- function(rcc,numSD = 1){
  
  #### INPUT: rcc - input from rcc
  ####      numSD - number of standard deviations to calibrate the LOD
  #### OUTPUT: flag for limit of detection
  
  counts = rcc$Code_Summary
  negControls = as.numeric(counts$Count[grepl('NEG',counts$Name)])
  LOD = mean(negControls) + numSD*sd(negControls)
  return(LOD)
  
}



positiveLinQC <- function(rcc){

  #### INPUT: rcc - input from rcc
  #### OUTPUT: flag for linearity for positive controls
  
  
	counts = rcc$Code_Summary
	posControls = as.numeric(counts$Count[grepl('POS_',counts$Name)])
	known = c(128,128/4,128/16,128/64,128/256,128/(256*4))
	r2 = summary(lm(sort(posControls)~sort(known)))$r.squared
	if(!(r2 > .95) | is.na(r2)) {return('Flag')}
	if(r2 > .95) {return('No flag')}

}



positiveLinQC_value <- function(rcc){
  
  #### INPUT: rcc - input from rcc
  #### OUTPUT: value for linearity for positive controls
  
  
  counts = rcc$Code_Summary
  posControls = as.numeric(counts$Count[grepl('POS_',counts$Name)])
  known = c(128,128/4,128/16,128/64,128/256,128/(256*4))
  r2 = summary(lm(sort(posControls)~sort(known)))$r.squared
  
  return(r2)
  
}



makeRLEplot <- function(data,metadata,id){
  
  #### INPUT: data - matrix of expressions with genes on rows and samples on columns
  ####        metadata - matrix of metadata with a column that corresponds to the colnames of data
  ####        id - colname of sample ids
  #### OUTPUT: ggplot2 RLE plot
  
  data = data - apply(data,1,median)
  stack = stack(data)
  colnames(stack)[1] = id
  stackPlot = merge(stack,metadata,by=id)
  colnames(stackPlot)[1:2] = c('Sample','values')
  rle_plots = ggplot(data = stackPlot,aes(x = Sample,y = values, color = ER_status)) +
    geom_boxplot(coef = 6) + theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=24),
          plot.title = element_text(size = 30),
          legend.title=element_text(size=20),
          legend.text=element_text(size=20),
          strip.text = element_text(size=24),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = "grey", 
                                      fill = NA, size = .1),
          legend.position = 'bottom',
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + xlab('Sample') +
    ylab('Median deviation of log expression') + ylim(c(-4,4))
  return(rle_plots)
  
}


volcanoplot_neo <- function (fit, coef = 1L, style = "p-value", highlight = 0L, 
                             names = fit$genes$ID, hl.col = "black", title = NULL, xlab = "Log2 Fold Change", 
                             ylab = NULL, pch = 16, cex = 0.7, legend = TRUE, ...) 
  {
    library(ggrepel)
    library(ggplot2)
    library(dplyr)
    if (!"package:ggrepel" %in% search()) {
      stop("ggplot2 is not loaded in the environment.")
    }
    if (!"package:ggplot2" %in% search()) {
      stop("ggplot2 is not loaded in the environment.")
    }
    if (!"package:dplyr" %in% search()) { 
      stop("dplyr is not loaded in the environment.")
    }
    if (!is(fit, "MArrayLM")) 
      stop("fit must be an MArrayLM")
    x <- as.matrix(fit$coef)[, coef]
    style <- match.arg(tolower(style), c("p-value", "b-statistic"))
    if (style == "p-value") {
      if (is.null(fit$p.value)) 
        stop("No p-values found in linear model fit object")
      y <- as.matrix(fit$p.value)[, coef]
      y <- -log10(y)
      if (is.null(ylab)) 
        ylab = "-log10(P-value)"
    }
    else {
      if (is.null(fit$lods)) 
        stop("No B-statistics found in linear model fit object")
      y <- as.matrix(fit$lods)[, coef]
      if (is.null(ylab)) 
        ylab = "Log Odds of Differential Expression"
    }
    df <- data.frame(x = x, y = y, names = names)
    
    # Mark genes with increased or decreased expression
    df$label <- ifelse(x > 1, paste0(names, "\n"), 
                       ifelse(x < -1, paste0(names, "\n"), ""))
    
    # Determine the upregulate
    df$diffexpressed <- "No"
    df$diffexpressed[df$x > 2 & df$y > 1.301] <- "Up"
    df$diffexpressed[df$x < -2 & df$y > 1.301] <- "Down"
    
    # Determine x-axis limits and breaks
    x_max <- max(abs(df$x)) + 0.5  # Find the maximum absolute value for x
    x_limits <- c(-ceiling(x_max), ceiling(x_max))
    x_breaks <- seq(-ceiling(x_max), ceiling(x_max), by = 1)  # Set breaks with interval of 2
    
    top_genes <- data.frame()
    
    # Highlight top significant genes
    if (highlight > 0) {
      # Subset for x > 0 and y maximum
      top_genes_up <- df %>% filter(x > 2) %>% arrange(desc(y)) %>% head(highlight)
      
      # Subset for x < 0 and y maximum
      top_genes_down <- df %>% filter(x < -2) %>% arrange(desc(y)) %>% head(highlight)
      
      # Combine both subsets
      top_genes <- rbind(top_genes_up, top_genes_down)
    
    }
    
    # Create the ggplot
    p <- ggplot(df, aes(x = x, y = y, label = label)) +
      geom_point(aes(color = diffexpressed), size = cex, shape = pch) +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "No" = "grey")
                         , labels = c("Up" = "Upregulated", "Down" = "Downregulated", 
                                      "No" = "Not significant")) +
      labs(x = xlab, y = ylab, color = "Expression Status", title = title) +
      theme_minimal() +
      geom_vline(xintercept = c(-2, 2), linetype = "dashed", col = "#A9A9A9") +
      geom_hline(yintercept = c(1.301), linetype = "dashed", col = "#A9A9A9") +
      scale_x_continuous(limits = x_limits, breaks = x_breaks) +
      scale_y_continuous(limits = c(0,10), breaks = seq(0,10,by=2)) +
      theme_bw() +
      theme(title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 12, face = "bold"), 
            legend.text = element_text(size = 10, face = "bold"),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold")) +
      geom_text_repel(data = top_genes, max.overlaps = Inf,
                      nudge_x = 0.02,
                      nudge_y = 0.02)
    
    if(!legend){
      p <- p + theme(legend.position = "none")
    }else{
      p <- p + theme(legend.position = "right")
    }
    
    invisible(p)
    
    return(p)
    
    
  }



