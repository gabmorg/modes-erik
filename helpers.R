# helpers.R
# data objects and plots for the shiny app

library(readr)
library(reshape)
library(RColorBrewer)
library(ggplot2)  
library(dplyr)


#####
# The helper functions
#####
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
gg_color_hue <- function(n) {
  ## ggplot default rainbow colors from: https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette 
  # hues = seq(15, 375, length = n + 1)
  hues = seq(25, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

######
# The variables and data objects
#####

tum.cols <- gg_color_hue(21)
reg.cols <- brewer.pal(n = 6,name = "Dark2")

# if testing script:
# r2.added.list <- readRDS("./myapp/data/r2.added.list.rds")

# if launching app.R:
# r2.added.list <- readRDS("./myapp/data/r2.added.list.rds")
# mapped_top10_lncrna_coefs <- readRDS("./myapp/data/mapped_top10_lncrna_coefs.rds")
# mapped_top10_mirna_coefs <- readRDS("./myapp/data/mapped_top10_mirna_coefs.rds")
# mapped_top10_tf_coefs <- readRDS("./myapp/data/mapped_top10_tf_coefs.rds")

cancer_names <- names(r2.added.list)

# r2.added.list <- readRDS("data/r2.added.list.rds")

#####
# The plots
#####

# mygene = string
# cond.vars = int vector, length 1 to 6 (optional)
# uncond.vars = int vector, length 1 to 6 (optional)

# TODO: rewrite so the ggplot object is made on front-end, 
# and this function returns only the right dataframe (does all the filtering)
cond.vars = c(1,3,4) # conditions
uncond.vars = c(1,3,4) # exclusions

genePlot <- function(mygene = "TRPM1", 
                     cond.vars = NULL, 
                     uncond.vars = NULL) {
  
  gene.list <- lapply(r2.added.list, function(x){
    if(mygene %in% names(x)){
      newdf <- x[[mygene]]
      rownames(newdf) <- 1:nrow(newdf)
      return(newdf)
    } else{
      return(NA)
    }
  })
  
  # Default displayed data
  full.gene.df <- Reduce(rbind, gene.list[names(gene.list[is.na(gene.list) == FALSE])])
  gene_spec_cancers = names(gene.list[is.na(gene.list) == FALSE])
  full.gene.df$cancer <- rep(gene_spec_cancers,each = nrow(gene.list[is.na(gene.list) == FALSE][[1]]))
  
  
  defaultPlot <- ggplot(full.gene.df, aes(x = var.added.f, y = r2.growth, fill = cancer)) + 
    geom_boxplot() + xlab("Regulator") + ylab("R2 added") + ylim(c(-1,1)) + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")) + 
    scale_fill_manual(values = tum.cols)
  
  if(!is.null(cond.vars)) {
    # needs testing
    cond.data <- c(1:length(cond.vars))
    cond.data <- lapply(cond.vars,
                        function(x) {
                          print(c(x, which(cond.vars == x)))
                          i <- which(cond.vars == x)
                          cond.data[i] <- as.vector(grep(as.character(x), 
                                                         full.gene.df$base.model))})
    
    cond.gene.df <- full.gene.df[Reduce(intersect, cond.data), ]
    
    condPlot <- ggplot(cond.gene.df, aes(x = var.added.f, y = r2.growth, fill = cancer)) +
      geom_boxplot() + xlab("Regulator") + ylab("R2 added") + ylim(c(-1,1)) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold")) +
      scale_fill_manual(values = tum.cols)
    
    return(condPlot)
  }
  
  else if(!is.null(uncond.vars)) {
    
    uncond.data <- c(1:length(uncond.vars))
    uncond.data <- lapply(uncond.vars,
                          function(x) {
                            print(c(x, which(uncond.vars == x)))
                            i <- which(uncond.vars == x)
                            uncond.data[i] <- as.vector(grep(as.character(x), 
                                                             full.gene.df$base.model,
                                                             invert = TRUE))})
    
    uncond.gene.df <- full.gene.df[Reduce(intersect, uncond.data), ]
    
    uncondPlot <- ggplot(uncond.gene.df, aes(x = var.added.f, y = r2.growth, fill = cancer)) +
      geom_boxplot() + xlab("Regulator") + ylab("R2 added") + ylim(c(-1,1)) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold")) +
      scale_fill_manual(values = tum.cols)
    
    return(uncondPlot)
  }
  
  else {
    return(defaultPlot)
  }
  
}

# regulator type: one of mirna, tf, lncrna
coefPlot <- function(regulator_type, mygene, plotType) {
  
  # get all the data for a single gene
  # colour code by the cancer
  # divide all values by 21 (number of cancers)
  gene_df <- NULL
  
  if(regulator_type == "lncrna") {
    mygene_coefs_list <- unlist(mapped_top10_lncrna_coefs[[mygene]])
    mygene_coefs_df <- data.frame(stringsAsFactors = F, 
                                  cancer.reg = names(mygene_coefs_list),
                                  weighted_average = mygene_coefs_list/21 #TODO: confirm this calc
    )
    mygene_coefs_df <- tidyr::separate(mygene_coefs_df, "cancer.reg", 
                                       into = c("cancer", "regulator_id", "ensg_version_id"),
                                       remove=T, sep = "\\.")
    gene_df <- mygene_coefs_df
    
  }
  
  else if(regulator_type == "mirna") {
    mygene_coefs_list <- unlist(mapped_top10_mirna_coefs[[mygene]])
    mygene_coefs_df <- data.frame(stringsAsFactors = F, 
                                  cancer.reg = names(mygene_coefs_list),
                                  weighted_average = mygene_coefs_list/21 #TODO: confirm this calc
    )
    mygene_coefs_df <- tidyr::separate(mygene_coefs_df, "cancer.reg", 
                                       into = c("cancer", "regulator_id"),
                                       remove=T, sep = "\\.")
    gene_df <- mygene_coefs_df
    
  }
  
  else if(regulator_type == "tf") {
    mygene_coefs_list <- unlist(mapped_top10_tf_coefs[[mygene]])
    mygene_coefs_df <- data.frame(stringsAsFactors = F, 
                                  cancer.reg = names(mygene_coefs_list),
                                  weighted_average = mygene_coefs_list/21 #TODO: confirm this calc
    )
    mygene_coefs_df <- tidyr::separate(mygene_coefs_df, "cancer.reg", 
                                       into = c("cancer", "regulator_id", "ensg_version_id"),
                                       remove=T, sep = "\\.")
    gene_df <- mygene_coefs_df
    
  }
  
  else {
    errorCondition("invalid regulator ID: tf, lncrna, mirna")
  }
  
  coef_barplot <- ggplot(gene_df, aes(fill = cancer, y = weighted_average, x = regulator_id)) + 
    geom_bar(position=plotType, stat="identity") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")) + 
    scale_fill_manual(values = tum.cols) +
    ggtitle(paste(mygene, "top 10", regulator_type)) 
  
  return(coef_barplot)
}

####### coefPlot testing code #######
# 
# scyl_tf_coefs <- unlist(mapped_top10_tf_coefs$SCYL3)
# scyl_tf_coefs <- data.frame(stringsAsFactors = F,
#                              cancer.reg = names(scyl_tf_coefs),
#                              weighted_average = scyl_tf_coefs/21 #TODO: confirm this calc
#                               )
# scyl_tf_coefs <- tidyr::separate(scyl_tf_coefs, "cancer.reg",
#                                   into = c("cancer", "regulator_id", "ensg_version_id"),
#                                    remove=T, sep = "\\.")
# 
# coef_stackedplot <- ggplot(scyl_tf_coefs, aes(fill = cancer, y = weighted_average, x= regulator_id)) + 
#   geom_bar(position="stack", stat="identity") +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold")) + 
#   scale_fill_manual(values = tum.cols) +
#   ggtitle(paste("SCYL3", "top 10", "regulator_type"))
