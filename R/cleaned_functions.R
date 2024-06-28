downloadableDT2 <- function(mytable, rownames=NULL, pageLength=10,...) {
  require(DT)
  if("pvalue" %in% colnames(mytable)){
    datatable(data = arrange(mytable, pvalue),
              rownames = rownames, ...,
              extensions = 'Buttons',
              options = list(
                dom = "Blfrtip",
                buttons = list("copy", list(
                  extend = "collection",
                  buttons = c("csv", "excel", "pdf"),
                  text = "Download"
                )), # end of buttons customization
                # customize the length menu
                lengthMenu = list( c(10, 20, -1), # declare values
                                   c(10, 20, "All") # declare titles
                ), # end of length Menu customization
                pageLength = pageLength)) # end of options
  }else{
    datatable(data = mytable,
              rownames = rownames, ...,
              extensions = 'Buttons',
              options = list(
                order = list(list(1, 'desc'), list(2, 'desc')), #comment out this line if you want to sort by the order of the table
                dom = "Blfrtip",
                buttons = list("copy", list(
                  extend = "collection",
                  buttons = c("csv", "excel", "pdf"),
                  text = "Download"
                )), # end of buttons customization
                # customize the length menu
                lengthMenu = list( c(10, 20, -1), # declare values
                                   c(10, 20, "All") # declare titles
                ), # end of length Menu customization
                pageLength = pageLength)) # end of options
  }
}

################################################################################

barplot_metrics <- function(mydf, metric = NULL, sort1 = NULL, sort2 = NULL, colorcol = NULL) {
  p <- mydf %>%
    # Arrange the data by timefactor
    arrange(.data[[sort1]], .data[[sort2]]) %>%
    # Modify the sample factor levels to follow the order of the sorts
    mutate(sample = fct_inorder(sample)) %>%
    # Create the plot
    ggplot(aes_string(x = "sample", y = metric, fill = colorcol)) +
    geom_bar(stat = "identity") +
    theme_bw() %+replace%
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size=6))
  return(p)
}

################################################################################

compute_cv <- function(x) sd(x) / mean(x)

################################################################################

pcaplot <- function(df, covar, pc1 = "PC1", pc2 = "PC2") {
  p <- df %>%
    ggplot(aes_string(x = pc1, y = pc2, color = covar, label = "sample")) +
    geom_point(size = 3, alpha=0.7) +
    xlab(paste0(pc1, " (percent variance explained ", round(pca_var[1], 2) * 100, "%)")) +
    ylab(paste0(pc2, " (percent variance explained ", round(pca_var[2], 2) * 100, "%)"))+
    theme(aspect.ratio=1) +
    theme(text = element_text(size = 16, ))+
    theme(legend.position="right")
  return(p)
}

################################################################################

plot_enhanced_volcano <- function(res_tib){
  p1<-EnhancedVolcano(res_tib,
                      lab = res_tib$symbol,
                      x='log2FC',
                      y='padj',
                      labSize = 4,
                      drawConnectors = F,
                      pCutoff = 5e-02,
                      subtitle = NULL,
                      title = NULL,
                      legendLabels=c('Not Sig','Sig Log2FC','Sig P-Value',
                                     'Sig P-Value & Log2FC'))
  return(p1)
}

################################################################################

plot_ma <- function(res_tib){
  plot_df <-
    res_tib %>%
    mutate(sig = case_when(padj < 0.05 & log2FC <= -0.585  ~ "padj < 0.05 & log2FC <= -0.585",
                           padj < 0.05 & log2FC >= 0.585 ~ "padj < 0.05 & log2FC >= 0.585",
                           TRUE ~ "padj > 0.05"),
           sig = factor(sig, levels = c("padj < 0.05 & log2FC <= -0.585", "padj > 0.05", "padj < 0.05 & log2FC >= 0.585")),
           abs_log2FC = abs(log2FC),
           log2fc_rank = rank(-abs_log2FC, ties.method = "first"),
           label = case_when(sig != "not sig" ~ symbol))

  top <- plot_df %>% filter(sig != "padj > 0.05") %>% top_n(-20, log2fc_rank)

  plot_df %>%
    ggplot(aes(x = log(baseMean), y = log2FC, color = sig)) +
    geom_point(alpha = 0.8, size = 0.5) +
    scale_color_manual(values=c(gg_color_hue(3)[3], "grey60", gg_color_hue(1)[1])) +
    labs(x = "log(baseMean)", color = "") +
    geom_text_repel(data = top, aes(x=log(baseMean), y=log2FC, label = label), size = 3, show.legend = F, color = "black", max.overlaps = 20) +
    theme_bw() +
    theme(legend.position = "top") + guides(color = guide_legend(override.aes = list(size = 3)))
}

################################################################################

display_de_genes3 <- function(mydat, gene_anns, filter = TRUE, distinct = FALSE, pval_cutoff = 0.05, lfc_cutoff = 0, LRT = FALSE) {
  mydat <- merge(geneAnns[, c(1,2, 3, 4,5, 6)], mydat, by.x = "geneid", by.y = "symbol", all.x = F, all.y = T) %>%
    dplyr::rename("log2FoldChange" = "log2FC") %>%
    filter(is.finite(padj))

  if (filter) {
    mydat <- mydat %>% filter(baseMean != 0, pvalue <= pval_cutoff & abs(log2FoldChange) >= lfc_cutoff)
  }
  if (distinct) {
    mydat <- mydat %>%
      distinct(geneid, .keep_all = TRUE)
  }
  format_cols <- c("log2FoldChange", "pvalue", "padj", "lfcSE", "baseMean")


  mydat %>%
    arrange(padj) %>% # This doesn't change the order of the table, it's sorted by the first and second columns per the downloadableDT function
    downloadableDT(filter = "bottom") %>%
    formatSignif(columns = format_cols, digits = 3)
}


################################################################################
gsea_results_go_kegg <- function(res_tib, geneset = "GO", myorgdb = "org.Hs.eg.db", org = "hsa", name, pval = 0.25, filter = NULL) {
  require(BiocParallel)
  require(parallel)

  geneset_name <- toupper(geneset)
  filename <- paste0("/gsea_res/", name, " ", geneset_name)
  dir.create(paste0(saveDir, "/gsea_res/"), showWarnings = F)

  if (!is.null(filter)) {
    res_tib <- res_tib %>% filter(padj < filter)
  }

  if (geneset_name == "GO") {
    # res_tib <-
    res_tib %>%
      left_join(., geneAnns, by = c("symbol" = "geneid"))     %>%
      select(log2FC, symbol)     %>%
      group_by(symbol) %>%
      summarize(log2FC = mean(log2FC)) %>%
      filter(symbol != "")
    generank <- makeGeneRankStat(res_tib$log2FC, res_tib$symbol)
    gseaRes <- gseGO(generank,
                     OrgDb = myorgdb,
                     keyType = "ALIAS",
                     ont = "BP",
                     pvalueCutoff = pval
    )
    gseaRes_table <- tibble(gseaRes@result)

    topPathwaysUp <- gseaRes_table %>%
      filter(NES > 0) %>%
      arrange(p.adjust, -NES) %>%
      head(10) %>%
      pull(ID)
    topPathwaysDown <- gseaRes_table %>%
      filter(NES < 0) %>%
      arrange(p.adjust, NES) %>%
      head(10) %>%
      pull(ID)
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

    gseaRes_table <- gseaRes_table %>%
      rename(
        pval = pvalue,
        pathway = Description,
        padj = p.adjust,
        ES = enrichmentScore,
        size = setSize
      )

    gseaRes_table <- gseaRes_table %>% mutate(link = paste0("<a href='https://amigo.geneontology.org/amigo/term/", ID, "'>", "Gene ontology info", "</a>"))
  } else if (geneset_name == "KEGG") {
    res_tib <- res_tib %>%
      left_join(., geneAnns, by = c("symbol" = "geneid")) %>%
      select(log2FC, entrezgene_id) %>%
      group_by(entrezgene_id) %>%
      summarize(log2FC = mean(log2FC)) %>%
      filter(entrezgene_id != "")
    generank <- makeGeneRankStat(res_tib$log2FC, res_tib$entrezgene_id)
    gseaRes <- gseKEGG(generank,
                       organism = org,
                       keyType = "kegg",
                       pvalueCutoff = pval
    )

    gseaRes_table <- tibble(gseaRes@result)

    topPathwaysUp <- gseaRes_table %>%
      filter(NES > 0) %>%
      arrange(p.adjust, -NES) %>%
      head(10) %>%
      pull(ID)
    topPathwaysDown <- gseaRes_table %>%
      filter(NES < 0) %>%
      arrange(p.adjust, NES) %>%
      head(10) %>%
      pull(ID)
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

    gseaRes_table <- gseaRes_table %>%
      rename(
        pval = pvalue,
        pathway = Description,
        padj = p.adjust,
        ES = enrichmentScore,
        size = setSize
      )

    gseaRes_table <- gseaRes_table %>% mutate(link = paste0("<a href='https://www.kegg.jp/entry/", ID, "'>", "KEGG info", "</a>"))
  }
  res <- list(
    gseares = gseaRes,
    gsea_table = gseaRes_table,
    topPathwaysUp = topPathwaysUp,
    topPathwaysDown = topPathwaysDown
  )

  saveRDS(res, file = file.path(saveDir, paste0("/gsea_res/", name, " ", geneset_name, ".RDS")))
  return(res)
}

################################################################################

plot_top_gsea <- function(gseaRes, top = "up"){
  gseares <- gseaRes$gseares
  topPathwaysUp <- gseaRes$topPathwaysUp
  topPathwaysDown <- gseaRes$topPathwaysDown
  if(top == "up"){
    if(length(topPathwaysUp) == 0){
      text = "No Pathways to Display"
      ggplot() +
        annotate("text", x = 4, y = 25, size=8, label = text) +
        theme_void()
    }else{
      pp <- lapply(topPathwaysUp, function(i) {

        anno <- gseares[i, c("NES", "pvalue", "p.adjust")] %>% mutate(NES=round(NES, 3),
                                                                      pvalue=formatC(pvalue, format = "e", digits = 2),
                                                                      p.adjust=formatC(p.adjust, format = "e", digits = 2))
        lab <- paste0(names(anno), "=",  anno, collapse="\n")

        gseaplot(gseares, i, gseares[i, 2], by="preranked") + xlab(NULL) +ylab(NULL) +
          annotate("text", 0, gseares[i, "enrichmentScore"] * .9, label = lab,
                   y= 1, x =Inf, vjust=0, hjust=1, size=2.75) + theme_minimal() %+replace% theme(title = element_text(size = 6))
      })
      cowplot::plot_grid(plotlist=pp, ncol=2)
    }
  }
  else if(top == "down"){
    if(length(topPathwaysDown) == 0){
      text = "No Pathways to Display"
      ggplot() +
        annotate("text", x = 4, y = 25, size=8, label = text) +
        theme_void()
    }else{
      pp <- lapply(topPathwaysDown, function(i) {

        anno <- gseares[i, c("NES", "pvalue", "p.adjust")] %>% mutate(NES=round(NES, 3),
                                                                      pvalue=formatC(pvalue, format = "e", digits = 2),
                                                                      p.adjust=formatC(p.adjust, format = "e", digits = 2))
        lab <- paste0(names(anno), "=",  anno, collapse="\n")

        gseaplot(gseares, i, gseares[i, 2], by="preranked") + xlab(NULL) +ylab(NULL) +
          annotate("text", x=-Inf,y=-0.25,hjust=0,vjust=1, label = lab, size=2.75, hjust=0) + theme_minimal() %+replace% theme(title = element_text(size = 6))
      })
      cowplot::plot_grid(plotlist=pp, ncol=2)
    }
  }
}

################################################################################

downloadableDT <- function(mytable, rownames = NULL, pageLength = 10, ...) {
  require(DT)
  if ("pvalue" %in% colnames(mytable)) {
    col_num <- which(colnames(mytable) == "pvalue")
    datatable(
      data = mytable,
      rownames = rownames, ...,
      extensions = "Buttons",
      options = list(
        order = list(list(col_num, "asc"), list(1, "asc")),
        dom = "Blfrtip",
        buttons = list("copy", list(
          extend = "collection",
          buttons = c("csv", "excel", "pdf"),
          text = "Download"
        )), # end of buttons customization
        # customize the length menu
        lengthMenu = list(
          c(10, 20, -1), # declare values
          c(10, 20, "All") # declare titles
        ), # end of length Menu customization
        pageLength = pageLength
      )
    ) # end of options
  } else {
    datatable(
      data = mytable,
      rownames = rownames, ...,
      extensions = "Buttons",
      options = list(
        order = list(list(1, "desc"), list(2, "desc")), # comment out this line if you want to sort by the order of the table
        dom = "Blfrtip",
        buttons = list("copy", list(
          extend = "collection",
          buttons = c("csv", "excel", "pdf"),
          text = "Download"
        )), # end of buttons customization
        # customize the length menu
        lengthMenu = list(
          c(10, 20, -1), # declare values
          c(10, 20, "All") # declare titles
        ), # end of length Menu customization
        pageLength = pageLength
      )
    ) # end of options
  }
}

################################################################################

gsea_results2 <- function(res_tib, geneset = hallmarks, name, filter = NULL) {
  require(BiocParallel)
  require(parallel)

  geneset_name <- deparse(substitute(geneset))
  filename <- paste0("/gsea_res/", name, " ", geneset_name)
  dir.create(paste0(saveDir, "/gsea_res/"), showWarnings = F)

  if (file.exists(paste0(saveDir, filename, ".RDS"))) {
    res <- read_rds(paste0(saveDir, filename, ".RDS"))
  } else {
    if (!is.null(filter)) {
      res_tib <- res_tib %>% filter(padj < filter)
    }
    generank <- makeGeneRankStat(res_tib$log2FC, res_tib$symbol)

    gseares <- GSEA(generank,
                    minGSSize = 5,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    TERM2GENE = geneset,
                    eps=0
    )

    gseares_table <- tibble(gseares@result)

    topPathwaysUp <- gseares_table %>%
      filter(NES > 0) %>%
      arrange(p.adjust, -NES) %>%
      head(10) %>%
      pull(Description)
    topPathwaysDown <- gseares_table %>%
      filter(NES < 0) %>%
      arrange(p.adjust, NES) %>%
      head(10) %>%
      pull(Description)
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

    gseares_table <- gseares_table %>%
      rename(
        pval = pvalue,
        pathway = Description,
        padj = p.adjust,
        ES = enrichmentScore,
        size = setSize
      )
    if (geneset_name == "keggs") {
      gseares_table <- gseares_table %>%
        left_join(select(keggs_info, id, term), by = c("pathway" = "term")) %>%
        distinct() %>%
        mutate(ID = id) %>%
        select(-id)
      gseares_table <- gseares_table %>% mutate(link = paste0("<a href='https://www.kegg.jp/entry/", ID, "'>", "KEGG info", "</a>"))
    } else if (geneset_name == "hallmarks") {
      gseares_table <- gseares_table %>% mutate(link = paste0("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", pathway, ".html", "'>", "mSigDB \n info", "</a>"))
    }
    res <- list(
      gseares = gseares,
      gsea_table = gseares_table,
      topPathwaysUp = topPathwaysUp,
      topPathwaysDown = topPathwaysDown
    )

    saveRDS(res, file = file.path(saveDir, paste0("/gsea_res/", name, " ", geneset_name, ".RDS")))
  }
  return(res)
}

################################################################################

sig_genes <- function(res_tib, direction = "up"){
  if(direction == "up"){
    res <- res_tib %>% filter(pvalue < 0.05 & log2FC >= 0.585) %>% distinct(symbol, .keep_all = T)
  }else if(direction == "down"){
    res <- res_tib %>% filter(pvalue < 0.05 & log2FC <= -0.585) %>% distinct(symbol, .keep_all = T)
  }else if(direction == "both"){
    res_down <- res_tib %>% filter(pvalue < 0.05 & log2FC <= -0.585) %>% distinct(symbol, .keep_all = T)
    res_up <- res_tib %>% filter(pvalue < 0.05 & log2FC >= 0.585) %>% distinct(symbol, .keep_all = T)
    res <- rbind(res_down, res_up) %>% distinct(symbol, .keep_all = T)
  }
  return(res)
}


################################################################################

plot_euler <- function(euler){
  require(eulerr)
  plot(euler(euler),
       quantities = list(type = c("counts", "percent"), fontsize=8
       ),
       lty = 1:3,
       labels = list(fontsize = 8),
       shape = "ellipse")
}

################################################################################

calculate_fig_height <- function(mat) {
  base_height <- 4  # Minimum height in inches
  per_row_height <- 0.2  # Additional height per row in inches
  height <- base_height + nrow(mat) * per_row_height
  return(height)
}

################################################################################


makeGeneRankStat <- function(stat, names) {
  #Create a list of statistics, with gene names as list names
  gene_rank <- stat
  names(gene_rank) <- names

  #Remove genes with no name, or those with NAs, NANs and Infs for statistics
  gene_rank <- gene_rank[which(!is.na(gene_rank))]
  gene_rank <- gene_rank[which(!is.nan(gene_rank))]
  gene_rank <- gene_rank[which(names(gene_rank) != "")]
  gene_rank[which(gene_rank == "-Inf")] <- min(gene_rank[which(gene_rank != "-Inf")], na.rm=T) -
    (0.01*min(gene_rank[which(gene_rank != "-Inf")], na.rm=T))
  gene_rank[which(gene_rank == "Inf")] <- max(gene_rank[which(gene_rank != "Inf")], na.rm=T) +
    (0.01*max(gene_rank[which(gene_rank != "Inf")], na.rm=T))

  #Remove duplicates with less-extreme values, first by ordering
  #the list by most-extreme absolute value, then removing
  #the second instance of any duplicate, with R's standard duplicated function
  gene_rank <- gene_rank[order(abs(gene_rank), decreasing=T)]
  gene_rank <- gene_rank[!duplicated(names(gene_rank))]

  #Set order, largest to smallest
  gene_rank <- gene_rank[order(gene_rank, decreasing=T)]
  return(gene_rank)
}

################################################################################


getsigmarker <- function(pval) {
  if(pval < 0.001)
    return("***")
  if(pval < 0.01) return("**")
  if(pval < 0.05) return("*")
  if(pval >= 0.05) return("")
  return("")
}
