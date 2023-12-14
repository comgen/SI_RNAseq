
#a collection of functions created to analyse the RNAseq data for  correlation and stochiometric imbalance
#fixing default function names
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
relocate <- dplyr::relocate
colnames <- base::colnames
rownames <- base::rownames

symbol_filterLow <- function(x) { 
  protein_ensg <- grch38 %>% filter(biotype=="protein_coding") %>% select(ensgene) #this step is dependent on the package annotables, which may not be updated
  x <- x[base::rownames(x) %in% protein_ensg$ensgene,]
  x <- x[rowSums(x > 1) >= length(base::colnames(x))/2, ] #filter out genes with >50% <1 TPM
  #replace ensembl IDs with gene symbols
  x <- cbind(symbol = grch38$symbol[ match(rownames(x), grch38$ensgene) ],
             x)
  x$symbol[duplicated(x$symbol)] <- NA #mark the duplicated gene names as NA
  x <- tibble::as_tibble(x, rownames = "gene_id") %>%  
    mutate(symbol = coalesce(symbol, gene_id)) #insert ensembl id if NA
  x <- as.data.frame(x[,-1]) #remove the ensemble ids
  rownames(x) <- x$symbol #make the gene names rownames
  x <- x[,-1] #remove the gene name column
  x
}


#####functions to create heatmap of gene expression correlation and order the genes according to another set
cluster.heatmap <- function(TPM_subset){
  TPM_prot_t <- t(TPM_subset)
  #the default color palette used in corrplot() 
  col2 <- colorRampPalette(c( "#053061", "#2166AC", "#4393C3",
                              "#92C5DE","#D1E5F0", "#FFFFFF",  
                              "#FDDBC7",  "#F4A582", "#D6604D", 
                              "#B2182B", "#67001F"))
 
  M <- cor(TPM_prot_t, method = "pearson")
  
  #Correlation is a measure of similarity so we’ll turn it to into a measure of dissimilarity before passing it to the as.dist
  spellman.dist <- as.dist(1 -M)
  
  #We first generate the hierarchical clustering, use the “average" method
  spellman.tree <- hclust(spellman.dist, method="average")
  spellman.dend <- as.dendrogram(spellman.tree) # create dendrogram object
  
  heatmap <- heatmap.2(M, 
                       Rowv = ladderize(spellman.dend), 
                       Colv = ladderize(spellman.dend), 
                       dendrogram = "column", 
                       #revC = TRUE,  #have previously reversed the columns, but this caused some confusion when I later wanted to force the order of other heatmaps.
                       trace = "none", 
                       density.info = "none",
                       breaks= seq(-1,1,length.out=16),
                       col = col2,
                       key = TRUE,
                       cexRow= 0.5,
                       cexCol = 0.5,
                       labRow = base::rownames(M), labCol = base::colnames(M)) 
}

#the default color palette used in corrplot() 
col2 <- colorRampPalette(c( "#053061", "#2166AC", "#4393C3",
                            "#92C5DE","#D1E5F0", "#FFFFFF",  
                            "#FDDBC7",  "#F4A582", "#D6604D", 
                            "#B2182B", "#67001F"))

heatmap.compare <- function(TPM_subset, order){
  TPM_prot_t <- t(TPM_subset)
  M <- cor(TPM_prot_t, method = "pearson")
 # M <- M[rownames(M)[order], colnames(M)[order]]
  M <- M[match(order, rownames(M)), match( order, colnames(M))]
  
  #Correlation is a measure of similarity so we’ll turn it to into a measure of dissimilarity before passing it to the as.dist
  spellman.dist <- as.dist(1 -M)
  heatmap <- heatmap.2(M, 
                       Rowv = FALSE, 
                       Colv = FALSE, 
                       dendrogram = "none", 
                       trace = "none", 
                       density.info = "none",
                       breaks= seq(-1,1,length.out=16),
                       col = col2, 
                       key = FALSE,
                       revC = TRUE,
                       cexRow = 0.5,
                       #cexCol = 0.5,
                       lwid=c(0.1,4), lhei=c(0.1,4),
                       labRow = rownames(M), 
                       labCol =  colnames(M)) 
}



##### ROC plot using grid
#geom and function written by Koundinya Desiraju @koundy 
GeomROCplot <- ggproto(
  
  "GeomROCplot", Geom,
  required_aes = c("x", "y"),
  
  default_aes = aes(colour = "black",
                    size = 1,
                    linetype = 1,
                    alpha = 1,
                    fill = NA,
                    diagonal =TRUE,
                    diagonal_colour = 'black',
                    diagonal_alpha = 0.8,
                    diagonal_linetype = 2,
                    diagonal_width = 0.5 ),
  
  draw_key = draw_key_abline,
  
  draw_group = function(data, panel_scales, coord) {
    
    n <- nrow(data)
    if (n <= 2) return(grid::nullGrob())
    
    coords <- coord$transform(data, panel_scales)
    
    first_row <- coords[1, , drop = FALSE]
    
    if(first_row$diagonal){
      
      grid::gList(
        
        grid::linesGrob(
          x= coords$x,
          y=coords$y,
          default.units = "native",
          gp = grid::gpar(
            
            col = scales::alpha(
              first_row$colour,
              first_row$alpha),
            
            lwd = first_row$size * .pt,
            lty = first_row$linetype)),
        
        grid::linesGrob(
          x=c(0,1),
          y=c(0,1),
          default.units = "native",
          gp = grid::gpar(
            
            col = scales::alpha(
              first_row$diagonal_colour, 
              first_row$diagonal_alpha),
            
            lwd = first_row$diagonal_width * .pt,
            lty = first_row$diagonal_linetype
          )
        )
      )
    }
    
    else{
      grid::linesGrob(
        x= coords$x,
        y=coords$y,
        default.units = "native",
        gp = grid::gpar(
          col = scales::alpha(
            first_row$colour, 
            first_row$alpha),
          
          lwd = first_row$size * .pt,
          lty = first_row$linetype)
      )
    }
  }
)


geom_roc_plot <- function(
  
  mapping = NULL, 
  data = NULL,
  stat = "identity",
  position = "identity", 
  na.rm = FALSE, 
  show.legend = NA, 
  inherit.aes = TRUE, ...) {
  
  
  layer(
    geom = GeomROCplot, 
    mapping = mapping,
    data = data, 
    stat = stat,
    position = position, 
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
  
}

