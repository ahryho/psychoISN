library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag) #devtools::install_github("mattflor/chorddiag")

#' Title
#'
#' @param net 
#' @param cpg_loc 
#' @param snp_loc 
#' @param thrsh 
#'
#' @return
#' @export
#'
#' @examples
CustomChordDiagram <- function(net, cpg_loc, snp_loc, thrsh = 0.8){
  
  # Get CpGs' and SNPs' locations for only features present in the network "net"
  net_cpg_loc <- cpg_loc[CpG_ID %in% colnames(net)] %>% mutate(ID = CpG_ID) %>% dplyr::select(-CpG_ID)
  net_snp_loc <- snp_loc[SNP %in% colnames(net)] %>% mutate(ID = SNP) %>% dplyr::select(-SNP)
  net_loc     <- rbind(net_cpg_loc, net_snp_loc)[match(rownames(net), ID)]
  
  # Take only QTLs
  qtl_net     <- net[rownames(net) %in% net_snp_loc$ID, colnames(net) %in% net_cpg_loc$ID] %>% 
    as.matrix() %>% as.data.frame()
  
  # Transform net to long format 
  data_long <- qtl_net %>% 
    rownames_to_column %>%
    gather(key = 'key', value = 'value', -rowname)
  
  # Add  CpGs' and SNPs' locations and
  # change columns' names to make it more readable: column "key" is for  CpGs' chr, column "rowname" is for SNPs chr
  data_long <- left_join(data_long, net_loc, by = c("key" = "ID")) %>% mutate(cpg_chr = chr) %>% dplyr::select(-chr)
  data_long <- left_join(data_long, net_loc[, c("ID", "chr")], by = c("rowname" = "ID")) %>% mutate(snp_chr = chr) %>% dplyr::select(-chr)
  
  # Subset "net" on the given threshold and 
  # count the number of associations
  data_plt <- data_long[abs(data_long$value) >= thrsh, c("cpg_chr", "snp_chr", "value")] %>% setDT() %>% count(snp_chr, cpg_chr)
  data_plt <- data_plt[, lapply(.SD, as.numeric)]
  data_plt <- data_plt[order(data_plt$cpg_chr),]
  
  # Set-up parameters for plot
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 2, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
  par(mar = rep(0, 4))
  
  # Prepare color palette
  nr_color <- max(length(unique(data_plt$snp_chr)), length(unique(data_plt$cpg_chr)))
  cpalette <- viridis(nr_color, alpha = 1, begin = 0, end = 1, option = "G")
  cpalette <- cpalette[sample(1:nr_color)]
  
  # Base plot
  cdm <- chordDiagram(
    x = data_plt, 
    grid.col = cpalette,
    transparency = 0.25,
    directional = 1,
    direction.type = c("arrows", "diffHeight"), 
    diffHeight  = -0.04,
    annotationTrack = "grid", 
    annotationTrackHeight = c(0.05, 0.1),
    link.arr.type = "big.arrow", 
    link.sort = TRUE, 
    link.largest.ontop = TRUE)
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    if(abs(CELL_META$cell.start.degree - CELL_META$cell.end.degree) > 3) {
      sn = CELL_META$sector.index
      i_chr = as.numeric(gsub("(C|R)_", "", sn))
      circos.text(CELL_META$xcenter, CELL_META$ycenter, i_chr, col = "white", 
                  font = 2, cex = 0.7, adj = c(0.5, 0.5), niceFacing = TRUE)
      # xlim = CELL_META$xlim
      # breaks = seq(0, xlim[2], by = 4e5)
      # circos.axis(major.at = breaks, labels = paste0(breaks/1000, "KB"), labels.cex = 0.5)
    }
  }, bg.border = NA)
}
