library(raster)
library(picante)
library(ape)
library(dplyr)
phylo_branch_matrix <- function(site_sp_matrix, tree_file) {
  sp_names <- colnames(site_sp_matrix)
  site_species_mat <- site_sp_matrix

  # Make branch by site matrix
  branches <- cbind(tree_file$edge, data.frame(edge_length=tree_file$edge.length))
  ext_branches <-branches[which(branches[, 2] <= length(tree_file$tip.label)), ]
  labs <- data.frame(tip.label=tree_file$tip.label)
  ext_branches.n <- merge(ext_branches, labs, by.x = '2', by.y = 0)
  int_branches <-branches[which (!branches[, 2] %in% 1:length(tree_file$tip.label)),]
  
  # Make species by node matrix
  dimnames <-list(ext_branches.n$"tip.label", int_branches$"2")
  tip_branch_mat <-matrix(0, nrow(ext_branches.n), nrow(int_branches), dimnames = dimnames)
  idx_first_branch <- min(int_branches$"2")
  for (brch in int_branches$"2") {
    tips.nums <- picante::internal2tips(tree_file, brch)
    tip_branch_mat[tips.nums, brch - idx_first_branch + 1] <- 1
  }
  
  # Reorder tip.branch to match the site.species
  tip_branch_mat.order <-tip_branch_mat[match(sp_names, rownames(tip_branch_mat)), ]
  site_branch <- site_species_mat %*% tip_branch_mat.order
  site_branch[site_branch[, 1] > 0, ]
  site_branch_ones <-(site_branch & TRUE) * 1 # This step turn matrix with probability values to presence-absence !
  
  # Replace speceis names by tips numbers
     newnames <- fortify(tree_file) %>% dplyr::select(node, label) %>% arrange(label) %>% na.omit() %>% data.frame()
     colnames(site_sp_matrix)
     site_sp_matrix <- site_sp_matrix[,newnames$label]
     colnames(site_sp_matrix) <- newnames$node
  
     site_sp_matrix <- data.frame(ncell=rownames(site_sp_matrix), site_sp_matrix)
     site_branch_ones <- data.frame(ncell=rownames(site_branch_ones), site_branch_ones)
     all_dist <- dplyr::left_join(site_sp_matrix, site_branch_ones, by='ncell')
  
     # Write raster of branches
     return(all_dist)
 }

  