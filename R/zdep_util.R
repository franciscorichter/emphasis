get_extant <- function(tm, tree) {
    origin <- setdiff(c(tree$extant$parent, tree$extinct$parent),
                      c(tree$extant$child, tree$extinct$child))
    # obtain extant tree from full tree
    if (nrow(tree$extinct) > 0) {
        if (tm == 0) {
            tm <- 1e-10
        }
        tree$extant$clade <- NULL
        extinct <- tree$extinct[tree$extinct$brts < tm &
                               tree$extinct$t_ext < tm, ]
        extinct <- extinct[order(extinct$brts), ]
        extant <- rbind(tree$extant[tree$extant$brts < tm, ],
                       tree$extinct[tree$extinct$brts < tm &
                                    tree$extinct$t_ext >= tm, -4])
        extant <- extant[order(extant$brts), ]
        extinct$t_ext <- NULL
        if (nrow(extinct) > 0) {
            i <- 1
            while (i <= nrow(extant)) {
                if ((sum(extant$parent[i] == extant$child) == 0) &&
                    (extant$parent[i] != origin)) {
                  ind <- which(extant$parent[i] == extinct$child)
                  ind2 <- which(extinct$parent == extant$parent[i] &
                               extinct$brts < extant$brts[i])
                  child <- extant$child[i]
                  new_extant_row <- extinct[ind, ]
                  new_extant_row[3] <- child
                  new_extinct_row <- extant[i, c(1, 3, 2)]
                  extant[i, ] <- new_extant_row
                  extinct[ind, ] <- new_extinct_row
                  extinct$parent[ind2] <- child
                } else {
                  i <- i + 1
                }
            }
            extant$parent[1] <- origin
            extant$parent[2] <- extant$child[1]
        }
    } else {
        extant <- tree$extant
    }
    extant <- extant[order(extant$brts), ]
    if (extant$parent[2] == origin) {
        extant[c(1, 2), ] <- extant[c(2, 1), ]
    }
    extant <- extant[extant$brts <= tm, ]
    return(extant)
}

transf <- function(name_spe, vec) {
    which(vec == name_spe)
}

newick <- function(tree, CT) {
    n <- nrow(tree)
    child.nms <- as.character(tree$child)
    parent.nms <- as.character(tree$parent)
    species.nms <- unique(child.nms, parent.nms)
    n.species <- length(species.nms)
    CT <- rep(CT, n.species)
    for (i in seq(n, 1)) {
        nw <- paste("(",
                    parent.nms[i],
                    ":",
                    as.character(CT[which(species.nms == parent.nms[i])] -
                                     tree$brts[i]),
                    ",",
                    child.nms[i],
                    ":",
                    as.character(CT[which(species.nms == child.nms[i])] -
                                     tree$brts[i]), ")", sep = "")
        j <- which(parent.nms[i] == child.nms)
        rp <- which(parent.nms == child.nms[j])
        if (length(rp) > 0) {
            parent.nms[rp] <- nw
        }
        species.nms[which(species.nms == child.nms[j])] <- nw
        child.nms[j] <- nw
        CT[j] <- CT[j] -
            (CT[which(species.nms == parent.nms[i])] - tree$brts[i])
    }
    return(paste(child.nms[1], ";", sep = ""))
}

PDTT_plot <- function(tree) {
    ct <- max(tree$brts)
    times <- seq(0, ct, length.out = 1000)
    times <- sort(times, decreasing = TRUE)
    PD <- NULL
    for (i in seq_along(times)) {
        G <- GPD(times[i], tree[-1, ])

        PD <- rbind(PD,
                   data.frame(time = rep(times[i], nrow(G)),
                              P = colSums(G) / (nrow(G) - 1),
                              lineage = as.character(seq_len(nrow(G)))))
    }
    g1 <- ggplot2::ggplot(PD) +
          ggplot2::geom_line(ggplot2::aes_string(x = "time",
                                                 y = "P",
                                                 colour = "lineage",
                                                 alpha = 0.5)) +
          ggplot2::theme(legend.position = "none",
                         panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         panel.background = ggplot2::element_blank(),
                         axis.line = ggplot2::element_line(colour = "black"))
    return(g1)
}

etree2phylo <- function(etree) {
    ext <- get_extant(tm = etree$ct, tree = etree)
    nw <- newick(ext, CT = etree$ct)
    tr <- ape::read.tree(text = nw)
    return(tr)
}

phylo2etree <- function(phylo) {
    # transformation of ultrametric trees into data frame
    tree <- DDD::phylo2L(phylo)
    brts_dd <- tree[, 1]
    brts <- cumsum(-diff(c(brts_dd, 0)))

    tree <- list(extant = data.frame(brts = c(0, brts[-length(brts)]),
                                     parent = c(1, abs(tree[, 2][-1])),
                                     child = abs(tree[, 3])),
                 extinct = data.frame(brts = numeric(),
                                      parent = numeric(),
                                      child = numeric(),
                                      t_ext = numeric()),
        ct = brts_dd[1])
    tree$extant[3:nrow(tree$extant), "parent"] <-
        tree$extant[3:nrow(tree$extant), "parent"] + 1
    tree$extant$child <- tree$extant$child + 1
    tree$extant$parent[1:2] <- 1
    tree$extinct$parent <- tree$extinct$parent + 1
    tree$extinct$child <- tree$extinct$child + 1
    # NOTE: This next line is a hack
    tree$extant[2, 2] <- 2
    class(tree) <- "etree"
    return(tree)
}

GPD <- function(tree, tm) {
    # input: an ultramedric tree defined by a data.frame with columns brts,
    # parent, child the first two rows
    # are
    n <- nrow(tree)
    child_nms <- as.character(tree$child)
    parent_nms <- as.character(tree$parent)
    species_nms <- child_nms
    gpd <- matrix(0, ncol = n, nrow = n)
    dimnames(gpd) <- list(species_nms, species_nms)
    species <- as.list(1:n)
    for (i in seq(n, 2)) {
        p_set <- species[[which(species_nms == parent_nms[i])]]
        c_set <- species[[which(species_nms == child_nms[i])]]
        gpd[p_set, c_set] <- tm - tree$brts[i]
        species[[which(species_nms == parent_nms[i])]] <- c(p_set, c_set)
    }
    gpd <- gpd + t(gpd)
    return(gpd)
}



# more utilities (emphasis)

n_from_time <- function(tm, tree) {
    # return N at tm.
    if (tm == 0) {
        N <- 2
    } else {
        extended_tree <- extend_tree(tree)

        n <- cumsum(extended_tree$event) + cumsum(extended_tree$event - 1) + 1
        brts <- extended_tree$brts
        if (tm == 0)
            tm <- 1e-15
        N <- n[max(which(brts < tm))]
    }
    return(N)
}

n_for_all_bt <- function(tree) {
    brts <- c(tree$extant$brts, tree$extinct$brts, tree$extinct$t_ext)
    brts <- c(0, sort(brts[brts != 0]))
    n <- sapply(brts, n_from_time, tree)
    return(n)
}

extend_tree <- function(tree) {
    if (is.null(tree$extinct)) {
        tree <- list(extant = tree,
                     extinct = data.frame(brts = NULL,
                                         t_ext = NULL))
    }
    extended_tree <- data.frame(brts = c(tree$extant$brts,
                                         tree$extinct$brts,
                                         tree$extinct$t_ext),
                                event = c(rep(1, nrow(tree$extant)),
                                          rep(1, nrow(tree$extinct)),
                                          rep(0, nrow(tree$extinct))))
    extended_tree <- extended_tree[order(extended_tree$brts), ]
    extended_tree <- rbind(extended_tree,
                           data.frame(brts = tree$ct,
                                     event = 2))
    if (extended_tree$brts[2] == 0) {
        extended_tree <- extended_tree[-1, ]
    }
    return(extended_tree)
}



phylodiversity <- function(tree, tm) {
    # input: an ultrametric tree
    if (!is.null(tree$extant)) {
        tree <- get_extant(tree, tm)
    } else {
        tree <- tree[tree$brts < tm, ]
    }
    dt <- diff(c(tree$brts[-1], tm))
    return(sum(dt * (2:(length(dt) + 2 - 1))))
}

get_current_species <- function(tm, tree) {
    species <- c(tree$extant$child[tree$extant$brts < tm],
                 tree$extinct$child[tree$extinct$t_ext > tm])
    return(species)
}

# time calculation
get.time <- function(time,
                     mode = "sec") {
    dif <- proc.time() - time
    ti <- as.numeric(dif[3])
    if (mode == "min")
        ti <- ti / 60
    if (mode == "hou")
        ti <- ti / 3600
    return(ti)
}

#' function to rescale the branching times of a phylogenetic tree
#' @param phy phy object
#' @param new_crown_age new crown age of the phylogenetic tree, set to 1 by 
#' default
#' @return rescaled phylogeny as a phy object
#' @export
rescale_tree <- function(phy, new_crown_age = 1) {
    ltable <- DDD::phylo2L(phy)
    ltable[, 1] <- new_crown_age * ltable[, 1] / max(ltable[, 1])
    new_phy <- DDD::L2phylo(ltable)
    return(new_phy)
}




# Compute (GPD) for an ultrametric phylogenetic tree
GPD_improved <- function(tree, tm) {
  # Check if the tree data frame has the required columns
  if (!all(c("brts", "parent", "child") %in% names(tree))) {
    stop("The input tree data frame must have columns: 'brts', 'parent', 'child'")
  }
  
  # Initialize variables
  num_rows <- nrow(tree)
  child_names <- as.character(tree$child)
  parent_names <- as.character(tree$parent)
  species_names <- child_names
  gpd_matrix <- matrix(0, ncol = num_rows, nrow = num_rows)
  dimnames(gpd_matrix) <- list(species_names, species_names)
  species_indices <- as.list(1:num_rows)
  
  # Loop through the tree nodes to fill the GPD matrix
  for (i in seq_len(num_rows)) {
    parent_set <- species_indices[[which(species_names == parent_names[i])]]
    child_set <- species_indices[[which(species_names == child_names[i])]]
    gpd_matrix[parent_set, child_set] <- tm - tree$brts[i]
    species_indices[[which(species_names == parent_names[i])]] <- c(parent_set, child_set)
  }
  
  # Make the GPD matrix symmetric
  gpd_matrix <- gpd_matrix + t(gpd_matrix)
  
  return(gpd_matrix)
}


# Compute the phylogenetic distance matrix for an ultrametric tree
computePhyloDistanceMatrix <- function(tree, tm) {
  # Validate that the tree data frame has the required columns
  if (!all(c("brts", "parent", "child") %in% names(tree))) {
    stop("The input tree data frame must have columns: 'brts', 'parent', 'child'")
  }
  
  # Initialize variables
  num_rows <- nrow(tree)
  child_names <- as.character(tree$child)
  parent_names <- as.character(tree$parent)
  species_names <- child_names
  phylo_distance_matrix <- matrix(0, ncol = num_rows, nrow = num_rows)
  dimnames(phylo_distance_matrix) <- list(species_names, species_names)
  species_indices <- as.list(1:num_rows)
  
  # Loop through the tree nodes to fill the phylogenetic distance matrix
  for (i in seq_len(num_rows)) {
    parent_set <- species_indices[[which(species_names == parent_names[i])]]
    child_set <- species_indices[[which(species_names == child_names[i])]]
    phylo_distance_matrix[parent_set, child_set] <- tm - tree$brts[i]
    species_indices[[which(species_names == parent_names[i])]] <- c(parent_set, child_set)
  }
  
  # Make the phylogenetic distance matrix symmetric
  phylo_distance_matrix <- phylo_distance_matrix + t(phylo_distance_matrix)
  
  return(phylo_distance_matrix)
}


## Function from the DDD package 
L2phylo <- function (L, dropextinct = T){
  L = L[order(abs(L[, 3])), 1:4]
  age = L[1, 1]
  L[, 1] = age - L[, 1]
  L[1, 1] = -1
  notmin1 = which(L[, 4] != -1)
  L[notmin1, 4] = age - L[notmin1, 4]
  if (dropextinct == T) {
    sall = which(L[, 4] == -1)
    tend = age
  }
  else {
    sall = which(L[, 4] >= -1)
    tend = (L[, 4] == -1) * age + (L[, 4] > -1) * L[, 4]
  }
  L = L[, -4]
  linlist = cbind(data.frame(L[sall, ]), paste("t", abs(L[sall, 
                                                          3]), sep = ""), tend)
  linlist[, 4] = as.character(linlist[, 4])
  names(linlist) = 1:5
  done = 0
  while (done == 0) {
    j = which.max(linlist[, 1])
    daughter = linlist[j, 3]
    parent = linlist[j, 2]
    parentj = which(parent == linlist[, 3])
    parentinlist = length(parentj)
    if (parentinlist == 1) {
      spec1 = paste(linlist[parentj, 4], ":", linlist[parentj, 
                                                      5] - linlist[j, 1], sep = "")
      spec2 = paste(linlist[j, 4], ":", linlist[j, 5] - 
                      linlist[j, 1], sep = "")
      linlist[parentj, 4] = paste("(", spec1, ",", spec2, 
                                  ")", sep = "")
      linlist[parentj, 5] = linlist[j, 1]
      linlist = linlist[-j, ]
    }
    else {
      linlist[j, 1:3] = L[which(L[, 3] == parent), 1:3]
    }
    if (nrow(linlist) == 1) {
      done = 1
    }
  }
  linlist[4] = paste(linlist[4], ":", linlist[5], ";", sep = "")
  phy = ape::read.tree(text = linlist[1, 4])
  tree = ape::as.phylo(phy)
  return(tree)
}