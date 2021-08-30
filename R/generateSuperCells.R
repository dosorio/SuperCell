generateSuperCells <- function(X){
  N.c <- ncol(X)
  n.var.genes = min(1000, nrow(X))
  gamma = 10
  k.knn = 5
  do.scale = TRUE
  n.pc = 10
  fast.pca = TRUE
  do.approx = FALSE
  approx.N = 20000
  use.nn2 = TRUE
  seed = 12345
  igraph.clustering = c("walktrap", "louvain")
  return.singlecell.NW = TRUE
  return.hierarchical.structure = TRUE
  block.size = 10000

  n.var.genes <- min(n.var.genes, nrow(X))
  if(N.c > 50000){
    set.seed(seed)
    idx         <- sample(N.c, 50000)
    gene.var    <- apply(X[,idx], 1, var)
  } else {
    gene.var    <- apply(X, 1, var)
  }
  
  genes.use   <- names(sort(gene.var, decreasing = TRUE))[1:n.var.genes]

  if(do.approx & approx.N >= N.c){
    do.approx <- FALSE
    warning("N.approx is larger or equal to the number of single cells, thus, an exact simplification will be performed")
  }
  
  if(do.approx & (approx.N < round(N.c/gamma))){
    approx.N <- round(N.c/gamma)
    warning(paste("N.approx is set to N.SC", approx.N))
  }
  
  if(do.approx & ((N.c/gamma) > (approx.N/3))){
    warning("N.approx is not much larger than desired number of super-cells, so an approximate simplification may take londer than an exact one!")
  }
  
  if(do.approx){
    set.seed(seed)
    approx.N            <- min(approx.N, N.c)
    presample           <- sample(1:N.c, size = approx.N, replace = FALSE)
    presampled.cell.ids <- colnames(X)[sort(presample)]
    rest.cell.ids       <- setdiff(colnames(X), presampled.cell.ids)
  } else {
    presampled.cell.ids <- colnames(X)
    rest.cell.ids       <- c()
  }
  
  X.for.pca            <- Matrix::t(X[genes.use, presampled.cell.ids])
  if(do.scale){ X.for.pca            <- scale(X.for.pca) }
  X.for.pca[is.na(X.for.pca)] <- 0
  
  if(is.null(n.pc[1]) | min(n.pc) < 1){stop("Please, provide a range or a number of components to use: n.pc")}
  if(length(n.pc)==1) n.pc <- 1:n.pc
  
  if(fast.pca & (N.c < 1000)){
    warning("Normal pca is computed because number of cell is low for irlba::irlba()")
    fast.pca <- FALSE
  }
  
  if(!fast.pca){
    PCA.presampled          <- prcomp(X.for.pca, rank. = max(n.pc), scale. = F, center = F)
  } else {
    PCA.presampled          <- irlba::irlba(X.for.pca, nv = max(n.pc, 25))
    PCA.presampled$x        <- PCA.presampled$u %*% diag(PCA.presampled$d)
    PCA.presampled$rotation <- PCA.presampled$v
  }
  
  build_knn_graph <- function(X, k = 5, from = c("dist", "coordinates"), use.nn2 = TRUE, return_neighbors_order = F, dist_method = "euclidean", cor_method = "pearson", p = 2){
    av.methods <- c("dist", "coordinates")
    method <-  pmatch(from[1], av.methods)
    if(is.na(method)){
      stop(paste("Unlnown method", from, "Available methods are", paste(av.methods, collapse = ", ")))
    }
    
    
    if (method == 2){ # from coordinates
      if(use.nn2){
        if(dist_method != "euclidean"){
          stop(paste0("Fast nn2 function from RANN package is used, so",
                      dist_method, "distnce is not acceptable.
                  To use nn2 method, please, choose eucleadian distance.
                  If you want to use", dist_method, "distance, please set parameter use.nn2 to FALSE"))}
        return(build_knn_graph_nn2(X = X, k = k))
      } else {
        av.dist      <- c("cor", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
        dist_method_ <-  pmatch(dist_method, av.dist)
        if(is.na(dist_method_)){
          stop(paste("Unknown distance method:", dist_method, "Available dist methods are", paste(av.dist, collapse = ", ") ))
        }
        if(dist_method_ == 1){
          #print("cor")
          #print(cor_method)
          av.cor_methods <- c("pearson", "kendall", "spearman")
          cor_method_    <- pmatch(cor_method, av.cor_methods)
          if(is.na(cor_method_)){
            stop(paste("Unknown cor method:", cor_method, "Available cor methods are", paste(av.cor_methods)))
          }
          X <- as.dist(as.matrix(1 - cor(t(X), method = cor_method)))
        } else {
          X <- dist(X, method = dist_method)
        }
      }
    } else {
      if(use.nn2 == TRUE){
        stop("Method nn2 cannot be applied to distance, to use fast nn2 method, please provide coordinates rather than distance
             and set parameter from to coordinates")
      }
      return(knn_graph_from_dist(D = X, k = k, return_neighbors_order = return_neighbors_order))
    }
    
    
    ### now X is distance in any case
    return(knn_graph_from_dist(D = X, k = k, return_neighbors_order = return_neighbors_order))
    
  }
  
  
  
  build_knn_graph_nn2 <- function(X, k = min(5, ncol(X))){
    nn2.res <- RANN::nn2(data = X, k = k)
    nn2.res <- nn2.res$nn.idx
    
    adj.knn       <- split(nn2.res, rep(1:nrow(nn2.res), times = ncol(nn2.res))) # get adj list
    
    graph.knn     <- igraph::graph_from_adj_list(adj.knn,  duplicate = F, mode = "all")
    
    graph.knn     <- igraph::simplify(graph.knn, remove.multiple = T)
    igraph::E(graph.knn)$weight <- 1
    
    return(res <- list(graph.knn = graph.knn))
    
  }
  
  
  knn_graph_from_dist <- function(D, k = 5, return_neighbors_order = T){
    
    ##print("Start knn_graph_from_dist")
    if(!is.matrix(D) & class(D) != "dist"){
      stop("D (X) mast be a matrix or dist!")
    }
    
    if(class(D) != "dist"){
      D <- as.dist(D)
    }
    
    
    
    N        <- (1 + sqrt(1+8*length(D)))/2 # number of cells
    
    if (k >= N)
      stop("Not enought neighbors in data set!")
    if (k < 1)
      stop("Invalid number of nearest neighbors, k must be >= 1!")
    
    row <- function(i, N){
      return(c(if(i>1) D[(i-1)+c(0:(i-2))*(N - 1 - c(1:(i-1))/2)],
               NA,
               if(i < N) D[((i-1)*(N-1) - ((i-1)*(i-2)/2) + 1) : (((i-1)*(N-1) - ((i-1)*(i-2)/2) + 1) + N-i-1)]))
    }
    
    neighbors <- t(sapply(1:N, function(i) {order(row(i,N))[1:k]}))
    
    adj.knn <- split(neighbors, rep(1:nrow(neighbors), times = ncol(neighbors)))
    
    
    graph.knn     <- igraph::graph_from_adj_list(adj.knn,  duplicate = F, mode = "all")
    graph.knn     <- igraph::simplify(graph.knn, remove.multiple = T)
    igraph::E(graph.knn)$weight <- 1
    
    if(return_neighbors_order){
      res <- list(graph.knn = graph.knn,
                  order = neighbors)
    } else {res <- list(graph.knn = graph.knn)}
    
    return(res)
  }
  
  sc.nw <- build_knn_graph(X = PCA.presampled$x[,n.pc], k = k.knn, from = "coordinates", use.nn2 = use.nn2, dist_method = "euclidean")
  
  #simplify
  
  k   <- round(N.c/gamma)
  
  if(igraph.clustering[1] == "walktrap"){
    g.s              <- igraph::cluster_walktrap(sc.nw$graph.knn)
    g.s$membership   <- igraph::cut_at(g.s, k)
    
  } else if(igraph.clustering[1] == "louvain") {
    warning(paste("igraph.clustering =", igraph.clustering, ", gamma is ignored"))
    g.s    <- igraph::cluster_louvain(sc.nw$graph.knn)
    
  } else {
    stop(paste("Unknown clustering method (", igraph.clustering, "), please use louvain or walkrtap"))
  }
  
  membership.presampled        <- g.s$membership
  names(membership.presampled) <- presampled.cell.ids
  
  SC.NW                        <- igraph::contract(sc.nw$graph.knn, membership.presampled)
  SC.NW                        <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")
  
  if(do.approx){
    
    PCA.averaged.SC      <- as.matrix(Matrix::t(supercell_GE(t(PCA.presampled$x[,n.pc]), groups = membership.presampled)))
    X.for.roration       <- Matrix::t(X[genes.use, rest.cell.ids])
    
    
    
    if(do.scale){ X.for.roration <- scale(X.for.roration) }
    X.for.roration[is.na(X.for.roration)] <- 0
    
    
    membership.omitted   <- c()
    if(is.null(block.size) | is.na(block.size)) block.size <- 10000
    
    N.blocks <- length(rest.cell.ids)%/%block.size
    if(length(rest.cell.ids)%%block.size > 0) N.blocks <- N.blocks+1
    
    
    if(N.blocks>0){
      for(i in 1:N.blocks){ # compute knn by blocks
        idx.begin <- (i-1)*block.size + 1
        idx.end   <- min(i*block.size,  length(rest.cell.ids))
        
        cur.rest.cell.ids    <- rest.cell.ids[idx.begin:idx.end]
        
        PCA.ommited          <- X.for.roration[cur.rest.cell.ids,] %*% PCA.presampled$rotation[, n.pc] ###
        
        D.omitted.subsampled <- proxy::dist(PCA.ommited, PCA.averaged.SC) ###
        
        membership.omitted.cur        <- apply(D.omitted.subsampled, 1, which.min) ###
        names(membership.omitted.cur) <- cur.rest.cell.ids ###
        
        membership.omitted   <- c(membership.omitted, membership.omitted.cur)
      }
    }
    
    membership.all       <- c(membership.presampled, membership.omitted)
    membership.all       <- membership.all[colnames(X)]
  } else {
    membership.all       <- membership.presampled[colnames(X)]
  }
  
  membership       <- membership.all
  groups <- membership
  
  weights = NULL
  do.median.norm = FALSE
  
  if(ncol(X) != length(groups)){
    stop("Length of the vector groups has to be equal to the number of cols in matrix ge")
  }
  N.SC <- max(groups)
  supercell_size <- as.vector(table(groups))
  j <- rep(1:N.SC, supercell_size) # column indices of matrix M.AV that, whene GE.SC <- ge %M.AV%
  
  goups.idx  <- plyr:::split_indices(groups)
  i <- unlist(goups.idx) # row indices of matrix M.AV that, whene GE.SC <- ge %M.AV%
  
  if(is.null(weights)){
    M.AV <- Matrix::sparseMatrix(i = i, j = j)
    GE.SC <- X %*% M.AV
    GE.SC <- sweep(GE.SC, 2, supercell_size, "/")
  } 
  return(GE.SC)
}
