# modified from CliquePercolation package

cpAlgorithm2 <- function(W, k){
  
  ## W is matrix of undirected network
  ###check whether W is a qgraph object
  ###if not check whether matrix is symmetric and convert to qgraph object
  ###error message if k is not larger than 2
  if (k < 3) {
    stop("k must be larger than 2, because this is the first reasonable clique size.")
  }

  
  ###function to determine list of cliques for unweighted networks
  unweighted <- function(W_unweighted){
    
    #transform to igraph object to use specific functions
    W_i_unweighted <- igraph::graph_from_adjacency_matrix(W_unweighted,
                                                          mode = "undirected", weighted = NULL)
    
    #extract cliques
    cliques_unweighted <- igraph::cliques(W_i_unweighted, min = k, max = k) %>% lapply(as.vector)
    
    #return cliques
    return(cliques_unweighted)
  }
  
  
    

  
  #function to derive communities, shared nodes, and isolated nodes
  results_cp <- function(W, cliques, labels){
    
    #loop to compare all cliques with each other
    #if cliques share k-1 nodes, a vector is created stating their indices
    #if there are not at least two cliques, there can be no edge; thus, the edge list is empty
    if (length(cliques) > 1) {
      edges <- list()
      for (i in 1:(length(cliques) - 1)) {
        for (j in (i + 1):length(cliques)) {
          if (length(Reduce(intersect, list(cliques[[i]],cliques[[j]]))) == (k - 1)) {
            edges[[length(edges) + 1]] <- c(i,j)
          }
        }
      }
    } else {edges <- list()}

    #if there is more than one clique and at least one edge...
    #create list of vectors with each vector being a community; this has a number of steps
    #a weights matrix is created with cliques as nodes and existence of k-adjacency between them as edge
    #then, components (connected subgraphs) of this graph are extracted from igraph object
    #components are therefore the communities
    #each clique is then put into its respective community and the nodes of each community are extracted
    if (length(cliques) > 1 & length(edges) > 0) {
      W_comm <- matrix(c(0), nrow = length(cliques), ncol = length(cliques),
                       byrow = TRUE)
      for (i in 1:length(edges)) {
        W_comm[edges[[i]][1],edges[[i]][2]] <- 1
      }
      W_i_comm <- igraph::graph_from_adjacency_matrix(W_comm, mode = "upper")
      members <- igraph::components(W_i_comm)$membership
      split <- split(cliques, members)
      communities <- lapply(split, function(x) sort(unique(unlist(x))))
    }
    #if there is at least one clique but no edge...
    #communities are the existing cliques
    if (length(cliques) > 0 & length(edges) == 0) {
      communities <- cliques
    }
    #if there are no cliques (and therefore also no edges)...
    #communities list is empty
    if (length(cliques) == 0) {
      communities <- list() 
    }
    
    #create vector of nodes that do not belong to a community
    #if there are no isolated nodes, create empty variable
    isolated <- subset(1:nrow(W),
                       subset = !(1:nrow(W)%in%unique(unlist(communities))))
    if (length(isolated) == 0) {isolated <-  c()}
    
    #if there is more than one community...
    #get shared nodes
    #if there are no shared nodes, create empty variable
    if (length(communities) > 1) {
      shared <- list()
      count <- 1
      for (i in 1:(length(communities) - 1)) {
        for (j in (i + 1):length(communities)) {
          shared[[count]] <- Reduce(intersect, list(communities[[i]],communities[[j]]))
          count <- count + 1
        }
      }
      shared <- sort(unique(unlist(shared)))
      if (length(shared) == 0) {shared <- c()}
    }
    
    #if there are communities...
    #create list of communities with labels from qgraph object instead of node numbers
    if (length(communities) > 0) {
      communities_labels <- list()
      for (i in 1:length(communities)) {
        communities_labels[[i]] <- labels[communities[[i]]]
      }
    }
    #if there are no communities...
    if (length(communities) == 0) {
      communities_labels <- list()
    }
    
    #create vector of isolated nodes with labels from qgraph object instead of node numbers
    #if there are no isolated nodes, create empty variable
    isolated_labels <- labels[isolated]
    if (length(isolated_labels) == 0) {isolated_labels <- c()}
    
    #if there is more than one community...
    #create vector with shared nodes with labels from qgraph object instead of node numbers
    if (length(communities) > 1) {
      shared_labels <- labels[shared]
    }
    #if there is less than two communities
    if (length(communities) < 2) {
      shared <- c()
      shared_labels <- c()
    }
    #if there are no shared nodes, create empty variable
    if (length(shared_labels) == 0) {shared_labels <- c()}
    
    #return community lists, shared nodes vectors, and isolated nodes vectors
    return(list(communities,communities_labels,
                shared,shared_labels,
                isolated,isolated_labels))
  }
  
  #extract weights matrix
  #Wmat <- qgraph::getWmat(W)
  Wmat <- W
  labels <- as.vector(rownames(Wmat))
  
  #run corresponding functions for respective method
  
  cliques_unweighted <- unweighted(Wmat)
  results_unweighted <- results_cp(Wmat, cliques_unweighted, labels)
  names(results_unweighted) <- c("list.of.communities.numbers","list.of.communities.labels",
                                   "shared.nodes.numbers","shared.nodes.labels",
                                   "isolated.nodes.numbers","isolated.nodes.labels")
  returned_object <- c(results_unweighted,list(k = k))

  
  class(returned_object) <- "cpAlgorithm2"
  
  return(returned_object)
  
}

#cpAlgorithm2(W, 3)
