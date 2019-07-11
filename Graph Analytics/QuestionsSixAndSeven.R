library(igraph)
library(datastructures)

find_longest_path <- function(graph){
     maxDistance <- 0
     for(i in V(graph)){
         for(j in V(graph)){
             if(i == j){ next }
             distance <- max(distances(graph, v = i, to = j))
             if(distance > maxDistance){
                 maxDistance <- distance
             }
         }
     }
     return(maxDistance)
}

contains_multiple_longest_paths <- function(graph, longestPath){
    for( i in V(graph)){
      for(j in V(graph)){
        if(i ==j){ next }
        if(max(distances(graph, v = i, to = j)) == longestPath){ return(TRUE)}
      }
    }
    return(FALSE)
}

determine_neighborhood_size <- function(graph, vertex){
    #Level one neighbors
    neighborhoodMembers <- vector()
    #Add the vertex to the neighborhood
    neighborhoodMembers <- c(neighborhoodMembers, names(vertex))
    levelOneNeighbors <- neighbors(graph, vertex)
    l1Reachability <- vector()
    l2Reachability <- vector()
    l3Reachability <- vector()
    #Loop through all neighbors of vertex
    for(levelOneNeighbor in levelOneNeighbors){
        neighborhoodMembers <- c(neighborhoodMembers, levelOneNeighbor)
        l1Reachability <- c(l1Reachability, levelOneNeighbor)
        #Loop through all Level 2 Neighbors
        for(levelTwoNeighbor in neighbors(graph, levelOneNeighbor)){
            if(! is.element(levelTwoNeighbor, neighborhoodMembers)){
                   neighborhoodMembers <- c(neighborhoodMembers, levelTwoNeighbor)
            }
            if(! is.element(levelTwoNeighbor, l2Reachability)){
              l2Reachability <- c(l2Reachability, levelTwoNeighbor)
            }
            #Loop through all Level 3 Neighbors
            for(levelThreeNeighbor in neighbors(graph, levelTwoNeighbor)){
                if(! is.element(levelThreeNeighbor, l3Reachability)){
                  l3Reachability <- c(l3Reachability, levelThreeNeighbor)
                }
                if(! is.element(levelThreeNeighbor, neighborhoodMembers)){
                     neighborhoodMembers <- c(neighborhoodMembers, levelThreeNeighbor)
                }
            }
        }
    }
    
    vertexNeighborhood <- induced.subgraph(graph, neighborhoodMembers)
    diam <- diameter(vertexNeighborhood)
    #Return an object containing Level One Neighbors, Level Two Neighbors,
    # Level Three Neighbors, the neighborhood, and the size of the neighborhood
    return(list(l1Reachability=l1Reachability,
                l2Reachability=l2Reachability,
                l3Reachability=l3Reachability,
                vertexNeighborhood=V(vertexNeighborhood),
                neighborhoodSize=diam,
                vertex=vertex))
}

pretty_print_neighbors <- function(neighorList){
  neighborString <- paste(neighorList)
  if(grepl("c", neighborString)){
    return(substr(neighborString,3,nchar(neighborString)-1))
  } else{
    return(neighborString)
  }
}

build_neighborhood_matrix <- function(neighborhoods){
    #Retrieve the neighborhood sizes for each vertex
    # The fibonacci heap from the datastructures package only offers min-heaps,
    # so we subtract from a large number to get the key.
    maximumInteger <- .Machine$integer.max
    neighbohoodSizes <- vector()
    for(n in neighborhoods){
      neighbohoodSizes <- c(neighbohoodSizes, maximumInteger - n[[5]])
    }
    #Construct a fibonacci heap, so that we can get the twenty largest nodes by neighborhood
    # size in O(log n). Note that it's O(log n) because we're only calling delete-min 20 times.
    fheap <- fibonacci_heap("numeric")
    fheap <- insert(fheap, neighbohoodSizes, neighborhoods)
    
    nodeIDs = vector()
    neighbors <- vector()
  
    #Build a matrix to store our nodes 
    neighborhoodMatrix <- matrix(nrow=20,ncol=3)
    for(i in 1:20){
      largeNeighborhood <- pop(fheap)
      neighborhoodMatrix[i,1] <- pretty_print_neighbors(largeNeighborhood[[1]][1])
      neighborhoodMatrix[i,2] <- pretty_print_neighbors(largeNeighborhood[[1]][2])
      neighborhoodMatrix[i,3] <- pretty_print_neighbors(largeNeighborhood[[1]][3])
      
      #Store neighborhood information to find common nodes later
      verticesOfLargeNeighborhood <- largeNeighborhood[[1]][4]
      
      neighbors <-  c(neighbors,  verticesOfLargeNeighborhood)
      nodeIDs <- c(nodeIDs, largeNeighborhood[[1]][6])
    }
    rownames(neighborhoodMatrix) <- nodeIDs
    
    #Find common nodes...R doesn't have a built-in hash mechanism, so this is painful
    printedNodes <- vector()
    for(j in 1:20){
      #Get a list of vectors for 
      neighborList <- names(neighbors[[j]])
      #Loop through every node in neighborhood j
      for(neighbor in neighborList){
        commonNodes <- c(nodeIDs[[j]])
        #Loop through all other neighborhoods
        for(k in 1:20){
          if(j == k){ next;}
          if(is.element(neighbor, names(neighbors[[k]]))){
            if(! is.element(nodeIDs[[k]], commonNodes)){
              commonNodes <- c(commonNodes, nodeIDs[[k]])
            }
          }
        }
        if(length(commonNodes) > 1){
          if(! is.element(neighbor, printedNodes)){
            print(sprintf("Node %1s is shared by %2s", neighbor, paste(commonNodes)))
            printedNodes <- c(printedNodes, neighbor)
          }
        }
      }
    }
    
    #Return matrix of 20 nodes with their reachability to the 3rd level
    return(neighborhoodMatrix)
}


#Set working directory
setwd('/Users/kevintyler/Documents/Programming/csci6444')

#Read in table
print("Reading table")
edges <- read.table('roadNet-CA.txt')

#Convert to matrix
edgeMatrix <- as.matrix(edges)

#Extract vectors from matrix
v1 <- edgeMatrix[,1]
v2 <- edgeMatrix[,2]

print("Building relations from table")
#Build relations
relations <- data.frame(from=v1,to=v2)

#Construct graph
print("Constructing graph from relations")
roadGraph <- graph.data.frame(relations,directed=TRUE)

#Remove nodes of degree < 6 
#TODO: Replace variable mentions below
#srg6 <- induced.subgraph(simplifiedRoadGraph, which(degree(simplifiedRoadGraph) >= 6))
print("Calculating mean degree")
roadGraphDegree <- degree(roadGraph)
meanDegree <- ceiling(mean(roadGraphDegree))
print("Pruning nodes from graph")
print(sprintf("Number of nodes before prune: %s", vcount(roadGraph)))
roadGraph <- delete_vertices(roadGraph, which(roadGraphDegree <= 8))
print(sprintf("Number of nodes after prune: %s", vcount(roadGraph)))

neighborhoods <- list()
counter <- 1
for(vertex in V(roadGraph)){
  neighborhoods[[counter]] <- determine_neighborhood_size(roadGraph,vertex)
  counter <- counter + 1
}
  
neighborhoodMatrix <- build_neighborhood_matrix(neighborhoods)
print(neighborhoodMatrix)


  
#We want the "central person" in the graph, so we'll find the most central nodes
# with respect to alpha-centrality.

#Alpha-centrality
print(sprintf("Most important node (i.e.: Alpha-centrality) of road graph: %s", 
              names(which.max(alpha_centrality(roadGraph, alpha=0.9)))))

#Find the largest clique of subgraph i
print(sprintf("Largest clique of road graph %s", largest_cliques(roadGraph)[[1]]))

#Find the ego of subgraph i
egoSize <- ego_size(roadGraph)
egoMaxIndex <- which.max(egoSize)
egoNodeSize <- egoSize[egoMaxIndex]
egoNodeLabel <- names(roadGraph[[egoMaxIndex]])
multipleEgosOfMaximumSize <- (length(which(egoSize == egoNodeSize)) > 1)
print(sprintf("Ego of road graph: Node %1s, which is of size %2s", 
              egoNodeLabel, egoNodeSize))
print(sprintf("Are there multiple node with the highest ego? %s", multipleEgosOfMaximumSize))

#Find betweenness centrality 
print(sprintf("Betweeness centrality of road graph: %s", 
              names(which.max(betweenness(roadGraph)))))

#Find the power centrality
print(sprintf("Power centrality of road graph: %s", 
              names(which.max(power_centrality(roadGraph, exponent = .7)))))

# Is there more than one person with the most degrees?
print(sprintf("Are there >1 nodes of maximum degree for road graph: %s", 
              length(which(degree(roadGraph) == max(degree(roadGraph)))) > 1))

#Are there multiple cliques?
print(sprintf("Are their mulitiple cliques for subgraph %s",  
              clique_num(roadGraph) > 1))
  
reducedRoadGraph <- graph.data.frame(relations,directed=TRUE)
reducedRoadGraph <- delete_vertices(roadGraph, which(degree(reducedRoadGraph) <= 12))

#Find the longest path 
longestPath <- find_longest_path(reducedRoadGraph)
print(sprintf("Longest path of reduced road graph: %s", longestPath))

#Are there multiple longest paths?
print(sprintf("Are their mulitiple longest paths for reduced road graph %s", 
              contains_multiple_longest_paths(reducedRoadGraph, longestPath)))