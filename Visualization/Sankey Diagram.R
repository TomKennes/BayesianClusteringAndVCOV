# Sankey Diagram

install.packages("networkD3")
install.packages("matrixcalc")
library("networkD3")
library("matrixcalc")

install.packages("riverplot")
library("riverplot")
set.seed(9999)
simulate.clusterdata <- function(N = 100, crises = c(), size.crises = length(crises), groups = 6, K = 20, S = 0.1, cluster = F){
  base = c()
  for(i in 1:groups){
    if(i == 1){
      base = rnorm(N, i, S)
    } else {
      base = cbind(base, rnorm(N, i, S))
    }
  }
  x = c()
  if(!cluster){
    for(i in 1:K){
      if(i == 1){
        x = base[,1] + rnorm(N, 0, S/10)
      } else {
        x = cbind(x, base[,(i %% groups+ 1)] + rnorm(N, 0, S/100))
      }
    }
  } else {
    for(i in 1:K){
      if(i == 1){
        x = base[,1] 
      } else {
        x = cbind(x, base[,(i %% groups + 1)])
      }
    }
  }
  
  if(length(crises) > 0){
    for(i in 1:length(crises)){
      if(crises[i] > size.crises){
        x[(crises[i] - size.crises):(crises[i] + size.crises), ] = rnorm(ncol(x)*(2*size.crises + 1), mean(x), 0.01)
      }
    }
  }
  
  if(length(crises) == 2){
    x[(crises[1] + 2):(crises[2] - 2),] = matrix(sample(x[(crises[1] + 2):(crises[2] - 2),]), nrow = nrow(x[(crises[1] + 2):(crises[2] - 2),]))
  }
    
  return(x)
}


# I need to start with zeros so I cannot work around that in this setting.
sankey <- function(clusters, notepadding = 50, notewidth = 10, fontsize = 10, var_names = as.character(seq(1, ncol(clusters), by = 1))){
  # Set up from-to-node
  counter = 1
  from_to_node <- c(0,0,0,1,0,0)
  for(i in 1:ncol(clusters)){
    from_to_node <- rbind(from_to_node, c(0,clusters[1,i], counter, 1,i,1))
    counter = counter + 1
  }
  colnames(from_to_node) <- c("from", "to", "node", "weight", "t", "series")
  
  for(i in 2:nrow(clusters)){
    for(j in 1:ncol(clusters)){
      from_to_node <- rbind(from_to_node, c(clusters[i-1,j], clusters[i,j],counter,1,j,i))
      counter = counter + 1
    }
  }
  from_to_node <- from_to_node[2:nrow(from_to_node),]
  
  #from_to_node
  # Make corresponding adjustments
  # Combine nodes 
  n = 1
  while(n < nrow(from_to_node)){
    for(i in n:(n + ncol(clusters) - 2)){
      for(j in (i+1):(n + ncol(clusters) - 1)){
        #print(paste("i: ", i, " j: ", j))
        if(from_to_node[i,2] == from_to_node[j,2]){
          from_to_node[j,3] = from_to_node[i,3]
          if(from_to_node[i,1] == from_to_node[j,1]){
            from_to_node[i,4] <- from_to_node[i, 4] + 1
            from_to_node[j,4] <- from_to_node[j, 4] + 1
          }
        }
      }
    }
    n = n + ncol(clusters)
  }
  
  from_to_node <- from_to_node[,c(1,3:ncol(from_to_node))]
  #from_to_node
  
  for(i in (ncol(clusters)+1):nrow(from_to_node)){
    from_to_node[i,1] <- from_to_node[i - ncol(clusters),2]
  }
  #from_to_node
  # Scale down the to-nodes
  current = 1
  counter = 1
  for(i in seq(1,nrow(from_to_node),by = ncol(clusters))){
    subset <- from_to_node[i:(i+ncol(clusters) - 1),2]
    m <- max(subset)
    while(current <= m){
      if(sum(subset == current) == 0){
        current = current+ 1
      } else {
        subset[subset == current] <- counter
        counter = counter + 1
        current = current + 1
      }
    }
    from_to_node[i:(i+ncol(clusters) - 1),2] <- subset
  }
  #from_to_node
  
  
  for(i in 1:nrow(from_to_node)){
    from_to_node[i,3] = (from_to_node[,1] == from_to_node[i,1]) %*% (from_to_node[,2] == from_to_node[i,2])
  }
  #from_to_node
  
  for(i in (ncol(clusters) + 1):nrow(from_to_node)){
    from_to_node[i,1] <- from_to_node[i - ncol(clusters),2]
  }
  #from_to_node
  if(ncol(clusters) > 10){
    from_to_node[from_to_node[,1] == 0, 3] <- from_to_node[from_to_node[,1] == 0, 3]/100
  }
  
  p = max(from_to_node[,2])
  
  nodes = data.frame("name" = 
                       rep("", (p + 1)))
  group <- as.character(from_to_node[,4])
  tmp <- from_to_node[,1:3]
  tmp <- vec(t(tmp))
  links = as.data.frame(matrix(tmp,# The third number is the value of the node
    byrow = TRUE, ncol = 3))
  names(links) = c("source", "target", "value")
  occs <- rep("", nrow(links))
  
  for(i in 1:nrow(links)){
    occs[i] = var_names[from_to_node[i,4]]
  }
  links = cbind(links, occs)
  names(links) = c("source", "target", "value", "group")
  
  
    return(sankeyNetwork(Links = links, 
                       Nodes = nodes,
                       Source = "source", 
                       Target = "target",
                       Value = "value", 
                       #NodeID = "name",
                       LinkGroup = "group",
                       fontSize= fontsize, 
                       nodeWidth = notewidth, 
                       nodePadding = notepadding,
                       #,NodeGroup = "group",
                       sinksRight = T
                       ))
  }





# Creating Simulated Data
x = simulate.clusterdata(cluster = F, S = 0.1, crises = c(15,35), N = 50, K = 10, groups = 3)
colnames(x) = paste("x", seq(1,ncol(x), by = 1), sep = "")
write.csv(x, "simulated data clustering.csv")

x = round(x)
write.csv(x, "simulated data clustering exact.csv")
p = sankey(x, notewidth = 2, notepadding = 5)
p


# Analysing clustered results
setwd("C:/Users/Tom Kennes/Desktop/Research Master/Thesis/Core/Functions and Codes")
clusters <- read.table("results for clustering.txt")
clusters <- clusters[8:nrow(clusters),]
for(i in 1:((ncol(clusters) - 1)/4)){
  p <- sankey(clusters[(2 + 4*(i - 1)):(5 + 4*(i - 1)),], notepadding = 50, notewidth = 1)
  print(p)
}




