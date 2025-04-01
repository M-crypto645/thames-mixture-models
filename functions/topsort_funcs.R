# implements the recursive version of the topsort algorithm
alltopsorts_recursion = function(n, adj_list){
  if(length(adj_list) == 0){
    Dlist = 1:n
  } else{
    countmat = sapply(1:n, function(s)
      sapply(adj_list, function(a) a[2]==s))
    if(is.matrix(countmat)){
      counts = colSums(countmat)
    } else{
      counts = countmat + 0
    }
    Dlist = which(counts==0)
  }
  PERMUTATIONRESULTS <<- list(numeric(n))
  alltopsorts(n, adj_list, Dlist, deleted_vertices=c())
  PERMUTATIONRESULTS = PERMUTATIONRESULTS[-length(PERMUTATIONRESULTS)]
  
  for(s in seq_along(PERMUTATIONRESULTS)[-1]){
    PERMUTATIONRESULTS[[s]][PERMUTATIONRESULTS[[s]]==0] = PERMUTATIONRESULTS[[s-1]][PERMUTATIONRESULTS[[s]]==0]
  }
  
  return(PERMUTATIONRESULTS)
}

alltopsorts = function(n, adj_list, Dlist, deleted_vertices){
  base = Dlist[length(Dlist)]
  repeatuntil = FALSE
  while(!repeatuntil){
    #browser()
    q = Dlist[length(Dlist)]
    deleted_vertices_new = c(deleted_vertices,q)
    if(length(adj_list)==0){
      counts = rep(0,n)
      adj_list_new = adj_list
    } else{
      adj_list_new = adj_list[!sapply(adj_list,function(a) a[1]==q)]
      # browser()
      if(length(adj_list_new)==0){
        counts = rep(0,n)
      } else{
        countmat = sapply(1:n, function(s)
          sapply(adj_list_new, function(a) a[2]==s))
        #if(length(countmat[[1]])!=0){
        #  countmat = countmat + 0
        #}
        if(is.matrix(countmat)){
          counts = colSums(countmat)
        } else{
          counts = countmat + 0
        }
      }
    }
    Dlist_new = ((1:n)*(counts==0))[-deleted_vertices_new]
    Dlist_new = Dlist_new[Dlist_new!=0]
    #browser()
    #print("before")
    rescopy = PERMUTATIONRESULTS
    rescopy[[length(rescopy)]][length(deleted_vertices_new)] <- q 
    PERMUTATIONRESULTS <<- rescopy
    #print("after")
    #print(paste0("output ",q," in column ",length(deleted_vertices_new)))
    if(length(deleted_vertices_new) == n-1){
      rescopy = PERMUTATIONRESULTS
      rescopy[[length(rescopy)]][length(deleted_vertices_new)] <- q 
      PERMUTATIONRESULTS <<- rescopy
      #print(paste0("output ",(1:n)[-deleted_vertices_new]," in column ",n))
      #print("newline")
      #browser()
      recopy = PERMUTATIONRESULTS
      rescopy[[length(rescopy)]][n] <- (1:n)[-deleted_vertices_new]
      rescopy[[length(rescopy)+1]] = numeric(n)
      PERMUTATIONRESULTS <<- rescopy
    } else{
      alltopsorts(n, adj_list_new, Dlist_new, deleted_vertices_new)
    }
    Dlist = c(q,Dlist[-length(Dlist)])
    #browser()
    repeatuntil = (base == Dlist[length(Dlist)])
  }
}

#n = 4
#res = alltopsorts_recursion(n, list(c(1,3),c(2,4)))
