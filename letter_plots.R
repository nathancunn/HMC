all.lets <- data.frame(key=1:54)
for(i in 1:54){
all.lets$letter[i] <- A[[i]]$letter
}
wordform <- function(x) {
  chars <- strsplit(x,"")
  key <- matrix(0,nrow=length(chars[[1]])) 
  for(i in 1:length(chars[[1]])) {
    if(any(all.lets$key[all.lets$letter==chars[[1]][i]])) {
    key[i,1] <- all.lets$key[all.lets$letter==chars[[1]][i]]
    } else key[i,1] <- 55
  }
  key
}

A[[1]]

plot(0:50,0:50,type="n",asp=1)
j=wordform("O")
for(i in 1:1000) {
    atts <- (length(A[[j]])-1)/3
    lambdas <- cumsum(A[[j]][2:(2+atts-1)])
    group <-  findInterval(runif(1),lambdas)+1
    points(rmvnorm(1,mean = A[[j]][1+atts+group][[1]], sigma = A[[j]][1+2*atts+group][[1]]),cex=0.75)
}

