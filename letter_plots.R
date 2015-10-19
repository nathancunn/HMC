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
for(i in 1:26) {
  atts <- (length(A[[1+2*i]])-1)/3
  for(j in 1:atts){
    A[[1+(i*2)]][1+atts+j][[1]] <- A[[1+(i*2)]][1+atts+j][[1]]+c(0,+2)
  }
}
  
WordPrint <- function(word,cols=NULL) {
  letters <- wordform(word)
  n <- length(letters)
  if(is.null(cols)) {cols=rep(1,n)}
  eq_probs <- seq(1/n,1,length=n)
  letter <- findInterval(runif(1),eq_probs)+1
  j <- letters[letter]
  displacement <- (0:(n-1))*18-20
  atts <- (length(A[[j]])-1)/3
  lambdas <- cumsum(A[[j]][2:(2+atts-1)])
  group <-  findInterval(runif(1),lambdas)+1
  mean.group <- A[[j]][1+atts+group][[1]]+c(displacement[letter],0)
  var.group <- A[[j]][1+2*atts+group][[1]]
  points(rmvnorm(1,mu = mean.group, sigma = var.group),cex=0.75,col=cols[letter])
}
plot(0:50,0:50,type="n",asp=1,bty="n",xaxt="n",yaxt="n",ann=F)
for(i in 1:10000) {
  WordPrint("HMC",cols=c(rep("#014FDC",4),rep("#01A401",2)))
}

