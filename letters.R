library(raster)
library(jpeg)

img <- readJPEG("../Downloads/Not-giant-enough-letter-w.jpg")

m <- as.raster(img, max = 255)
table(m)
m[m=="#000000"] <- 1
m[m=="#000001"] <- 1
m[m=="#000100"] <- 1
m[m=="#010001"] <- 0
m[m=="#010100"] <- 0
m[m=="#010101"] <- 0


m <- as.matrix(m)
write.csv(m,"m.csv")


m <- matrix(as.numeric(unlist(m)),nrow=nrow(m))

sum(m)
coords <- matrix(0,nrow=sum(m),ncol=2)
count=1
while(count<=sum(m)) {
for(i in 1:nrow(m)) {
  for(j in 1:ncol(m)) {
    if(m[i,j]==1) {
      coords[count,] <- c(550-j,550-i)
      count <- count+1
    }
  }
}
}

outqda <- qda(coords[,1:2], grouping = coords[,3])


coords <- matrix(0,nrow=length(m),ncol=3)
k <- 1
  for(i in 1:nrow(m)) {
    for(j in 1:ncol(m)) {
      coords[k,] <- c(550-j,550-i,m[i,j])
      k <- k+1
      }
    }

library(mixtools)

seqn <- floor(seq(1,nrow(m),length=64))
compressed <- matrix(0,64,64)
for(i in 1:63){
  for(j in 1:63){
 compressed[i,j] <- mean(m[seqn[i]:seqn[i+1],seqn[j]:seqn[j+1]])  
  }
}

compressed <- round(compressed)

coords <- matrix(0,nrow=length(compressed),ncol=3)
k <- 1
for(i in 1:nrow(compressed)) {
  for(j in 1:ncol(compressed)) {
    coords[k,] <- c(550-j,550-i,compressed[i,j])
    k <- k+1
  }
}


out <- mvnormalmixEM(coords[coords[,3]==1,1:2],k=4)


