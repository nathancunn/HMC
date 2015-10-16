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

maxdim <- 32
seqn <- floor(seq(1,nrow(m),length=maxdim))
compressed <- matrix(0,maxdim,maxdim)
for(i in 1:(maxdim-1)){
  for(j in 1:(maxdim-1)){
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
plot(coords[,1],coords[,2],col=coords[,3],asp=1)

out <- mvnormalmixEM(coords[coords[,3]==1,1:2],k=2)
#Creating each of the letters as images
setwd("Letters/")
all.letters <- c(letters,toupper(letters),"?","!")
for(i in 1:length(all.letters)){
jpeg(filename = paste(all.letters[i],".jpg",sep=""))
plot(0,0,pch=all.letters[i],xaxt="n",yaxt="n",ann=F, bty="n",cex=30)
dev.off()
}


out <- list()
files <- list.files()
maxdim <- 32 #allows for a 32*32 grid representing the letters
seqn <- floor(seq(1,nrow(m),length=maxdim))
ks <-c(2,4,5,3,4,6,4,4,4,4,4,4,3,3,5,5,3,3,2,2,3,3,3,3,2,2,5,4,3,3,4,6,4,3,4,7,
       2,5,5,5,3,2,3,3,2,2,4,4,2,2,3,3,3,3)
for(i in 32:length(files)) {
  # Read in the image and convert it into a matrix
  img <- readJPEG(files[i])
  m <- as.raster(img, max = 255)
  m[m=="#000000"] <- 1
  m[m=="#010101"] <- 0
  # Convert the matrix to numeric
  m <- matrix(as.numeric(unlist(m)),nrow=nrow(m))
  # Compressing the matrix to save on computation time
  compressed <- matrix(0,maxdim,maxdim)
  for(j in 1:(maxdim-1)){
    for(k in 1:(maxdim-1)){
      compressed[j,k] <- mean(m[seqn[j]:seqn[j+1],seqn[k]:seqn[k+1]])  
    }
  }
  # Rounding so we only have 1's and 0's
  compressed <- round(compressed)
  # Converting the matrix to a co-ordinates matrix
  coords <- matrix(0,nrow=length(compressed),ncol=3)
  l <- 1
  for(j in 1:nrow(compressed)) {
    for(k in 1:ncol(compressed)) {
      coords[l,] <- c(j,maxdim-k,compressed[j,k])
      l <- l+1
    }
  }
  gauss.mix <- mvnormalmixEM(coords[coords[,3]==1,1:2],k=ks[i])
  out[[i]] <- c(letter = all.letters[i],
                lambda = gauss.mix$lambda,
                mu = gauss.mix$mu,
                sigma = gauss.mix$sigma)
}
plot(coords[,1],coords[,2],col=coords[,3],asp=1)
