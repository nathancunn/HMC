# Code for using letter plots
setwd("HMC/")

# Load the Gaussians specifying the letters
letter.models <- source("Output/output")[[1]]

# Creating a key matching the letters to their location in letter.models
all.letters <- data.frame(key=1:54)
for (i in 1:54){
all.letters$letter[i] <- letter.models[[i]]$letter
}

# A function to map letters from a word according to the key
wordform <- function(word) {
  # Split the word into individual characters
  chars <- strsplit(word, "")[[1]]
  out <- matrix(0, nrow=length(chars)) 
  for(i in 1:length(chars)) {
    if(any(all.letters$key[all.letters$letter==chars[i]])) {
    out[i,1] <- all.letters$key[all.letters$letter==chars[i]]
    } else out[i,1] <- 55
  }
  out
}



# Below this is probably not necessary

for(i in 1:26) {
  atts <- (length(letter.models[[1]][[1+2*i]])-1)/3
  for(j in 1:atts){
    letter.models[[1+(i*2)]][1+atts+j][[1]] <- letter.models[[1+(i*2)]][1+atts+j][[1]]+c(0,+2)
  }
}
  
Wordprint <- function(word,lastValueH) {
  #TRY44<<-lastValueH
  letters <- wordform(word)
  n <- length(letters)
  eq_probs <- seq(1/n,1,length=n)
  letter <<- findInterval(runif(1),eq_probs)+1
  j <- letters[letter]
  displacement <- (0:(n-1))*18-20
  #for (l in 1:SampleSize) {
  #for (k in 1:n) {
  #j=letters[k]
  atts <- (length(A[[1]][[j]])-1)/3
  lambdas <- cumsum(A[[1]][[j]][2:(2+atts-1)])
  group <<-  findInterval(runif(1),lambdas)+1
  mean.group <- A[[1]][[j]][1+atts+group][[1]]+c(displacement[letter],0)
  var.group <- A[[1]][[j]][1+2*atts+group][[1]]
  #TRY<<-lastValueH
  #TRY1<<-letter
  #TRY2<<-group
  if (lastValueH[letter,group,1]==0 & lastValueH[letter,group,2]==0)
    {lastValueH[letter,group,]=as.vector(mean.group)}
  FinalOutput=HMCVer1(1,densityA = function(l) dmvnorm(l,as.vector(mean.group),var.group),q = as.vector(lastValueH[letter,group,]),M=diag(2),epsilon = 0.05,L = 60,densitydiff = function(x) {solve(var.group)%*%(as.matrix(x)-as.matrix(mean.group))}, burnin = 0)
  #points(rmvnorm(1,mu = mean.group, sigma = var.group),cex=0.75,col=cols[letter])
(list(FinalOutput,letter,group))
}
WordOutputAll= function(Word,SampleSize,cols=NULL){
  plot.new()
  FinalOutput=matrix(0,SampleSize,2);
  #plot.window(xlim=c(-10,50), ylim=c(0,20))
  if(is.null(cols)) {cols=rep(1,n)}
  length.of.string <-nchar(Word);
  LastValueG <<- array(0,c(length.of.string,6,2));
  for (i in 1:SampleSize) {
  OutputOnce<-Wordprint(Word,LastValueG)
  FinalOutput[i,]=OutputOnce[[1]]
  LastValueG[OutputOnce[[2]],OutputOnce[[3]],]=FinalOutput[i,]
  #if (i==1){
  #plot(FinalOutput[i,],cex=0.75,col="black")}
  #else {
  #points(FinalOutput[i,],cex=0.75,col="black")
  #}
  }
  FinalOutputMinX=min(FinalOutput[,1])
  FinalOutputMinY=min(FinalOutput[,2])
  FinalOutputMaxX=max(FinalOutput[,1])
  FinalOutputMaxY=max(FinalOutput[,2])
  #quartz()
  #plot(FinalOutput[1,],cex=0.75,xlim=c(FinalOutputMinX-2,FinalOutputMaxX+2),ylim=c(FinalOutputMinY-2,FinalOutputMaxY+2))
  #Sys.sleep(1)
  #for (i in 2:SampleSize){
  #points(FinalOutput[i,1],FinalOutput[i,2],cex=0.75,col="black",xlim=c(FinalOutputMinX-2,FinalOutputMaxX+2),ylim=c(FinalOutputMinY-2,FinalOutputMaxY+2)) }
  #Sys.sleep(0.01)
  FinalOutput
}
plot(0:50,0:50,type="n",asp=1,bty="n",xaxt="n",yaxt="n",ann=F)
for(i in 1:10000) {
  WordPrint("Hello",cols=c(rep("#014FDC",4),rep("#01A401",2)))
  #WordPrint("mMC",cols=c(rep("#014FDC",4),rep("#01A401",2)))
}

