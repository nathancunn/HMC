#Sort out all the letters 
all.lets <- data.frame(key=1:54)
for (i in 1:54){
all.lets$letter[i] <- A[[1]][[i]]$letter
}
#Split up the rest
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
#work out mixtures
for (i in 1:26) {
  atts <- (length(A[[1]][[1+2*i]])-1)/3
  for(j in 1:atts){
    A[[1+(i*2)]][1+atts+j][[1]] <- A[[1+(i*2)]][1+atts+j][[1]]+c(0,+2)
  }
}
  
Wordprint <- function(word,lastValueH) {
  letters <- wordform(word)
  n <- length(letters)
  eq_probs <- seq(1/n,1,length=n)
  #Weights for letters
  letter <<- findInterval(runif(1),eq_probs)+1
  j <- letters[letter]
  displacement <- (0:(n-1))*18-20
  #Extract the words and customise the spacing 
  atts <- (length(A[[1]][[j]])-1)/3
  lambdas <- cumsum(A[[1]][[j]][2:(2+atts-1)])
  group <<-  findInterval(runif(1),lambdas)+1
  #Find the mixture in each letter
  mean.group <- A[[1]][[j]][1+atts+group][[1]]+c(displacement[letter],0)
  var.group <- A[[1]][[j]][1+2*atts+group][[1]]
  #Shift Mean and extract variance
  if (lastValueH[letter,group,1]==0 & lastValueH[letter,group,2]==0)
    {lastValueH[letter,group,]=as.vector(mean.group)}
  #If no starting value, use the mean
  FinalOutput=HMCVer1(1,densityA = function(l) dmvnorm(l,as.vector(mean.group),var.group),q = as.vector(lastValueH[letter,group,]),M=diag(2),epsilon = 0.05,L = 60,densitydiff = function(x) {solve(var.group)%*%(as.matrix(x)-as.matrix(mean.group))}, burnin = 0)
  #Use Metropolis Hasting algorithm, use the last position.
  #points(rmvnorm(1,mu = mean.group, sigma = var.group),cex=0.75,col=cols[letter]) 
(list(FinalOutput,letter,group))
}
WordOutputAll= function(Word,SampleSize){
  plot.new()
  FinalOutput=matrix(0,SampleSize,2);
  ColorsMemory=matrix(0,SampleSize,2)
  #Allocation
  #plot.window(xlim=c(-10,50), ylim=c(0,20))
  length.of.string <-nchar(Word);
  LastValueG <- array(0,c(length.of.string,6,2));
  #Allocation
  for (i in 1:SampleSize){
  OutputOnce<-Wordprint(Word,LastValueG)
  #Metropolis Hasting using last update
  FinalOutput[i,]=OutputOnce[[1]]
  #Final output
  ColorsMemory[i,]=c(OutputOnce[[2]],OutputOnce[[3]])
  LastValueG[OutputOnce[[2]],OutputOnce[[3]],]=FinalOutput[i,]
  #Remember last position for letter and mixture
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
  #Range Limit
  #quartz()
  plot(FinalOutput[1,],cex=0.75,xlim=c(FinalOutputMinX-2,FinalOutputMaxX+2),ylim=c(FinalOutputMinY-2,FinalOutputMaxY+2),col=rgb(ColorsMemory[1,1]/length.of.string,ColorsMemory[1,2]/6,0.2))
  #Plot first plot
  #Sys.sleep(1)
  for (i in 2:SampleSize){
  points(FinalOutput[i,1],FinalOutput[i,2],cex=0.75,col=rgb(ColorsMemory[i,1]/length.of.string,ColorsMemory[i,2]/6,0.2),xlim=c(FinalOutputMinX-2,FinalOutputMaxX+2),ylim=c(FinalOutputMinY-2,FinalOutputMaxY+2)) }
  #Plot mixture and letter in different colours
  #Sys.sleep(0.01)
  FinalOutput
  Size=dim(unique(FinalOutput))
  cat("Total Samples after Burn in:",SampleSize,"\nMean Of Samples=",apply(FinalOutput,2,mean),"\nCovariance=",cov(FinalOutput), "\nRejectionRate=",1-(Size[1]/SampleSize))
}

plot(0:50,0:50,type="n",asp=1,bty="n",xaxt="n",yaxt="n",ann=F)
for(i in 1:10000) {
  WordPrint("Hello",cols=c(rep("#014FDC",4),rep("#01A401",2)))
  #WordPrint("mMC",cols=c(rep("#014FDC",4),rep("#01A401",2)))
}

