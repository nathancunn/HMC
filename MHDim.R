MetropolisHastingC = function(TamplesN, densityInt, sigma, initial,thin) {
 # now=Sys.time()
  #Multivariate Normal with Samples, integration, (sigma in this case for 2-dim),initial conditions, thinning (1 if no thinning)
  SamplesN=thin*TamplesN;
  #Works out how much thinning
  OldPos=initial;
  OutputSamples= matrix(0,SamplesN,length(initial));
 # pb <- txtProgressBar(min = 0, max = SamplesN, style = 3);
  #Allocation
  Dimensions=length(initial);
  for (i in 1:SamplesN)
{
    Proposal<- as.vector(OldPos)+as.vector(runif(Dimensions,-sigma,sigma))
    #Proposal Value using uniform proposal (2-dimesional))
      #rmvnorm(1,as.vector(OldPos),sigma);
    ProposedValue=densityInt(as.vector(Proposal));
    #Proposed Value Density
    OldValue=densityInt(as.vector(OldPos));
    #Old Value Density
    if (runif(1)<= (ProposedValue/OldValue))
      #MH rejection
{
      OldPos=Proposal;
    }
OutputSamples[i,]=OldPos
#setTxtProgressBar(pb, i)
  }
if (TamplesN==1)
{
  FinalOutputSamples=OutputSamples[thin,]
}
else {
#Size=dim(unique(OutputSamples))
FinalOutputSamples=OutputSamples[seq(1,dim(OutputSamples)[1],thin),] }
#Thinning
#plot(FinalOutputSamples[,1],FinalOutputSamples[,2],type='l',col=2,xlab="X",ylab="Y",xlim=c(-2,2),ylim=c(-2,2))
#points(FinalOutputSamples[,1],FinalOutputSamples[,2],pch=4)
#cat("rejectionRate=",1-(Size[1]/SamplesN))
#after=Sys.time()-now;
#print(after)
return(FinalOutputSamples)
#close(pb)
}
MetropolisHastingC(1000,function(j){dnorm(j,0,1)},1,0,10)

MetropolisHastingC(25,function(j){dmvnorm(j,c(0,0),matrix(c(1,0.95,0.95,1),2,2))},2.1,c(0,0),20)
MetropolisHastingC(25,function(j){dmvnorm(j,c(0,0),matrix(c(1,0.95,0.95,1),2,2))},0.4,c(0,0),20)

all.lets <- data.frame(key=1:54)
for (i in 1:54){
  all.lets$letter[i] <- A[[1]][[i]]$letter
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
  atts <- (length(A[[1]][[1+2*i]])-1)/3
  for(j in 1:atts){
    A[[1+(i*2)]][1+atts+j][[1]] <- A[[1+(i*2)]][1+atts+j][[1]]+c(0,+2)
  }
}

WordprintMH <- function(word,lastValueH,thin) {
  #Set up the word 
  letters <- wordform(word)
  n <- length(letters)
  eq_probs <- seq(1/n,1,length=n)
  letter <<- findInterval(runif(1),eq_probs)+1
  j <- letters[letter]
  displacement <- (0:(n-1))*18-20
  #Extract the words and customise the spacing 
  atts <- (length(A[[1]][[j]])-1)/3
  lambdas <- cumsum(A[[1]][[j]][2:(2+atts-1)])
  group <<-  findInterval(runif(1),lambdas)+1
  mean.group <- A[[1]][[j]][1+atts+group][[1]]+c(displacement[letter],0)
  var.group <- A[[1]][[j]][1+2*atts+group][[1]]
  #Shift Mean and extract variance
  if (lastValueH[letter,group,1]==0 & lastValueH[letter,group,2]==0)
  {lastValueH[letter,group,]=as.vector(mean.group)}
  #If no starting value, use the mean
  FinalOutput=MetropolisHastingC(1, function(l) {dmvnorm(l,as.vector(mean.group)-c(3,3),var.group)},sqrt(sqrt(max(as.vector(var.group)))), as.vector(lastValueH[letter,group,]),thin) 
  #Use Metropolis Hasting algorithm, use the last position.
  #points(rmvnorm(1,mu = mean.group, sigma = var.group),cex=0.75,col=cols[letter]) 
  (list(FinalOutput,letter,group)) 
}
WordOutputAllMH= function(Word,SampleSize,thin=1){
  #Takes in a string, sample size.
  now=Sys.time()
  plot.new()
  FinalOutput=matrix(0,SampleSize,2);
  ColorsMemory=matrix(0,SampleSize,2)
  length.of.string <-nchar(Word);
  LastValueG <- array(0,c(length.of.string,6,2));
  #Allocation
  for (i in 1:SampleSize) {
    OutputOnce<-WordprintMH(Word,LastValueG,thin)
    #Metropolis Hasting using last update
    FinalOutput[i,]=OutputOnce[[1]]
    #Final output
    ColorsMemory[i,]=c(OutputOnce[[2]],OutputOnce[[3]])
    #Allocation of Color
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
  #Computes Limit
  #quartz()
  plot(FinalOutput[1,],cex=0.75,xlim=c(FinalOutputMinX-2,FinalOutputMaxX+2),ylim=c(FinalOutputMinY-2,FinalOutputMaxY+2),col=rgb(ColorsMemory[1,1]/length.of.string,ColorsMemory[1,2]/6,0.2),xlab="",ylab="",main="MH")
  #Plot first point
  #Sys.sleep(1)
  for (i in 2:SampleSize){
  points(FinalOutput[i,1],FinalOutput[i,2],cex=0.75,col=rgb(ColorsMemory[i,1]/length.of.string,ColorsMemory[i,2]/6,0.2),xlim=c(FinalOutputMinX-2,FinalOutputMaxX+2),ylim=c(FinalOutputMinY-2,FinalOutputMaxY+2)) }
  #Plot mixture and letter different colours.   
  #Sys.sleep(0.01)
  Size=dim(unique(FinalOutput))
  after=Sys.time()-now
  print(after)
  cat("Total Samples after Burn in:",SampleSize,"\nMean Of Samples=",apply(FinalOutput,2,mean),"\nCovariance=",cov(FinalOutput), "\nRejectionRate=",1-(Size[1]/SampleSize))
  return(FinalOutput)
}

JJ=MetropolisHastingC(1000,function(j){dmvnorm(j,rep(0,150),diag(seq(from=0.02,to=1,length=150)^2))},0.025, rep(0,150),60)
plot(MatrixMake[1],var(JJ[,1]),xlim=c(0,1),ylim=c(0,1),xlab="Real Variance of coordinate",ylab="Sample Variance",main="Metropolis Hasting")
for (t in 2:150) {
  points(MatrixMake[t],var(JJ[,t]))
}
lines(MatrixMake,MatrixMake)
