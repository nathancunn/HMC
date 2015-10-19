MetropolisHastingC = function(TamplesN, densityInt, sigma, initial,thin) {
  SamplesN=thin*TamplesN;
  OldPos=initial;
  OutputSamples= matrix(0,SamplesN,length(initial));
  for (i in 1:SamplesN)
{
    Proposal<- as.vector(OldPos)+c(runif(1,-sigma,sigma),runif(1,-sigma,sigma))
      #rmvnorm(1,as.vector(OldPos),sigma);
    ProposedValue=densityInt(as.vector(Proposal));
    OldValue=densityInt(as.vector(OldPos));
    if (runif(1)<= (ProposedValue/OldValue))
{
      OldPos=Proposal;
    }
OutputSamples[i,]=OldPos
  }
if (TamplesN==1)
{
  FinalOutputSamples=OutputSamples[thin,]
}
else {
FinalOutputSamples=OutputSamples[seq(1,dim(OutputSamples)[1],thin),] }
return(FinalOutputSamples)
}
MetropolisHastingC(1000,function(j){dnorm(j,0,1)},1,0,10)

MetropolisHastingC(1000,function(j){dmvnorm(j,c(0,0),matrix(c(1,0.9,0.9,1),2,2))},diag(2),c(0,0),10)


all.lets <- data.frame(key=1:54)
for(i in 1:54){
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

WordprintMH <- function(word,lastValueH) {
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
  FinalOutput=MetropolisHastingC(1, function(l) {dmvnorm(l,as.vector(mean.group)-c(3,3),var.group)},sqrt(sqrt(max(as.vector(var.group)))), as.vector(lastValueH[letter,group,]),1) 
  #points(rmvnorm(1,mu = mean.group, sigma = var.group),cex=0.75,col=cols[letter])
  (list(FinalOutput,letter,group))
}
WordOutputAllMH= function(Word,SampleSize,cols=NULL){
  plot.new()
  FinalOutput=matrix(0,SampleSize,2);
  #plot.window(xlim=c(-10,50), ylim=c(0,20))
  if(is.null(cols)) {cols=rep(1,n)}
  length.of.string <-nchar(Word);
  LastValueG <<- array(0,c(length.of.string,6,2));
  for (i in 1:SampleSize) {
    OutputOnce<-WordprintMH(Word,LastValueG)
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
  for (i in 2:SampleSize){
  points(FinalOutput[i,1],FinalOutput[i,2],cex=0.75,col="black",xlim=c(FinalOutputMinX-2,FinalOutputMaxX+2),ylim=c(FinalOutputMinY-2,FinalOutputMaxY+2)) }
  #Sys.sleep(0.01)
  Size=dim(unique(FinalOutput))
  cat("Total Samples after Burn in:",SampleSize,"\nMean Of Samples=",apply(FinalOutput,2,mean),"\nCovariance=",cov(FinalOutput), "\nRejectionRate=",1-(Size[1]/SampleSize))
  return(FinalOutput)
}
