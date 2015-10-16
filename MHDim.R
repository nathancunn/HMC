MetropolisHastingC = function(TamplesN, densityInt, sigma, initial,thin) {
  SamplesN=thin*TamplesN;
  OldPos=initial;
  OutputSamples= matrix(0,SamplesN,length(initial));
  for (i in 1:SamplesN)
{
    Proposal<- rmvnorm(1,as.vector(OldPos),sigma);
    ProposedValue=densityInt(as.vector(Proposal));
    OldValue=densityInt(as.vector(OldPos));
    if (runif(1)<= (ProposedValue/OldValue))
{
      OldPos=Proposal;
    }
OutputSamples[i,]=OldPos
  }
FinalOutputSamples=OutputSamples[seq(1,dim(OutputSamples)[1],thin),]
return(FinalOutputSamples)
}
MetropolisHastingC(1000,function(j){dnorm(j,0,1)},1,0,10)

MetropolisHastingC(1000,function(j){dmvnorm(j,c(0,0),matrix(c(1,0.9,0.9,1),2,2))},diag(2),c(0,0),10)