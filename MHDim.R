MetropolisHastingC = function(SamplesN, densityInt, sigma, initial) {
  OldPos=initial;
  OutputSamples= matrix(0,SamplesN,length(initial));
  for (i in 1:SamplesN)
{
    Proposal= rmvnorm(1,OldPos,as.matrix(sigma));
    ProposedValue=densityInt(as.vector(Proposal));
    OldValue=densityInt(OldPos);
    if (runif(1)<= (ProposedValue/OldValue))
{
      OldPos=Proposal;
    }
OutputSamples[i,]=OldPos
  }
return(OutputSamples)
}
MetropolisHastingC(1000,function(j){dnorm(j,0,1)},1,0)