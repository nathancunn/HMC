Samples=matrix(0,100,2);
for (i in 1:100) {
S=runif(1);
if (S<=0.29){
  Samples[i,]=rmvnorm(1,c(524,535.4),matrix(c(5.87,-9,-9,30.3),2,2))
}
if (S>0.29 & S<=0.5){
  Samples[i,]=rmvnorm(1,c(531.5,534.9),matrix(c(3.25,5.03,5.03,27),2,2))
}

if (S>0.5 & S<=0.7){
  Samples[i,]=rmvnorm(1,c(536.7,534.5),matrix(c(2.66,-4.45,-4.45,26.8),2,2))
}

if (S>0.7 & S<=1){
  Samples[i,]=rmvnorm(1,c(543.2,535.26),matrix(c(5.574,8.46,8.46,30.6),2,2))
}
}