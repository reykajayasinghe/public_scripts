###Reyka Jayasinghe
###November 24th, 2015
###rjayasin@genome.wustl.edu

#reference: http://davetang.org/muse/2013/01/30/sequence-logos-with-r/

###Script to make PWM for splice site consensus sequences and plot

###Install and Load appropriate packages:
source("http://bioconductor.org/biocLite.R")
library(grid)
library(seqLogo)
biocLite("seqLogo")
#Load data
Sequenceinformation<- read.table("~/Desktop/context/acceptor_same", quote="\"", comment.char="")
row.names(Sequenceinformation)=Sequenceinformation$V1
Sequenceinformation$V1 <- NULL

  
df <- as.data.frame(t(Sequenceinformation))

#define function that divides the frequency by the row sum i.e. proportions
proportion <- function(x){
  rs <- sum(x);
  return(x / rs);
}

#create position weight matrix
pwm <- apply(df, 1, proportion)
pwm <- makePWM(pwm)
seqLogo(pwm)
