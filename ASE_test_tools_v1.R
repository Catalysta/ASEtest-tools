setwd("~/Jobs/STSI")

#p-val is NA for anything less than 8 total count
#output colu,m of pvals and adj pvals using Benjamini Hoschberg

#read in and plot data
dat<-read.table("Data/HG00096.1.M_111124_6_out.txt", header = TRUE)
plot(dat$refCount,dat$altCount, main = "Reference Count vs. Alternate Count",
     xlab = "Reference Count", ylab = "Alternative Count")
abline(a = 0, b = 1, col = "blue")

#read in and plot second set of data
dat2<-read.table("Data/HG00097.7.M_120219_2_out.txt", header = TRUE)
plot(dat$refCount,dat$altCount, main = "Reference Count vs. Alternate Count",
     xlab = "Reference Count", ylab = "Alternative Count")
abline(a = 0, b = 1, col = "blue")

#The following function performs a binomial test on two columns of count data (eh1 and eh2)
#taken from a raw .txt file of ASE data. Filepath for the raw file must be entered in
#quotation marks (""). The default binomial parameter to test is 0.5, while the default FDR is 0.05.
#The function will return 
#the entire binomial test output.
ASE.binom.test<-function(filepath, output_columns, eh1 = "refCount", eh2 = "altCount", prob = 0.5, FDR = 0.05){
  dat<-read.table(filepath, header = TRUE)
  output<-dat[,output_columns]
  eh1.zero.ind <- which(dat[eh1] == 0)
  eh2.zero.ind <- which(dat[eh2] == 0)
  #dat[eh1.zero.ind,eh1]<-.1 (possibly don't need these lines since no zeros in total counts greater than 8)
  #dat[eh2.zero.ind,eh2]<-.1
  for (i in 1:nrow(dat)){
    if (dat$totalCount[i]>8){
      test<-binom.test(c(dat[i,eh1],dat[i,eh2]), p = prob)
      output$p.val[i]<-test$p.value
    }
    else{output$p.val[i]<-NA}
    #if (output$p.val[i]<0.05){
    #  output$outlier[i]<-TRUE
    #}
    #else{output$outlier[i]<-FALSE}
  }
  #Carry out Benjamini-Hochberg procedure to get adjusted p-values
  ordered.pvals<-sort(output$p.val)
  m = length(ordered.pvals)
  for (i in 1:nrow(dat)){
    output$adj.pval[i]<-ordered.pvals[m]
    for (j in m:1){
      if (output$adj.pval[i]>=ordered.pvals[j]){
        output$adj.pval[i]<-ordered.pvals[j]*(m/j)
      }
      else{break}
    }
  }
  return(output)
}
i = 77435
eh1 = "refCount"
eh2 = "altCount"
prob = .5


#Look for outliers based on result of tests
for (i in 1:nrow(dat)){
  dat$prop[i]<-dat$altCount[i]/(dat$refCount[i]+dat$altCount[i])
  if ((dat$prop[i]==0) | (dat$prop[i]==1)){
    dat$outlier[i]<-TRUE
  }
  else{dat$outlier[i]<-FALSE}
}

#sub all zeros in count by .1

#Make new plot to show proportions outside of the confidence interval
plot(dat$refCount, dat$altCount, log="xy", main = "Reference Count vs. Alternate Count",
     xlab = "Reference Count", ylab = "Alternative Count",
     col = ifelse(dat$outlier==TRUE,'red','black'), pch = 19)
abline(a = 0, b = 1, col = "blue")


