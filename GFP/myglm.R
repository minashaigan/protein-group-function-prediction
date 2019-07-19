library(MASS)
library(arm);

read_data = read.table("features.txt",header=TRUE)
attach(read_data)
data = read_data;

#FAM_pos FAM_neg neiE_pos neiE_neg funsim_pos funsim_neg label

glmfit = bayesglm(label~ neiE_pos+neiE_neg+funsim_pos+funsim_neg+FAM_pos+FAM_neg,family=binomial(link=logit), data=data)


coeff = glmfit$coefficients;
x = as.matrix(coeff);
paramC = x[1:dim(x)[1],1];
paramC[is.na(paramC)] <- 0;

params = paramC;

write.table(params,"params.txt")
