library(MASS);
library(arm);
#library(glm2);

read_data = read.table("features.txt",header=TRUE)
attach(read_data)
#prot = read_data$prot
#label = read_data$label
#data = scale(read_data[,3:dim(read_data)[2]-1])     
#prot = as.data.frame(prot)
#label = as.data.frame(label)
#data1 = as.data.frame(scale(read_data[,3:dim(read_data)[2]-1]))	
#data = cbind(prot,data1,label)
data = read_data;


feature <- data[,3:dim(data)[2]-1];
feature <- as.matrix(feature);

params = as.matrix(read.table("params.txt"));
param_current = params[2:dim(params)[1], 1];
intercept_current = params[1, 1];

#current_model = dget("model.rdata");

glmfit = bayesglm(label~ neiE_pos+neiE_neg+funsim_pos+funsim_neg+FAM_pos+FAM_neg,family=binomial(link=logit), data=data);
#glmfit = bayesglm(label~ neiE_pos+neiE_neg+funsim_pos+funsim_neg,family=binomial(link=logit), data=data)
#glmfit = bayesglm(label~ neiE_pos+neiE_neg,family=binomial(link=logit), data=data)
#cov = vcov(glmfit)

coeff = glmfit$coefficients;
params_new = as.matrix(coeff);
#params_new = as.matrix(mvrnorm(1, params_new,cov))
param_new = params_new[2:dim(params_new)[1],1];
param_new[is.na(param_new)] <- 0;
intercept_new = params_new[1, 1];


#param_new = mvrnorm(1, param_new,diag(length(param_new)))

#l_cur = as.matrix(predict.glm(current_model,data));
l_cur = (feature %*% param_current);	# + intercept_current;
lor_cur = exp(l_cur);
lor_cur = lor_cur/(1+lor_cur);
lor_cur[lor_cur == 1] <- 0.99;

#l_new = as.matrix(predict.glm(new_model,data));
l_new = (feature %*% param_new);	# + intercept_new;
lor_new = exp(l_new);
lor_new = lor_new/(1+lor_new);
lor_new[lor_new == 1] <- 0.99;
 
ones <- as.matrix(data$label);
ll_cur <- sum(log(lor_cur[which(ones == 1),1]));	#log likelihood
ll_cur <- ll_cur + sum(log(1-lor_cur[which(ones == 0),1]));
#ll_cur <- ll_cur + sum(log(lor_cur[which(ones == 0),1]));	#my change
ll_new <- sum(log(lor_new[which(ones == 1),1]));	#log likelihood
ll_new <- ll_new + sum(log(1-lor_new[which(ones == 0),1]));
#ll_new <- ll_new + sum(log(lor_new[which(ones == 0),1]));	#my change


acceptance_prob = min(0, ll_new-ll_cur);		#bcrf
#acceptance_prob = min(1, exp(ll_new)/exp(ll_cur));	#our
uniform_number = log(runif(1, min=0, max=1));

#if(is.na(acceptance_prob))
#{
#	acceptance_prob <- 0;
#}

if(acceptance_prob >= uniform_number)
{
	write(paste(acceptance_prob, ll_new-ll_cur, uniform_number, "accept",sep="_"), file = "metropolis.txt", append = TRUE);
	write.table(coeff,"params.txt");

#	param = as.matrix(param_new);
#	write(as.vector(param), file = "all_params.txt", append = TRUE, ncolumns = dim(param)[1]);
#	dput(new_model,file="model.rdata")
}else
{
	write(paste(acceptance_prob, ll_new-ll_cur, uniform_number, "reject",sep="_"), file = "metropolis.txt", append = TRUE);
#	param = as.matrix(param_current);
#	write(as.vector(param), file = "all_params.txt", append = TRUE, ncolumns = dim(param)[1]);
}
