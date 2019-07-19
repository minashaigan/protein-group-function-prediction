library(MASS);
#library(glm2);

read_data = read.table("features.txt",header=TRUE)
attach(read_data)
#prot = read_data$prot
#label = read_data$label
#prot = as.data.frame(prot)
#label = as.data.frame(label)
#data1 = as.data.frame(scale(read_data[,3:dim(read_data)[2]-1]))	
#data = cbind(prot,data1,label)
data = read_data;

#model = dget("model.rdata")

feature <- data[,3:dim(data)[2]-1];
feature <- as.matrix(feature);
params = as.matrix(read.table("params.txt"));
param = params[2:dim(params)[1], 1];
intercept = params[1, 1];

#l = as.matrix(predict.glm(model,data));

l = (feature %*% param)	+ intercept;
log_odd_ratio = exp(l);
prob1 = log_odd_ratio/(1+log_odd_ratio);
prob1 = as.matrix(pmax(0,prob1));

#if(is.na(prob1))
#{
#	write.table(params, "params_na.txt");
#	write.table(data, "data_na.txt");
#}
# write prob1 of each iteration of MH
write(as.vector(prob1), file = "prob1.txt", append = TRUE, ncolumns = dim(prob1)[1]);

threshold = as.matrix(runif(dim(prob1)[1], min=0, max=1));
threshold = prob1 - threshold;
pass <- as.matrix(which(threshold>=0));

label <- matrix(0,dim(prob1)[1],1)
label[pass, 1] <- 1;

final <- cbind(as.matrix(data[,1]), label);
write.table(final, file = "upd_label.txt", row.names = FALSE, col.names = FALSE, quote = FALSE);

