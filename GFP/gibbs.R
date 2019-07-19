args<-commandArgs(TRUE);
iter = as.numeric(args[1]);
burn = ceiling(iter*0.2); 

prob = read.table("prob1.txt"); 
lag_prob = prob[burn:dim(prob)[1],]
lag_prob = as.matrix(lag_prob);		#for case where number of test data = 1
mean_prob = as.matrix(colMeans(lag_prob));

#make final inference based on mean_prob

threshold = as.matrix(runif(dim(mean_prob)[1], min=0, max=1));
threshold = mean_prob - threshold;
pass <- as.matrix(which(threshold>=0));

label <- matrix(0,dim(mean_prob)[1],1)
label[pass, 1] <- 1;

data = read.table("features.txt",header=TRUE);
attach(data);

final <- cbind(as.matrix(data[,1]), label);
final_prob <- cbind(as.matrix(data[,1]), as.matrix(mean_prob));
write.table(final, file = "upd_label.txt", row.names = FALSE, col.names = FALSE, quote = FALSE);
write.table(final_prob, file = "mean_prob", row.names = FALSE, col.names = FALSE, quote = FALSE);
