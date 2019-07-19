library(apcluster);

args<-commandArgs(TRUE);
iter = as.numeric(args[1]);
networkName = args[2];


# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
networkName <- trim(networkName)

snf = as.matrix(read.table(networkName, header = TRUE));
funsim = as.matrix(read.table(paste("funsim_matrix", iter, sep = "_")), header = TRUE);

mean = (snf + funsim) / 2;
model = apcluster(mean);
write.table(mean, paste("apcluster_matrix", iter, sep = "_"), append = FALSE);
#attributes(model);

ncluster = length(attributes(model)$clusters);
write(ncluster, file = paste("apcluster_result", iter, sep = "_"), append = FALSE);

for(i in 1:ncluster)
{
	cur_cluster = attributes(model)$clusters[i];
	table = as.table(cur_cluster[[1]]);
	write.table(table, paste("apcluster_result", iter, sep = "_"), append = TRUE);
	
}

