true = as.matrix(read.table(file = "true"));
pred = as.matrix(read.table(file = "pred"));

truex = true[7,];
predx = pred[1,]

tp = 0;
t = 0;
p = 0;
for(i in 1:length(truex))
{
	if(truex[i] == 1)
	{
		t = t + 1;
		if(predx[i] == 1)
		{
			tp = tp + 1;
		}
	}

	if(predx[i] == 1)
	{
		p = p + 1;
	}
}

prec = tp / p;
recall = tp / t;
