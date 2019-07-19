x = as.vector(scan(file = "tmp/group.txt", sep = " "));
y = as.vector(scan(file = "tmp/node.txt", sep = " "));

z = x | y;
final = as.numeric(z);

write(final, file = "tmp/or_result.txt", ncolumns = length(final));
