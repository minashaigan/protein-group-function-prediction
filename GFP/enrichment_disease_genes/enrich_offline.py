import sys, re
import numpy as np
import urllib2

# -------
# Qing Wei

gocountHash = {}

def parseArgs(args):
  """Parses arguments vector, looking for switches of the form -key {optional value}.
  For example:
    parseArgs([ 'main.py', '-p', 5 ]) = {'-p':5 }"""
  args_map = {}
  curkey = None
  for i in xrange(1, len(args)):
    if args[i][0] == '-':
      args_map[args[i]] = True
      curkey = args[i]
    else:
      assert curkey
      args_map[curkey] = args[i]
      curkey = None
  return args_map

def validateInput(args):
    args_map = parseArgs(args)
    #human organism id
    organism_id = 9606
    input_file = None
    
    if '-o' in args_map:
      organism_id = int(args_map['-o'])
    if '-f' in args_map:
      input_file = args_map['-f']
    return [organism_id, input_file]
    
def nCr(n, r):
    if r > n / 2.0:
        r = n - r
    ans = 1
    i = 1
    while i <= r:
        tmp = n - r + i
        ans = ans * tmp
        tmp = i
        ans = ans / tmp
        i += 1
    
    return ans
    
def f(k, N, m, n):
    a = nCr(m, k)
    b = nCr(N-m, n-k)
    c = nCr(N, n)
    
    res = a * b
    res = res / float(c)
    return res
    pass

def calculate_enrich(org, go, k, n):
    total = int(18687)       #number of annotated proteins in org = human
    m = 0               #number of proteins in org annotated with go
    if gocountHash.has_key(go):
        m = int(gocountHash[go])
    
#    print total, m
#    print go, m
    _sum = 0
    i = k
    
    while i <= n:
        _sum += f(i, total, m, n)
        i += 1
    
    return _sum

def main():
    usage = "python enrich_offline.py [-o organism_id] -f input"
    arguments = validateInput(sys.argv)
    organism_id, input_file = arguments
    gohash = {}
    proteinhash = {}
    pro_num = 0

    #load gocount

    f = open("data/Human_uniprot_annot_full+CustomSlimTermsIc0_online_GOcount.txt", "r")
    for line in f:
        go, count = line.split("\t")
        if not gocountHash.has_key(go):
            gocountHash[go] = count
    f.close();
    
    fout = open("GO-above-pval0.01.txt", "w");
    #print organism_id, input_file
    if input_file == None:
        print usage
    else:
        f = open(input_file)
        for line in f:
            tmp = line[:-1]
            parts = tmp.split(' ')
            if parts != None and len(parts) == 2:
                goterms = parts[1][:-1].split(',')
                #print goterms
                for go in goterms:
                    if len(go) == 10:
                        if go in gohash:
                            gohash[go] += 1
                        else:
                            gohash[go] = 1
                        if go in proteinhash:
                            proteinhash[go].append(parts[0])
                        else:
                            proteinhash[go] = []
                            proteinhash[go].append(parts[0])
            else:
                continue;
#                print "Error parsing the file!"
#                exit(0)
            pro_num += 1
            
    #print gohash
    #print proteinhash
    
    for go in gohash:
	pval = calculate_enrich(organism_id, go, gohash[go], pro_num);
        
    	if pval <= 0.01:
		print go + " " + str(pval)
		fout.write(go + ",")
       
    fout.close();
    f.close();
if __name__ == '__main__':
    main()



