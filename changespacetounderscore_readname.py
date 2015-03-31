import sys
fd=open(sys.argv[1])
output=open(sys.argv[2],'w')
columntochange=int(sys.argv[3])-1  # default is 6-1=5
lines=fd.xreadlines()
for line in lines:
	temp=line.strip().split('\t')
	temp=filter(None,temp)
	temp2=temp[columntochange].replace(' ','_')
	product=temp[:columntochange]
	product.append(temp2)
	product.extend(temp[columntochange+1:])
	output.writelines('\t'.join(product)+'\n')
fd.close()
output.close()