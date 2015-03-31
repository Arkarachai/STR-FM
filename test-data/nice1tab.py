import sys
fd=open(sys.argv[1])
lines=fd.readlines()
for line in lines:
    temp=line.strip().split()
    print '\t'.join(temp)