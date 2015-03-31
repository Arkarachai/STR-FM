import sys
# remove all read that have unmatch microsat
# check only one line at a time
def complement_base(read):
    collect=''
    for i in read:
        if i.upper()=='A':
            collect+='T'
        elif i.upper()=='T':
            collect+='A'
        elif i.upper()=='C':
            collect+='G'
        elif i.upper()=='G':
            collect+='C'
    return collect
   
def makeallpossible(read):
    collect=[]
    for i in range(len(read)):
        tmp= read[i:]+read[:i]
        collect.append(tmp)
        collect.append(complement_base(tmp))
    return collect


fd=open(sys.argv[1])
lines=fd.xreadlines()
firstcolumn=int(sys.argv[2])-1 #4
secondcolumn=int(sys.argv[3])-1 # 10
for line in lines:
    temp=line.strip().split('\t')
    temp=filter(None,temp)
    micro1=temp[firstcolumn]
    micro2=temp[secondcolumn]
    if micro1 in makeallpossible(micro2):
        print line.strip()