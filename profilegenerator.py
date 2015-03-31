import collections
import itertools
import sys

filename=sys.argv[1]
MOTIF=sys.argv[2]
MOTIFSIZE=len(MOTIF)
MaxDEPTH=int(sys.argv[3])
MINIMUMPROB=float(sys.argv[4])##1.0/(10**4)
MININUMCOUNT=1
fd=open(filename)
lines=fd.readlines()
countbymajorallele=collections.defaultdict(list)
for line in lines:
    temp=line.strip().split('\t')
    t_major=int(temp[0])
    t_count=int(temp[2])
    countbymajorallele[t_major].append(t_count)
fd.close()
sumbymajorallele=collections.defaultdict(int)
for t_majorallele in countbymajorallele.keys():
    sumbymajorallele[t_majorallele]=sum(countbymajorallele[t_majorallele])

fd=open(filename)
##fd=open('PCRinclude.mono.A.bymajorallele')
lines=fd.readlines()
allmajor=collections.defaultdict(list)
for line in lines:
    temp=line.strip().split()
    if int(temp[0])%MOTIFSIZE==0:
        if (int(temp[2])/(sumbymajorallele[int(temp[0])]*1.0))>=MINIMUMPROB:
            if int(temp[2])>=MININUMCOUNT:
                allmajor[int(temp[0])].append(int(temp[1]))
##print allmajor
allkey=allmajor.keys()
allkey.sort()
#print allkey
keycount=0
combinelist_collection=[]
for dummycount in range(len(allkey)-1):
    pair1,pair2=allkey[keycount],allkey[keycount+1]
    pair1list=allmajor[pair1]
    pair2list=allmajor[pair2]
    #print pair1list,pair2list
    pair1list.extend(pair2list)
    combinelist=list(set(pair1list))
    combinelist.sort()
    ##print combinelist
    combinelist_collection.append(tuple(combinelist))
    keycount+=1
combinelist_collection=list(set(combinelist_collection))
newcombinelist_collection=combinelist_collection[:]
#combinelist_collection=set(combinelist_collection)
for smallset1 in combinelist_collection:
    for smallset2 in combinelist_collection:
        if set(smallset1).issubset(set(smallset2)) and smallset1 != smallset2:
            newcombinelist_collection.remove(smallset1)
            break
##print combinelist_collection
    
for depth in range(2,MaxDEPTH+1):
    for member_list in newcombinelist_collection:
        for member in itertools.combinations_with_replacement(member_list,depth):
            print 'chr'+'\t'+','.join(map(str,member))+'\t'+MOTIF
                
    
