import sys
import collections
import math
SAMPLINGCOL=11
ALLELE1COL=7
ALLELE2COL=8
SIGNCOL=4
readprofileCOL=2
motifCOL=3
filaname=sys.argv[1]
fd=open(filaname)
lines=fd.readlines()
binomialcombine=collections.defaultdict(list)
for line in lines:
    temp=line.strip().split('\t')
    allelelist=[]
    allelelist.append(int(temp[ALLELE1COL-1]))
    allelelist.append(int(temp[ALLELE2COL-1]))
    allelelist.sort()
    #allelelist=map(str,allelelist)
    alleleave=str(allelelist[0])+'_'+str(allelelist[1])
    #alleleave=str(sum(allelelist)/2.0)
    ##alleleave=str(allelelist[0])+'_'+str(allelelist[1])
    totalcov=len(temp[readprofileCOL-1].split(','))
    motif=temp[motifCOL-1]
    samplingvalue=float(temp[SAMPLINGCOL-1])
    SIGN=1 
    binomialcombine[(totalcov,alleleave,motif)].append(SIGN*samplingvalue)
allkeys= binomialcombine.keys()
allkeys.sort()
##print allkeys
print 'read_depth'+'\t'+'allele'+'\t'+'heterozygous_prob'+'\t'+'motif'
for key in allkeys:
    ##templist=[str(key[0]),key[1],str(sum(binomialcombine[key])),key[2],str(map(str,(binomialcombine[key])))]
    templist=[str(key[0]),key[1],str(sum(binomialcombine[key])),key[2]]

    print '\t'.join(templist)
#print allkeys#,binomialcombine
    
    
        
