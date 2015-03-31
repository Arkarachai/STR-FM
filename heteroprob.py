### import libraries ###
import sys
import collections, math
import heapq
import itertools



### basic function ###
def permuterepeat(n,rlist):
    f = math.factorial
    nfac=f(n)
    rfaclist=[f(i) for i in rlist]
    for rfac in rfaclist:
        nfac=nfac/rfac
    return nfac

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
    
def averagelist(a,b,expectedlevelofminor):
    product=[]
    for i in range(len(a)):
        product.append((1-expectedlevelofminor)*a[i]+expectedlevelofminor*b[i])
  
    return product
        
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

def motifsimplify(base):
    '''str--> str
    '''
    motiflength=len(base)
    temp=list(set(ALLMOTIF[motiflength]).intersection(set(makeallpossible(base))))
    
    return temp[0]

def majorallele(seq):
    binseq=list(set(seq))  
    binseq.sort(reverse=True)   # highly mutate mode
    #binseq.sort()              # majority mode
    storeform=''
    storevalue=0
    for i in binseq:
        if seq.count(i)>storevalue:
            storeform=i
            storevalue=seq.count(i)
            
    return int(storeform)

### decide global parameter ###
COORDINATECOLUMN=1
ALLELECOLUMN=2
MOTIFCOLUMN=3
inputname=sys.argv[1]
errorprofile=sys.argv[2]
EXPECTEDLEVELOFMINOR=float(sys.argv[3])
if EXPECTEDLEVELOFMINOR >0.5:
	try:
		errorexpectcontribution=int('a')
	except Exception, eee:
		print eee
		stop_err("Expected contribution of minor allele must be at least 0 and not more than 0.5")
MINIMUMMUTABLE=0 ###1.2*(1.0/(10**8))  #http://www.ncbi.nlm.nih.gov/pubmed/22914163 Kong et al 2012


## Fixed global variable
ALLREPEATTYPE=[1,2,3,4]
ALLREPEATTYPENAME=['mono','di','tri','tetra']
monomotif=['A','C']
dimotif=['AC','AG','AT','CG']
trimotif=['AAC','AAG','AAT','ACC','ACG','ACT','AGC','AGG','ATC','CCG']
tetramotif=['AAAC','AAAG','AAAT','AACC','AACG','AACT','AAGC','AAGG','AAGT','AATC','AATG','AATT',\
'ACAG','ACAT','ACCC','ACCG','ACCT','ACGC','ACGG','ACGT','ACTC','ACTG','AGAT','AGCC','AGCG','AGCT',\
'AGGC','AGGG','ATCC','ATCG','ATGC','CCCG','CCGG','AGTC']
ALLMOTIF={1:monomotif,2:dimotif,3:trimotif,4:tetramotif}
monorange=range(5,60)
dirange=range(6,60)
trirange=range(9,60)
tetrarange=range(12,80)
ALLRANGE={1:monorange,2:dirange,3:trirange,4:tetrarange}

#########################################
######## Prob calculation sector ########
#########################################
def multinomial_prob(majorallele,STRlength,motif,probdatabase):
    '''int,int,str,dict-->int
    ### get prob for each STRlength to be generated from major allele
    '''
    #print (majorallele,STRlength,motif)
    prob=probdatabase[len(motif)][motif][majorallele][STRlength]
    return prob

################################################
######## error model database sector ###########
################################################

## structure generator
errormodeldatabase={1:{},2:{},3:{},4:{}}
sumbymajoralleledatabase={1:{},2:{},3:{},4:{}}
for repeattype in ALLREPEATTYPE:
    for motif in ALLMOTIF[repeattype]:
        errormodeldatabase[repeattype][motif]={}
        sumbymajoralleledatabase[repeattype][motif]={}
        for motifsize1 in ALLRANGE[repeattype]:
            errormodeldatabase[repeattype][motif][motifsize1]={}
            sumbymajoralleledatabase[repeattype][motif][motifsize1]=0
            for motifsize2 in ALLRANGE[repeattype]:
                errormodeldatabase[repeattype][motif][motifsize1][motifsize2]=MINIMUMMUTABLE
#print errormodeldatabase
## read database

## get read count for each major allele
fd=open(errorprofile)
lines=fd.readlines()
for line in lines:
    temp=line.strip().split('\t')
    t_major=int(temp[0])
    t_count=int(temp[2])
    motif=temp[3]
    sumbymajoralleledatabase[len(motif)][motif][t_major]+=t_count
fd.close()
##print sumbymajoralleledatabase

## get probability
fd=open(errorprofile)
lines=fd.readlines()
for line in lines:
    temp=line.strip().split('\t')
    t_major=int(temp[0])
    t_read=int(temp[1])
    t_count=int(temp[2])
    motif=temp[3]
    if sumbymajoralleledatabase[len(motif)][motif][t_major]>0:
        errormodeldatabase[len(motif)][motif][t_major][t_read]=t_count/(sumbymajoralleledatabase[len(motif)][motif][t_major]*1.0)
        #errormodeldatabase[repeattype][motif][t_major][t_read]=math.log(t_count/(sumbymajorallele[t_major]*1.0))
        
    #else:
    #    errormodeldatabase[repeattype][motif][t_major][t_read]=0
fd.close()
#print errormodeldatabase    
#print math.log(100,10)
#########################################
######## input reading sector ###########
#########################################



fd = open(inputname)
##fd=open('sampleinput_C.txt')
lines=fd.xreadlines()
for line in lines:
    i_read=[]
    i2_read=[]
    temp=line.strip().split('\t')
    i_coordinate=temp[COORDINATECOLUMN-1]
    i_motif=motifsimplify(temp[MOTIFCOLUMN-1])
    i_read=temp[ALLELECOLUMN-1].split(',')
    i_read=map(int,i_read)
    depth=len(i_read)
    heteromajor1=int(temp[6])
    heteromajor2=int(temp[7])

### calculate the change to detect combination (using error profile)
    heterozygous_collector=0  
    alist=[multinomial_prob(heteromajor1,x,i_motif,errormodeldatabase)for x in i_read]
    blist=[multinomial_prob(heteromajor2,x,i_motif,errormodeldatabase)for x in i_read]
      
    ablist=averagelist(alist,blist,EXPECTEDLEVELOFMINOR)
       
    if 0 in ablist:
        continue
    heterozygous_collector=reduce(lambda y, z: y*z,ablist )

### prob of combination (using multinomial distribution)
    frequency_distribution=[len(list(group)) for key, group in itertools.groupby(i_read)]
    ## print frequency_distribution
    expandbypermutation=permuterepeat(depth,frequency_distribution)

    print line.strip()+'\t'+str(heterozygous_collector)+'\t'+str(expandbypermutation)+'\t'+str(expandbypermutation*heterozygous_collector)+'\t'+str(depth)
