### import libraries ###
import sys
import collections, math
import heapq


    


### basic function ###
def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()
        
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
  ##(0.01-0.5)
MINIMUMMUTABLE=1.2*(1.0/(10**8))  #http://www.ncbi.nlm.nih.gov/pubmed/22914163 Kong et al 2012


## Fixed global variable
inputname=sys.argv[1]
errorprofile=sys.argv[2]
Genotypingcorrected=sys.argv[3]
EXPECTEDLEVELOFMINOR=float(sys.argv[4])
if EXPECTEDLEVELOFMINOR >0.5:
	try:
		expected_contribution_of_minor_allele=int('expected_contribution_of_minor_allele')
	except Exception, eee:
		print eee
		stop_err("Expected contribution of minor allele must be at least 0 and not more than 0.5")
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

#########################################
######## input reading sector ###########
#########################################
fdout=open(Genotypingcorrected,'w')

fd = open(inputname)

lines=fd.xreadlines()
for line in lines:
    i_read=[]
    i2_read=[]
    temp=line.strip().split('\t')
    i_coordinate=temp[COORDINATECOLUMN-1]
    i_motif=motifsimplify(temp[MOTIFCOLUMN-1])
    i_read=temp[ALLELECOLUMN-1].split(',')
    i_read=map(int,i_read)
    coverage=len(i_read)

### Evaluate 1 major allele ###    
    i_all_allele=list(set(i_read))
    i_major_allele=majorallele(i_read)
    f_majorallele=i_read.count(i_major_allele)
### Evaluate 2 major allele ### 
    if len(i_all_allele)>1:
        i2_read=filter(lambda a: a != i_major_allele, i_read)
        i_major2_allele=majorallele(i2_read)
        f_majorallele2=i_read.count(i_major2_allele)
        ### Evaluate 3 major allele ### 
        if len(i_all_allele)>2:
            i3_read=filter(lambda a: a != i_major2_allele, i2_read)
            i_major3_allele=majorallele(i3_read)
            f_majorallele3=i_read.count(i_major3_allele)
        ### No 3 major allele ### 
        elif len(i_all_allele)==2:
            i_major3_allele=i_major2_allele
    ### No 2 major allele ### 
    elif len(i_all_allele)==1:
        #i_major2_allele=majorallele(i_read)
        i_major2_allele=i_major_allele+len(i_motif)
        i_major3_allele=i_major2_allele
        #print line.strip()+'\t'+'\t'.join(['homo','only',str(i_major_allele),str(i_major_allele),'NA'])
        #continue
    else:
        print("no allele is reading")
        sys.exit()
    
## scope filter

#########################################
######## prob calculation sector ########
#########################################
    homozygous_collector=0
    heterozygous_collector=0

      
    alist=[multinomial_prob(i_major_allele,x,i_motif,errormodeldatabase)for x in i_read]
    blist=[multinomial_prob(i_major2_allele,x,i_motif,errormodeldatabase)for x in i_read]
    clist=[multinomial_prob(i_major3_allele,x,i_motif,errormodeldatabase)for x in i_read]
    
    ablist=averagelist(alist,blist,EXPECTEDLEVELOFMINOR)
    bclist=averagelist(blist,clist,EXPECTEDLEVELOFMINOR)
    aclist=averagelist(alist,clist,EXPECTEDLEVELOFMINOR)
    
    #print alist,blist,clist
    majora=sum([math.log(i,10) for i in alist])
    majorb=sum([math.log(i,10) for i in blist])    
    majorc=sum([math.log(i,10) for i in clist])
    homozygous_collector=max(majora,majorb,majorc)
    
    homomajor1=max([(majora,i_major_allele),(majorb,i_major2_allele),(majorc,i_major3_allele)])[1]
    homomajordict={i_major_allele:majora,i_major2_allele:majorb,i_major3_allele:majorc}
    
    majorab=sum([math.log(i,10) for i in ablist])
    majorbc=sum([math.log(i,10) for i in bclist])    
    majorac=sum([math.log(i,10) for i in aclist])
    heterozygous_collector=max(majorab,majorbc,majorac)
    bothheteromajor=max([(majorab,(i_major_allele,i_major2_allele)),(majorbc,(i_major2_allele,i_major3_allele)),(majorac,(i_major_allele,i_major3_allele))])[1]
    ##heteromajor1=max(bothheteromajor)
    ##heteromajor2=min(bothheteromajor)
    pre_heteromajor1=bothheteromajor[0]
    pre_heteromajor2=bothheteromajor[1]
    heteromajor1=max((homomajordict[pre_heteromajor1],pre_heteromajor1),(homomajordict[pre_heteromajor2],pre_heteromajor2))[1]
    heteromajor2=min((homomajordict[pre_heteromajor1],pre_heteromajor1),(homomajordict[pre_heteromajor2],pre_heteromajor2))[1]
    
    logratio_homo=homozygous_collector-heterozygous_collector
    
    if logratio_homo>0:
        fdout.writelines(line.strip()+'\t'+'\t'.join(['homo',str(logratio_homo),str(homomajor1),str(heteromajor1),str(heteromajor2)])+'\n')
    elif logratio_homo<0:
        fdout.writelines(line.strip()+'\t'+'\t'.join(['hetero',str(logratio_homo),str(homomajor1),str(heteromajor1),str(heteromajor2)])+'\n')
fd.close()
fdout.close()  
