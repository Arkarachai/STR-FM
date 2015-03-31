def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()
        
def info2require(X,L,F,r):
    '''infodepth,readlength,flanksize,repeatlength
    '''
    return int(math.ceil((X*L*1.0)/(L-(1*((2*F)+r-1)))))
    
def poissondef(meancov,specificcov):
    nominator=1.0*(meancov**specificcov)*(math.e**(-1*meancov))
    denominator=math.factorial(specificcov)
    return nominator/denominator

def require2recommend(needprob,mindepth):
    i=mindepth
    reverseneedprob=1-needprob
    sumprob=1
    while sumprob>reverseneedprob: #mean cov
        sumprob=0
        for j in range(0,mindepth): #specific cov
            sumprob+=poissondef(i,j)
        i+=1
        
    return i-1

import sys,math

repeatlength=int(sys.argv[1])
flanksize=int(sys.argv[2])#20
readlength=int(sys.argv[3])#100
infodepth=int(sys.argv[4])#5
probdetection=float(sys.argv[5])#0.90

if probdetection >1:
    try:
        probvalue=int('probvalue')
    except Exception, eee:
        print eee
        stop_err("Proportion of genome to have certain locus specific must be between 0 and 1")

print 'repeat_length'+'\t'+'read_length'+'\t'+'informative_read_depth''\t'+'=locus_specific_sequencing_depth'+'\t'+'=genome_wide_sequencing_depth'
t_requiredepth=info2require(infodepth,readlength,flanksize,repeatlength)
t_recomendseq=require2recommend(probdetection,t_requiredepth)
preplotlist=[repeatlength,readlength,infodepth,t_requiredepth,t_recomendseq]
plotlist=map(str,preplotlist)
print '\t'.join(plotlist)

#print info2require(infodepth,readlength,flanksize,repeatlength)
#print poissondef(10,3)
#print require2recommend(0.90,80)
#informative_read_depth
#required_seq_depth
#recommend_seq_depth