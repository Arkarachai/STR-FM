import sys
# remove all read that have impure microsat
# check only one line at a time


fd=open(sys.argv[1])
lines=fd.xreadlines()
##motifIx=int(sys.argv[2])
period=int(sys.argv[2])
tr_ref_seqIx=int(sys.argv[3])-1
##output=(sys.argv[4])
##fout=open(output,'w')
for line in lines:
    temp=line.strip().split('\t')
    temp=filter(None,temp)
    #motif=temp[motifIx]
    tr_ref_seq=temp[tr_ref_seqIx]
    ##period=len(motif)
    cand_motif=tr_ref_seq[:period]
    len_microsat=len(tr_ref_seq)
    expand_microsat_cand=cand_motif*(len_microsat/period) + cand_motif[:(len_microsat%period)]
    if tr_ref_seq == expand_microsat_cand:
    	print line.strip()
        ##print line.strip() >> fout