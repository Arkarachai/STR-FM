#!/usr/bin/env python
# pair_fetch_DNA_ff.py
# Function: filter microsat and flanking region by quality score;
# remove read with any base that has lower quality score than "quality_require" within "flanking_base" and convert from snoope to fastq
# Note that require flanking length need to be screen by Bob snoope script first

# Author: Arkarachai Fungtammasan
# Version 1.0.0 (15 July 2012)
# Input format: length_of_repeat[0] 	 left_flank_length[1]	right_flank_length[2]	repeat_motif[3]	hamming_distance[4]	read_name[5]	read_sequence[6]	read_quality[7]
# Output format: two fastq file. First file contain left flank. Second file contain right flank.
# Command: python pair_fetch_DNA_ff.py input.txt

import sys
from galaxy import eggs

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()
    
# read file name


	
filename=sys.argv[1]
L_filename=sys.argv[2]
R_filename=sys.argv[3]
quality_require=sys.argv[4]
flanking_base=sys.argv[5]
try:
	quality_require=int(quality_require)
	flanking_base=int(flanking_base)
except Exception, eee:
	print eee
	stop_err("Quality score cutoff and Length of flanking regions that require quality screening must be integer")
	
fd=open(filename)
fdd1=open(L_filename,'w')
fdd2=open(R_filename,'w')
lines=fd.xreadlines()
for line in lines:
    temp=line.strip().split('\t')
    temp=filter(None,temp)
    #get index
    left_flank=(0,int(temp[1]))
    microsat=(int(temp[1]),int(temp[1])+int(temp[0]))
    right_flank=(int(temp[1])+int(temp[0]),int(temp[1])+int(temp[0])+int(temp[2]))
    flag=0
    #filter length of left and right flank
    if (right_flank[1]-right_flank[0])<flanking_base:
    	continue
    if (left_flank[1]-left_flank[0])<flanking_base:
    	continue
    #filter quality score
    for i in temp[7][microsat[0]-flanking_base:microsat[1]+flanking_base]:
        if ord(i)<(quality_require+33):
            flag=1
        else:
            flag=flag
    #print out to seperated files
    if flag ==0:
        newname= temp[5]##+'_'+temp[3]+'_'+temp[0]
        fdd1.writelines('@'+newname+'\n')
        fdd2.writelines('@'+newname+'\n')
        fdd1.writelines(temp[6][left_flank[0]:left_flank[1]]+'\n')
        fdd2.writelines(temp[6][right_flank[0]:right_flank[1]]+'\n')
        fdd1.writelines('+'+newname+'\n')
        fdd2.writelines('+'+newname+'\n')
        fdd1.writelines(temp[7][left_flank[0]:left_flank[1]]+'\n')
        fdd2.writelines(temp[7][right_flank[0]:right_flank[1]]+'\n')

fd.close()
fdd1.close()
fdd2.close()


