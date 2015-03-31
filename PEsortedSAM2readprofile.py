#!/usr/bin/env python

import sys
from galaxy import eggs
import pkg_resources
pkg_resources.require( "bx-python" )
import bx.seq.twobit

##output columns: read_name chr prefix_start    prefix_end  TR_start    TR_end  suffix_start    suffix_end  TR_length   TR_sequence

samf = open(sys.argv[1],'r') #assumes sam file is sorted by readname
seq_path = sys.argv[2] #Path to the reference genome in 2bit format

##maxTRlength=int(sys.argv[4])
##maxoriginalreadlength=int(sys.argv[5])
maxTRlength=int(sys.argv[3])
maxoriginalreadlength=int(sys.argv[4])
outfile=sys.argv[5]
fout = open(outfile,'w')

twobitfile = bx.seq.twobit.TwoBitFile( file( seq_path ) )

skipped=0
while True:
    read = samf.readline().strip()
    if not(read): #EOF reached
        break
    if read[0] == "@":
        #print read
        continue
    mate = samf.readline().strip()
    if not(mate): #EOF reached
        break
    read_elems = read.split()
    mate_elems = mate.split()
    read_name = read_elems[0].strip()
    mate_name = mate_elems[0].strip()
    while True:
        if read_name == mate_name:
            break
        elif read_name != mate_name:
            #print >>sys.stderr, "Input SAM file doesn't seem to be sorted by readname. Please sort and retry."
            #break
            skipped += 1
            read = mate
            read_elems = mate_elems
            mate = samf.readline().strip()
            read_name = read_elems[0].strip()
            mate_name = mate_elems[0].strip()
            if not(mate): #EOF reached
                break
            mate_elems = mate.split()
    #extract XT:A tag
    #for e in  read_elems:
    #    if e.startswith('XT:A'):
    #        read_xt = e
    #for e in  mate_elems:
    #    if e.startswith('XT:A'):
    #        mate_xt = e
    #if 'XT:A:U' not in read_elems or 'XT:A:U' not in mate_elems:   #both read and it's mate need to be mapped uniquely
    #    continue
    read_chr = read_elems[2]
    read_start = int(read_elems[3])
    read_cigar = read_elems[5]
    if len(read_cigar.split('M')) != 2:     #we want perfect matches only..cigar= <someInt>M
        continue
    read_len = int(read_cigar.split('M')[0])
    mate_chr = mate_elems[2]
    mate_start = int(mate_elems[3])
    mate_cigar = mate_elems[5]
    if len(mate_cigar.split('M')) != 2:     #we want perfect matches only..cigar= <someInt>M
        continue
    mate_len = int(mate_cigar.split('M')[0])
    if read_chr != mate_chr:            # check that they were mapped to the same chromosome
        continue
    if abs(read_start - mate_start) > (maxoriginalreadlength+maxTRlength):
        continue
    if read_start < mate_start:
        pre_s = read_start-1
        pre_e = read_start-1+read_len
        tr_s = read_start-1+read_len
        tr_e = mate_start-1
        suf_s = mate_start-1
        suf_e = mate_start-1+mate_len
    else:
        pre_s = mate_start-1
        pre_e = mate_start-1+mate_len
        tr_s = mate_start-1+mate_len
        tr_e = read_start-1
        suf_s = read_start-1
        suf_e = read_start-1+read_len
    tr_len = abs(tr_e - tr_s)
    if tr_len > maxTRlength:
        continue
    if pre_e >= suf_s:  #overlapping prefix and suffix
        continue
    tr_ref_seq = twobitfile[read_chr][tr_s:tr_e]
    ##print >>fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(read_name,read_chr,pre_s,pre_e,tr_s,tr_e,suf_s,suf_e,tr_len,tr_ref_seq)
    fout.writelines('\t'.join(map(str,[read_name,read_chr,pre_s,pre_e,tr_s,tr_e,suf_s,suf_e,tr_len,tr_ref_seq]))+'\n')

print  "Skipped %d unpaired reads" %(skipped)
