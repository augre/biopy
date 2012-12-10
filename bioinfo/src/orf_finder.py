#! /usr/bin/env python
# -*- coding: Utf-8 -*-

f=open('/home/augre/sequence.fna')
sequence=f.read()
#find the biggining of orfs
#print sequence.count('ATG') #non overlapping
class Sequence:
    def __init__(self,seq):
        self.seq=seq
    def start_codons_fwd(self,codon='ATG'):
        "find start codons"
        pos_list,pos=[],0
        while(1):
            pos=self.seq.find(codon,pos)
            if pos==-1: break
            pos_list.append(pos)
            pos=pos+1
        self.s_cods_fwd=pos_list
        return pos_list
    def start_codons_bwd(self,codon='ATG'):
        self.seq=self.seq[::-1]
        lista=self.start_codons_fwd()
        self.s_cods_bwd=lista
        return lista
    def find_orfs_fwd(self):
        
def demo():
    s=Sequence(sequence)
    fwd=s.start_codons_fwd()
    bwd=s.start_codons_bwd()
    print len(bwd)

if __name__=="__main__":
    demo()
