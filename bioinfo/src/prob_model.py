from orf import *
from math import log

class ZEROTH:
    """
    Input is a list of character lists, intergenetic region or ORFs
    """
    def __init__(self, seq):
        self.seq=seq
        self.pa1,self.pa2,self.pa3=0,0,0
        self.pc1,self.pc2,self.pc3=0,0,0
        self.pt1,self.pt2,self.pt3=0,0,0
        self.pg1,self.pg2,self.pg3=0,0,0
        self.pa4,self.pc4,self.pt4,self.pg4=0,0,0,0
        
    def orf_fillup(self):
        """
        Fills up the index 1,2,3 variables with the number of each letter at each position
        """
        p=0
        for s in self.seq:            
            i=0
            self.pa1,self.pa2,self.pa3=0,0,0
            self.pc1,self.pc2,self.pc3=0,0,0
            self.pt1,self.pt2,self.pt3=0,0,0
            self.pg1,self.pg2,self.pg3=0,0,0
            while i < (len(s)-2):
                t1,t2,t3=s[i],s[i+1],s[i+2]
                if t1=='A':self.pa1=self.pa1+1
                if t2=='A':self.pa2=self.pa2+1
                if t3=='A':self.pa3=self.pa3+1
                
                if t1=='C':self.pc1=self.pc1+1
                if t2=='C':self.pc2=self.pc2+1
                if t3=='C':self.pc3=self.pc3+1
                
                if t1=='T':self.pt1=self.pt1+1
                if t2=='T':self.pt2=self.pt2+1
                if t3=='T':self.pt3=self.pt3+1
                
                if t1=='G':self.pg1=self.pg1+1
                if t2=='G':self.pg2=self.pg2+1
                if t3=='G':self.pg3=self.pg3+1
                i=i+3
            num=(len(s)/3)
            if self.pa1>0 and self.pa2>0 and self.pa3>0 and self.pc1>0 and self.pc2>0 and self.pc3>0 and self.pt1>0 and self.pt2>0 and self.pt3>0 and self.pg1>0 and self.pg2>0 and self.pg3>0:
                p=p+log(float(self.pa1)/num)+log(float(self.pa2)/num)+log(float(self.pa3)/num)+log(float(self.pc1)/num)+log(float(self.pc2)/num)+log(float(self.pc3)/num)+log(float(self.pt1)/num)+log(float(self.pt2)/num)+log(float(self.pt3)/num)+log(float(self.pg1)/num)+log(float(self.pg2)/num)+log(float(self.pg3)/num)
        self.orf_p=p
        
    def intgen_fillup(self):
        """
        Fills up the index 4 variables with the number of letters in the intergenetic regions
        """
        p=0
        for s in self.seq:
            self.pa4,self.pc4,self.pt4,self.pg4=0,0,0,0
            for b in s:
                if b=='A':self.pa4=self.pa4+1
                
                if b=='C':self.pc4=self.pc4+1
                
                if b=='T':self.pt4=self.pt4+1
                
                if b=='G':self.pg4=self.pg4+1
            num=len(s)
            if self.pa4>0 and self.pc4>0 and self.pt4>0 and self.pg4>0:
                p=p+log(float(self.pa4)/num) +log(float(self.pc4)/num)+log(float(self.pt4)/num)+log(float(self.pg4)/num)
            self.intgen_p=p
        
def demo():
    a1=ORF_FINDER(rf1)
    a1.orf_collector()
#    Fill up the counters for the
#    orfs of rf1
    z1o=ZEROTH(a1.orfs)
    z1o.orf_fillup()
    print '\nFirst reading frame ORFs log(P1(A))+.....:'
    print z1o.orf_p
#    fill up the counters for the
#    inter genetic regions of rf1
    z1i=ZEROTH(a1.intgens)
    z1i.intgen_fillup()
    print 'First reading frame inter genetic log(P4(A))+.....:'
    print z1i.intgen_p
    
    a1=ORF_FINDER(rf2)
    a1.orf_collector()
#    Fill up the counters for the
#    orfs of rf2
    z1o=ZEROTH(a1.orfs)
    z1o.orf_fillup()
    print '\nSecond reading frame ORFs log(P1(A))+.....:'
    print z1o.orf_p
#    fill up the counters for the
#    inter genetic regions of rf1
    z1i=ZEROTH(a1.intgens)
    z1i.intgen_fillup()
    print 'Second reading frame inter genetic log(P4(A))+.....:'
    print z1i.intgen_p
    
    a1=ORF_FINDER(rf3)
    a1.orf_collector()
#    Fill up the counters for the
#    orfs of rf3
    z1o=ZEROTH(a1.orfs)
    z1o.orf_fillup()
    print '\nThird reading frame ORFs log(P1(A))+.....:'
    print z1o.orf_p
#    fill up the counters for the
#    inter genetic regions of rf1
    z1i=ZEROTH(a1.intgens)
    z1i.intgen_fillup()
    print 'Third reading frame inter genetic log(P4(A))+.....:'
    print z1i.intgen_p
    
    
if __name__ == "__main__": 
    demo()
        
    
            