from orf import *
from math import log
from training import training_noncode, training_orfs

class ZEROTH:
    """
    Input is a list of character lists, intergenetic region or ORFs
    """
    def __init__(self, seq):
        self.seq,self.res,self.nonres=seq,[],[]
        self.freqs=[[0,0,0],     #P1,2,3(A)
                    [0,0,0],     #P1,2,3(C)
                    [0,0,0],     #P1,2,3(T)
                    [0,0,0]]     #P1,2,3(G)
        self.nonfreqs=[0,0,0,0]  #P4(A,C,T,G)

    def orf_fillup(self):
        """
        Fills up the index 1,2,3 variables with the number of each letter at each position
        """
        for s in self.seq:            
            i=0
            self.freqs=[[0,0,0],     #P1,2,3(A)
                        [0,0,0],     #P1,2,3(C)
                        [0,0,0],     #P1,2,3(T)
                        [0,0,0]]    #P1,2,3(G)
            while i < (len(s)):
                col=(i%12)%3
                if s[i]=='A':
                    self.freqs[0][col]=self.freqs[0][col]+1
                elif s[i]=='C':
                    self.freqs[1][col]=self.freqs[1][col]+1
                elif s[i]=='T':
                    self.freqs[2][col]=self.freqs[2][col]+1
                elif s[i]=='G':
                    self.freqs[3][col]=self.freqs[3][col]+1
                i=i+1
            self.res.append((len(s),self.freqs)) #length of the orf and numbers of appearances
            
    def noncoding_fillup(self):
        """
        Fills up the index 4 variables with the number of letters in the noncoding regions
        """
        self.indicator=1
        for s in self.seq:
            self.nonfreqs=[0,0,0,0]  #P4(A,C,T,G)
            for nuc in s:
                if nuc=='A': self.nonfreqs[0]=self.nonfreqs[0]+1
                elif nuc=='C': self.nonfreqs[1]=self.nonfreqs[1]+1
                elif nuc=='T': self.nonfreqs[2]=self.nonfreqs[2]+1
                elif nuc=='G': self.nonfreqs[3]=self.nonfreqs[3]+1
            self.nonres.append((len(s),self.nonfreqs))
    
    def log_up_p(self):
        """
        Logs up the probabilities, returns a list with the scores with matching indexes.
        """
        self.log_ups=[]
        i,j,p=0,0,0
        for s in self.res:
            i,p=0,0
            while i<4:
                while j<3:
                    if s[0]!=0 and float(s[1][i][j])!=0:
                        p+=log(float(s[1][i][j])/s[0])
                    else:
                        p=0
                    j+=1
                j=0
                i+=1
            self.log_ups.append(p)
    
    def pure_data(self):
        """
        No index matching after this, takes out the 0s, from the orf probabilities
        """
        self.pure_logs=[]
        for s in self.log_ups:
            if s!=0:
                self.pure_logs.append(s)
        
    def nonlog_up_p(self):
        """
        Logs up the probabilities, returns a list with the scores with matching indexes.
        """
        self.nonlog_ups=[]
        for s in self.nonres:
            i,p=0,0
            while i < 4:
                if s[0]!=0 and float(s[1][i])!=0:
                    p+=log(float(s[1][i])/s[0])
                else:
                    p=0
                i+=1
            self.nonlog_ups.append(p)
            
    def nonpure_data(self):
        """
        No index matching after this, takes out the 0s, from the noncoding probabilities
        """
        self.nonpure_logs=[]
        for s in self.nonlog_ups:
            if s !=0:
                self.nonpure_logs.append(s)
def print_prop(numlist):
    print "\n"
    print max(numlist)
    print min(numlist)
    print sum(numlist)/len(numlist)

def demo():
#    Reading frame 1
    a1=ORF_FINDER(rf1)
    a1.orf_collector()
#    Fill up the counters for the
#    orfs of rf1
    z1o=ZEROTH(a1.orfs)
    z1o.orf_fillup()
    z1o.log_up_p()
    z1o.pure_data()
    print z1o.res[0]
    print "Reading frame 1 open reading frames scores max, min, avg:"
    print_prop(z1o.pure_logs)
    
    #    Reading frame 2
    a2=ORF_FINDER(rf2)
    a2.orf_collector()
#    Fill up the counters for the
#    orfs of rf2
    z2o=ZEROTH(a2.orfs)
    z2o.orf_fillup()
    z2o.log_up_p()
    z2o.pure_data()
    print "\nReading frame 2 open reading frames scores max, min, avg:"
    print_prop(z2o.pure_logs)
    
    #    Reading frame 3
    a3=ORF_FINDER(rf3)
    a3.orf_collector()
#    Fill up the counters for the
#    orfs of rf3
    z3o=ZEROTH(a3.orfs)
    z3o.orf_fillup()
    z3o.log_up_p()
    z3o.pure_data()
    print "\nReading frame 3 open reading frames scores max, min, avg:"
    print_prop(z3o.pure_logs)
    
#    rf1 noncoding
    z1i=ZEROTH(a1.intgens)
    z1i.noncoding_fillup()
    z1i.nonlog_up_p()
    z1i.nonpure_data()
    print "\nReading frame 1 noncoding regions scores max, min, avg:"
    print_prop(z1i.nonpure_logs)
    
#    rf2 noncoding
    z2i=ZEROTH(a2.intgens)
    z2i.noncoding_fillup()
    z2i.nonlog_up_p()
    z2i.nonpure_data()
    print "\nReading frame 2 noncoding regions scores max, min, avg:"
    print_prop(z2i.nonpure_logs)

#    rf1 noncoding
    z3i=ZEROTH(a3.intgens)
    z3i.noncoding_fillup()
    z3i.nonlog_up_p()
    z3i.nonpure_data()
    print "\nReading frame 3 noncoding regions scores max, min, avg:"
    print_prop(z3i.nonpure_logs)    

#    Get the training set data
    tr=ZEROTH(training_orfs)
    tr.orf_fillup()
    tr.log_up_p()
    print "\nTraining set data orfs scores max, min, avg:"
    print_prop(tr.log_ups)
    
    nontr=ZEROTH(training_noncode)
    nontr.noncoding_fillup()
    nontr.nonlog_up_p()
    nontr.nonpure_data()
    print "\nTraining set data noncoding scores max, min, avg:"
    print_prop(nontr.nonpure_logs)


    
    
if __name__ == "__main__": 
    demo()
                
                
                
                    