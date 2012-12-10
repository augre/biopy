class ORF_FINDER:
    def __init__(self, seq, pos=0):
        self.seq=seq
        self.pos=pos
        self.length=len(seq)
        
    def find_start_codon(self):
        """
        Finds a start codon by not overlaping
        forward from the position pos
        """
        self.initpos=self.pos
        while self.pos<(len(self.seq)-2):
            triplet=str(self.seq[self.pos])+str(self.seq[self.pos+1])+str(self.seq[self.pos+2])
            if triplet=='ATG':
                self.intgen_end=self.pos
                return True
            self.pos=self.pos+3
            
    def find_stop_codon(self):
        """
        Finds a stop codon by not overlaping
        forward from the position pos
        """
        while self.pos<(len(self.seq)-2):
            triplet=str(self.seq[self.pos])+str(self.seq[self.pos+1])+str(self.seq[self.pos+2])
            if triplet=='TAA' or triplet=='TAG' or triplet=='TGA':
                return True
            self.pos=self.pos+3
            
    
    def find_orf(self):
        """
        Finds a start codon, looks for the nearest stop codon
        and returns the sequence in between
        """
        if self.find_start_codon():
            start=self.pos
            if self.find_stop_codon():
                stop=self.pos
                if len(self.seq[start+3:stop])>0:
                    self.orf=self.seq[start+3:stop]
                return True
            
    def find_intgen(self):
        """
        Finds the intergenetic regions
        """
        end=self.intgen_end
        start=self.initpos+3
        self.intergen=self.seq[start:end]
        
    def orf_collector(self):
        """
        Puts every ORF from the sequence to a list
        """
        self.orfs,self.intgens=[],[]
        while 1:
            if self.find_orf():
                self.orfs.append(self.orf)
                self.find_intgen()
                if len(self.intergen)>1:
                    self.intgens.append(self.intergen)
            else:
                return self.orfs
            
    def triplet_numbers(self):
        """
        Makes a list containing the number of the triplets for each orf
        """
        num=[]
        for orf in self.orfs:
            num.append(len(orf)/3)
        return num
    
    def avg_triplet_num(self):
        average=float(sum(self.triplet_numbers()))/len(self.triplet_numbers())
        return average
        
from Bio import SeqIO
#from Bio.Seq import Seq
record = SeqIO.read("good.fna", "fasta")    
my_seq=record.seq
#The DNA and its reverse complement
#in two lists
my_seq_reverse_comp=list(my_seq.complement()[::-1])
my_seq=list(record.seq)

#elkeszitem a reading freameket
rf1=my_seq
rf2=my_seq[1:]
rf3=my_seq[2:] 

   


def demo():
    print '\nFirst reading frame:'
    a=ORF_FINDER(rf1)
    orfs=a.orf_collector()
    print 'number of orfs:' 
    print len(orfs)
    print 'number of intergenetic regions:'
    print len(a.intgens)


#   number of triplets in orfs
    a.triplet_numbers()
    print 'average triplet number in orfs:'
    print a.avg_triplet_num() 
    
    print '\nSecond reading frame:'
    b=ORF_FINDER(rf2, 0)
    orfs=b.orf_collector()
    print 'number of orfs:' 
    print len(orfs)
    print 'number of intergenetic regions:'
    print len(b.intgens)
#    number of triplets in orfs
    b.triplet_numbers()
    print 'average triplet number in orfs:'
    print b.avg_triplet_num() 

    print '\nThird reading frame:'
    c=ORF_FINDER(rf3, 0)
    orfs=c.orf_collector()
    print 'number of orfs:' 
    print len(orfs)
    print 'number of intergenetic regions:'
    print len(c.intgens)
#    number of triplets in orfs
    c.triplet_numbers()
    print 'average triplet number in orfs:'
    print c.avg_triplet_num()
    
if __name__ == "__main__": 
    demo()        
            