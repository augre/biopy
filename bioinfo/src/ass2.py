import itertools
from Bio import SeqIO
#from Bio.Seq import Seq
record = SeqIO.read("francisella.fasta", "fasta")    
my_seq=record.seq

rf1= str(my_seq)

com=itertools.product(range(0,4), repeat=3)
com_list=list(com)




def seq_to_numbers(my_seq):
    output=[]
    for letter in my_seq:
        if letter=='A':
            output.append(0)
        elif letter=='C':
            output.append(1)
        elif letter=='T':
            output.append(2)
        elif letter=='G':
            output.append(3)
    return output


def window_size(my_seq):
    window_size=len(my_seq)/2.0
    if window_size%3!=0:
        window_size-=window_size%3
    return int(window_size)


def extend_the_seq(my_seq,end_pos):
    return my_seq+my_seq[0:end_pos]

        
        
def count_appearances(window,com):
    appear_list=[0]*64
    i=0
    while i<(len(window)-2):
        triplet=(window[i],window[i+1],window[i+2])
        pos=com.index(triplet)
        appear_list[pos]+=1
        i+=3
    return appear_list
    
def tally_the_seq(my_seq, jump_number=997):
    """
    """
    result=[]
    w_size=window_size(my_seq)
    start_pos=0
    end_pos=start_pos+w_size
    new_seq=extend_the_seq(list(my_seq),end_pos)
    while end_pos<=len(new_seq):
        window=new_seq[start_pos:end_pos]
        result.append((start_pos,count_appearances(window,com_list)))
        start_pos+=jump_number
        end_pos+=jump_number
    return result



seq_sample='ACGTTGCAATGCCAGT'

num_seq=seq_to_numbers(rf1)
#print len(num_seq),"==?",len(rf1)
res=tally_the_seq(num_seq)

f=open('result', 'w')
for item in res:
    f.write(str("%s\n" % ''.join(str(item[1]))))
f.close()


