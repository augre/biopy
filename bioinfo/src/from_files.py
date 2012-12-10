import itertools
com=itertools.product(range(0,4), repeat=3)
com_list=list(com)
print com_list[42]

def most_common(lst):
    return max(set(lst), key=lst.count)

f=open('result','r')

number_list=[]
for line in f:
    line =line[1:(len(line)-2)]
    number_list.append(line.split(","))

max_pos=[]
max_list=[]
for i in number_list:
    max_pos.append(i.index(max(i)))
    max_list.append(max(i))
    i[i.index(max(i))]=0
print max_list
print max(max_pos)

def twodstrtoint(twodstrlist):
    int_list=[]
    for i in twodstrlist:
        int_list.append(map(int, i))
    return int_list

int_list=twodstrtoint(number_list)


max_pos2=[]
max_list2=[]
for i in int_list:
    max_pos2.append(i.index(max(i)))
    max_list2.append(max(i))
    i[i.index(max(i))]=0
print max_list2
print max(max_pos2)


final=[]
for item in number_list:
    final+=item


def str_list_to_int(str_list):
    int_list=[]
    for item in str_list:
        int_list.append(int(item))
    return int_list

int_list= str_list_to_int(final)




# -*- coding: cp1250 -*-
import matplotlib.pyplot as plt
x=range(0,len(max_list))

def plot():
    plt.plot(int_list, 'ro')
    plt.xlabel('window places')
    plt.ylabel('# of appearances')
    plt.title('Carsonella ruddii blue:AAA green:TTT in 160 windows')
    plt.show()
#plot(x,int_list)
plot()