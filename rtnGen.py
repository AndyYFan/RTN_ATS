import numpy as np
from pylab import *

def rtnGen(timelist,kai):
    seq=zeros(101)

    sign=np.random.randint(0,2) #generate random 0,1

    tj=random.expovariate(kai)
    index=0

    for i in range(5):
        if tj<timelist[i]:
            sign=~sign+2
            tj=tj+random.expovariate(kai)
        seq[index]=sign
        index=index+1
    return seq
