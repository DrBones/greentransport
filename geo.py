# -*- coding: utf-8 -*-
"""
Import of geometry to build Hamiltonian

The X-axis is the row Y-Axis the column, X=Down; Y=Right
"""
import numpy
from PIL import Image
import matplotlib.pyplot as plt
import time
from itertools import tee,chain,izip,izip_longest
from scipy.sparse import lil_matrix,eye
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
from numpy.random import rand

#from itertools import islice
atlas = 'cross10x10onblack.bmp'
Confinement = 0
q = 1.6e-19
hbar = 1.0545e-34/q                         #eV.s
a = 2e-10                                  #mesh distance in meter
#m0 = 9.1e-31;                              #kg
m0 = 0.510e6/((3e8)**2)                     #restmass of electron in eV
mass = m0                                  #effective mass in eV
t = (hbar**2)/(2*mass*(a**2))
mu = 0.25                                    #electro-chemical potential in eV
kT = 0.0025                                  #Temperature * k_boltzmann in eV, 0.0025ev~30K
lambdaf = 10

def read_geometry(atlas):
    Img = Image.open(atlas)
    Arr = numpy.asarray(Img) 
    Contact = numpy.transpose(numpy.where(Arr == 178))
    Conductor = numpy.transpose(numpy.where(Arr > 0))
    return Contact,Conductor

def neighbour(iterable):
    iterator = iter(iterable)
    prev = None
    item = iterator.next()
    for next in iterator:
        yield (prev,item,next)
        prev = item
        item = next
    yield (prev,item,None)
 
def pairwise(seq):
    a,b = tee(seq)
    b.next()
    return izip_longest(a,b)
    
def compose_geo(Cond):
    Nodes = {}
    Count = 0
    for item,next_item in pairwise(Cond):
        try:
            if item[1] == next_item[1]-1:
                Nodes[tuple(item)] = [Count,Count+1,None]
            else:
                Nodes[tuple(item)] = [Count,None,None]
        except TypeError:
            Nodes[tuple(item)] = [Count,None,None]                
        if item[0]>1:
            try:
                Nodes[tuple(item - [1,0])][2] = Count
            except KeyError:
                pass
        Count +=1 
    return Nodes
 
def build_HD(Nodes,Confinement,t):
    HD = lil_matrix((len(Nodes),len(Nodes)),dtype=numpy.complex128)   
    HD.setdiag([Confinement+4*t]*len(Nodes))   
    for item in Nodes:
        if Nodes[item][1] != None:
            HD[Nodes[item][0],Nodes[item][1]] = -numpy.exp(1j*Nodes[item][0])*t
        else: pass
        if Nodes[item][2] == None: continue
        HD[Nodes[item][0],Nodes[item][2]] = -t  
    HD = HD + HD.conjugate().T
    #    HD = (HD.tocsr() + HD.tocsr().conjugate().T).tolil might be faster
    return HD

def build_EF(lambdaf,t,Nodes):
    EF = 2*t*(numpy.cos(numpy.pi/lambdaf))*eye(len(Nodes),len(Nodes),dtype=numpy.complex128, format='lil')
    return EF

def build_SIGMA(Nodes,Cont):
    SIGMA = lil_matrix((len(Nodes),len(Nodes)),dtype=numpy.complex128)
    
    return SIGMA
    
Cont,Cond = read_geometry(atlas)
Nodes = compose_geo(Cond)
HD = build_HD(compose_geo(read_geometry(atlas)[1]),0,t)
EF = build_EF(lambdaf,t,Nodes)
plt.imshow(numpy.asarray(HD.todense()).real)
plt.show()
#t=time.time()
#Nodes = compose_geo(Cond)
#print "The Process took", time.time()-t, "seconds"
#Nodes2 = compose_geo2(Cond)
#def main():
#    Cont,Cond = read_geometry(atlas)            
#    Nodes = compose_geo(Cond)  
#
#if __name__ == '__main__':
#    main()
#    

