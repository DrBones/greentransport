# -*- coding: utf-8 -*-
"""
Import of geometry to build Hamiltonian

The X-axis is the row Y-Axis the column, X=Down; Y=Right
"""
import numpy
from PIL import Image
import matplotlib.pyplot as plt
import time
from itertools import tee,izip_longest
from scipy.sparse import lil_matrix,eye
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
from numpy.random import rand

#from itertools import islice
def neighbour(iterable):
    iterator = iter(iterable)
    prev = None
    item = iterator.next()
    for next in iterator:
        yield (prev,item,next)
        prev = item
        item = next
    yield (prev,item,None)

class Parameters:
    """Contains all nedded physical constants and parameters"""

    atlas = 'cross10x10onblack.bmp'
    Confinement = 0
    q = 1.6e-19
    hbar = 1.0545e-34/q                         #eV.s
    a = 2e-10                                  #mesh distance in meter
    m0 = 9.1e-31;                              #kg
    m0 = 0.510e6/((3e8)**2)                     #restmass of electron in eV
    mass = m0                                  #effective mass in eV
    t = (hbar**2)/(2*mass*(a**2))
    mu = 0.25                                    #electro-chemical potential in eV
    kT = 0.0025                                  #Temperature * k_boltzmann in eV, 0.0025ev~30K
    lambdaf = 10


class Device:
    """Contains all methods relating to the device in particulr and
    implicit like the defining geometry and the building of the
    Hamiltonian HD"""


    def read_geometry(self, atlas):
        Img = Image.open(atlas)
        Arr = numpy.asarray(Img) 
        self.Contact = []
        self.Contact.append(numpy.array(numpy.where(Arr == 149)))
        self.Contact.append(numpy.array(numpy.where(Arr == 179)))
        self.Contact.append(numpy.array(numpy.where(Arr == 209)))
        self.Contact.append(numpy.array(numpy.where(Arr == 239)))
        self.Conductor = numpy.transpose(numpy.where(Arr > 0))
        return Contact,Conductor

     
    def pairwise(self, seq):
        a,b = tee(seq)
        b.next()
        return izip_longest(a,b)
        
    def compose_geo(self, Cond):
        Nodes = {}
        Count = 0
        for item,next_item in self.pairwise(Cond):
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
     
    def build_HD(self, Nodes, Potential=Parameters.Confinement,t=Parameters.t):
        HD = lil_matrix((len(Nodes),len(Nodes)),dtype=numpy.complex128)   
        HD.setdiag([Potential+4*t]*len(Nodes))   
        for item in Nodes:
            if Nodes[item][1] != None:
                HD[Nodes[item][0],Nodes[item][1]] = -numpy.exp(1j*Nodes[item][0])*t
            else: pass
            if Nodes[item][2] == None: continue
            HD[Nodes[item][0],Nodes[item][2]] = -t  
        #HD = HD + HD.conjugate().T
        HD = (HD.tocsr() + HD.tocsr().conjugate().T).tolil()# might be faster
        return HD

def build_EF(Nodes, lambdaf=Parameters.lambdaf , t=Parameters.t):
    EF = 2*t*(numpy.cos(numpy.pi/lambdaf))*eye(len(Nodes),len(Nodes),dtype=numpy.complex128, format='lil')
    return EF

class Contact:
    """Defines all the methods and attributes needed to evaluate the
    effects of Contacts via the self-energy"""

    def build_SIGMA(self, Nodes, contact):
        SIGMA = lil_matrix((len(Nodes), len(Nodes)), dtype=numpy.complex128)
        for ind_contact in contact:
            for contact_node in range(ind_contact.shape[1]):
              SIGMA[Nodes[tuple(ind_contact.T[contact_node])][0],
                    Nodes[tuple(ind_contact.T[contact_node])][0]] = 22
        return SIGMA

    def transverse_mode(self, contact):
        """calculates the value of the transverse mode's square at the point
        of the contact_node, essentially sin(pi)sin(pi) multiplied with
        appropriate amplitude and constats """
        pass

timeittook=time.time()   
#Cont,Cond = read_geometry(atlas)
#Nodes = compose_geo(Cond)
#HD = build_HD(compose_geo(read_geometry(atlas)[1]),0,t)
#EF = build_EF(lambdaf,t,Nodes)
#S = build_SIGMA(Nodes,Cont)
print time.time()-timeittook
#plt.imshow(numpy.asarray(HD.todense()).real)
#plt.show()

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
