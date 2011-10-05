import numpy as np
import scipy.linalg as sl

def heaviside(x):
    from scipy import sign
    heaviside=  ((sign(x)+1)/2)
    return heaviside

def Pot(x,y,factor=1):
    from scipy import cos,pi
    size_y = y.max()
    size_x = x.max()
    V_zero = 0.5;
    Ef = 1;
    delta = 1;
    xoffset = size_x/2
    yoffset = size_y/2
    x = x+xoffset
    y = y-yoffset
    yplus = (1-factor*size_y*cos(2*pi*x/size_x));
    yminus = -(1-factor*size_y*cos(2*pi*x/size_x));
    sum1 = ((y-size_y-yplus)/delta)**2*heaviside(y-yplus-size_y);
    sum2 = ((y+size_y-yminus)/delta)**2*heaviside(-(y-yminus+size_y));
    V = V_zero/2+Ef*(sum1+sum2);
    return V

def sphericalPot(x,y,shift=0,radius=1,scale=1):
    from scipy import real,sqrt
    size_x = x.max()
    size_y = y.max()
    Vbottom = 0
    x = x-size_x/2
    left_sphere = real(sqrt(radius**2-(x**2+(y-shift)**2)))*heaviside(y-shift)+real(sqrt(radius**2-x**2))*heaviside(-(y-shift))
    right_sphere = real(sqrt(radius**2-(x**2+(y-size_y+shift)**2)))*heaviside(-(y-size_y+shift))+real(sqrt(radius**2-x**2))*heaviside((y-size_y+shift))
    V = Vbottom +scale*(left_sphere+right_sphere)
    return V

def invit(A, eigenvalue, tolerance):
    # solve equation (A-mu*I)w = v**(k-1)
    n = A.shape[0]
    matrix = A-eigenvalue*np.eye(n)
    lup = sl.lu_factor(matrix)
    # initialize random v with norm  = 1
    v = np.random.randn(n)
    v /= sl.norm(v)
    v_old = np.zeros(n)
    while (sl.norm(v-v_old)> tolerance):
         v_old = v
         w  = sl.lu_solve(lup, v_old)
         v =  w/sl.norm(w)
    return v

def eigenvector_from_eigenvalue(A, eigenvalue):
    # solve equation (A-mu*I)w = v**(k-1) once
    n = A.shape[0]
    matrix = A-eigenvalue*np.eye(n)
    lup = sl.lu_factor(matrix)
    # initialize random v with norm  = 1
    v_init = np.random.randn(n)
    v_init /= sl.norm(v_init)
    w  = sl.lu_solve(lup, v_init)
    vector =  w/sl.norm(w)
    return vector

def all_elements_are_unique(array):
    for i in range(len(array)):
        for j in set(range(len(array)))-set([i]):
            if abs(array[i]-array[j]) < 1e-3: return False
            else: return True

def nodenames_of_contact(lead_graph_list, contact):
    pass
