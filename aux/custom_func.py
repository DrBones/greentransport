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
