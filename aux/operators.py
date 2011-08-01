def spindens(self,lrgm_out):
    from scipy import split,pi
    Gup, Gdown = split(lrgm_out.reshape(self.wafer.shape[0],self.wafer.shape[1]*2),2,axis=1)
    Sz = self.hbar/(4*pi*1j*self.a**2)*(Gup-Gdown)
    print 'max Spin Split: ', Sz.imag.max()-Sz.imag.min()
#Realteil scheint wahrscheinlicher, imag oszilliert wie bloed
    return Sz.real

def edens(self,lrgm_out):
    from scipy import split,pi,asarray
    if self.multi ==1:
        edensity = lrgm_out.reshape(asarray(self.wafer.shape)+[0,0])*2/(2*pi*self.a**2) #times 2 for spin
    if self.multi ==2:
        Gup, Gdown = split(lrgm_out.reshape(self.wafer.shape[0],self.wafer.shape[1]*2),2,axis=1)
        edensity = 1/(2*pi*self.a**2)*(Gup+Gdown)
    return edensity.real
