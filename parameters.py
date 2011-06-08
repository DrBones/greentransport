import scipy
class Parameters:
    """Contains all needed physical constants and parameters"""

    atlas = 'wirecentralcontact100x50.bmp'
    potential_drop = [-0.01, 0.01] #Potential Drop over legth of device
    q = 1.6e-19
    hbar = 1.0545e-34/q                                                 #eV.s
    a = 2e-10                                                           #mesh distance in meter
    #m0 = 9.1e-31;                                                       #kg
    m0 = 0.510e6/((3e8)**2)                                             #restmass of electron in eV
    mass = 0.25*m0                                                      #effective mass in eV
    t = (hbar**2)/(2*mass*(a**2))
    kT = 0.025                                                          #Temperature * k_boltzmann in eV, 0.0025ev~30K
    lambdaf = 10
    BField = 0
    zplus = 1j*1e-12
    Egrid = scipy.linspace(0.4,1.0,100)+zplus
    #Efermi = 2*t*(scipy.cos(scipy.pi/lambdaf))
    Efermi = 0.1
    dE = Egrid[1].real-Egrid[0].real
    mu_l = Efermi + (potential_drop[1] - potential_drop[0])/2                                                      #electro-chemical potential in eV
    mu_r = Efermi - (potential_drop[1] - potential_drop[0])/2
