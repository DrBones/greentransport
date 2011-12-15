class UserParams():
    """
    a [m]
    alpha [eV m]
    effmassfactor [1]
    Temp [Kelvin]
    lambdaf [a]
    BField [?]
    potential_drop [eV?]
    """
    def __init__(self,
                 a = 1e-9,
                 alpha = 20e-12,
                 band_bottom = 0,
                 effmassfactor = 0.026,
                 #effmassfactor = 0.068,
                 Temp = 4,
                 #lambdaf = 35,
                 lambdaf = 35,
                 BField = 0,
                 potential_drop = None):
        if potential_drop is None:
            self.potential_drop = [0,0]
        self.a = a
        self.alpha = alpha
        self.band_bottom = band_bottom
        self.effmassfactor = effmassfactor
        self.Temp = Temp
        self.lambdaf = lambdaf
        self.BField = BField

upar = UserParams()
