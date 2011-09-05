def lrgm(Ablock, sigma_in_l, sigma_in_r):
    """ Performs recursive algorithm (Svizhenko et.al) to calculate
    lesser green's function by the use of the diagonal and offdiagonal
    matrix-elements of the retarded green's function """
    from scipy import array, hstack
    from greensolver import rrgm
    ignored, grl, Gr = rrgm(Ablock)
    number_of_blocks = len(Ablock)
    gll = [grl[0,0] * sigma_in_l * grl[0,0].getH()]
    #len(block)-1 equals N-2 because of pythons way of building the range exclusive the endpoint
    for i in range(1,number_of_blocks-1):
        prev_lesser_greenfnc = grl[i,i] * (Ablock[i,i-1] * gll[i-1] * Ablock[i,i-1].getH()) * grl[i,i].getH()
        gll.append(prev_lesser_greenfnc)
    #N-1 step extra, belongs to for loop
    i +=1
    gll.append(grl[i,i] * (sigma_in_r + Ablock[i,i-1] * gll[i-1] * Ablock[i-1,i].conj()) * grl[i,i].getH())
    Gl = {}
    Gl[number_of_blocks-1,number_of_blocks-1] = gll[-1]
    for i in reversed(range(1,number_of_blocks)):
        Gl[i,i-1] = Gr[i,i] * Ablock[i,i-1] * gll[i-1] - Gl[i,i] * Ablock[i-1,i].getH() * grl[i-1,i-1].getH()
        Gl[i-1,i-1] = (gll[i-1] + grl[i-1,i-1] * (Ablock[i-1,i] * Gl[i,i] * Ablock[i-1,i].getH()) * grl[i-1,i-1].getH() + (gll[i-1] * Ablock[i,i-1].getH() * Gr[i-1,i].getH() + Gr[i-1,i] * Ablock[i,i-1] * gll[i-1]))
    less_green_diag = array([])
    for i in range(number_of_blocks):
        less_diag = array(Gl[i,i].diagonal()).reshape(-1)
        less_green_diag = hstack((less_green_diag, less_diag))
    return less_green_diag, grl
