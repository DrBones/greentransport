def rrgm(Ablock):
    """ Performs recursive algorithm (Svizhenko et. al) to calculate
    the retarded green's function, uses views on A, i.e. the block
    matrices of A, by Ablock """
    from scipy import array, hstack
    from collections import OrderedDict
    number_of_blocks = len(Ablock)
    grl = OrderedDict()
    grl[0,0] = Ablock[0,0].todense().I
    prev_greensfnc = grl[0,0]
    for i in range(1,number_of_blocks):
        prev_greensfnc = (Ablock[i,i]-Ablock[i, i-1] * prev_greensfnc * Ablock[i-1,i]).I
        grl[i,i] = prev_greensfnc
        grl[i,0] = grl[i,i]*Ablock[i,i-1]*grl[i-1,0]
        #Gr only needed in final calculations of lesser Gr, possible speed gain
    Gr = {}
    Gr[number_of_blocks-1,number_of_blocks-1] = grl[(number_of_blocks-1,number_of_blocks-1)]
    rev_iterator = reversed(range(1,number_of_blocks))
    for i in rev_iterator:
        Gr[i, i-1] = Gr[i,i] * Ablock[i,i-1] * grl[i-1,i-1]
        Gr[i-1, i] = grl[i-1,i-1] * Ablock[i-1,i] * Gr[i,i]
        Gr[i-1, i-1] = grl[i-1,i-1]+grl[i-1,i-1] * Ablock[i-1,i] * Gr[i, i-1]
    green_diag = array([])
    for i in range(number_of_blocks):
        diag = array(Gr[i,i].diagonal()).reshape(-1)
        green_diag = hstack((green_diag, diag))
    return green_diag, grl, Gr

