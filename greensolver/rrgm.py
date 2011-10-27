def rrgm(Ablock):
    """ Performs recursive algorithm (Svizhenko et. al) to calculate
    the retarded green's function, uses views on A, i.e. the block
    matrices of A, by Ablock """
    from numpy import isfinite,nan_to_num
    from scipy import array, hstack
    from collections import OrderedDict
    number_of_blocks = len(Ablock)
    grl = OrderedDict()
    grl[0,0] = Ablock[0,0].todense().I
    while not isfinite(grl[0,0]).all():
        print 'grl contains NaN'
        grl[0,0] = Ablock[0,0].todense().I
    prev_greensfnc = grl[0,0]
    for i in range(1,number_of_blocks):
        prev_greensfnc = Ablock[i,i]-Ablock[i, i-1] * prev_greensfnc * Ablock[i-1,i]
        ind = 0
        while not isfinite(prev_greensfnc).all():
            print ind,'prev_greensfnc contains NaNs at block',i
            prev_greensfnc = Ablock[i,i]-Ablock[i, i-1] * prev_greensfnc * Ablock[i-1,i]
            ind+=1
        prev_greensfnc_tmp = prev_greensfnc.I
        check_if_equal = nan_to_num(prev_greensfnc_tmp)
        while not isfinite(prev_greensfnc_tmp).all():
            print ind,'prev_greensfnc inverse contains NaNs at block',i
            prev_greensfnc_tmp = prev_greensfnc.I
            if (abs((prev_greensfnc_tmp-check_if_equal)/prev_greensfnc_tmp) < 1e-6).all():
                print "Don't do the while loop, use nan_to_num, it's the same"
        prev_greensfnc = prev_greensfnc_tmp
        grl[i,i] = prev_greensfnc
        grl[i,0] = -grl[i,i]*Ablock[i,i-1]*grl[i-1,0]
        #Gr only needed in final calculations of lesser Gr, possible speed gain
    Gr = {}
    Gr[number_of_blocks-1,number_of_blocks-1] = grl[(number_of_blocks-1,number_of_blocks-1)]
    rev_iterator = reversed(range(1,number_of_blocks))
    for i in rev_iterator:
        Gr[i, i-1] = -Gr[i,i] * Ablock[i,i-1] * grl[i-1,i-1]
        Gr[i-1, i] = -grl[i-1,i-1] * Ablock[i-1,i] * Gr[i,i]
        Gr[i-1, i-1] = grl[i-1,i-1]-grl[i-1,i-1] * Ablock[i-1,i] * Gr[i, i-1]
    green_diag = array([])
    for i in range(number_of_blocks):
        diag = array(Gr[i,i].diagonal()).reshape(-1)
        green_diag = hstack((green_diag, diag))
    return green_diag, grl, Gr

