def spy(mat):
    from summon import matrix
    from itertools import izip
    #from scipy.io import mmwrite
    #from os import system
    #mmwrite('tempmatrix', mat.real)
    #system("ruby zerobased.rb")
    #matrix.open_matrix('newmatrix.txt', m, format='imat')
    m = matrix.Matrix()
    mat = mat.real.tocoo()
    nrows, ncols = mat.shape
    nnz = mat.getnnz()
    imat = izip(mat.row, mat.col, mat.data)
    matrix.load_matrix(nrows, ncols, nnz, imat, m)
    viewer = matrix.MatrixViewer(m, title="Sparsity")
    viewer.show()
    viewer.draw_border(30,30)
    viewer.draw_border(60,60)
