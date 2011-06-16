class World:

    defaultatlas = 'wirecentralcontact100x50nopad.bmp'

    def __init__(self, atlas=defaultatlas):
        self.atlas = atlas
        self.q = 1.6e-19 #Coulombs
        self.hbar = 1.0545e-34/self.q
        self.m0 = 0.510e6/((3e8)**2)
        self.eps0 = 8.854e-12 # Vacuum permittivity F m**-1

        self.__read_geometry()
        self.__compose_nodes()
        self.__blocksizes_from_coords()

    def __read_geometry(self):
        from PIL import Image
        from scipy import where, asarray, array, transpose
        img = Image.open(self.atlas)
        arr = asarray(img)
        contacts = []
        contact_shades = [149, 179, 209, 239]
        for shade in contact_shades:
            a = array(where(arr == shade))
            if a.shape[1] == 0: continue
            contacts.append(a)
        active_coords = transpose(where(arr > 0))
        self.wafer = arr
        self.canvas = arr.shape
        self.contacts = contacts
        self.active_coords = active_coords

    def __compose_nodes(self):
        from collections import OrderedDict
        from customiterators import pairwise
        nodes = OrderedDict()
        count = 0
        for item,next_item in pairwise(self.active_coords):
            try:
                if item[1] == next_item[1]-1:
                    nodes[tuple(item)] = [count,count+1,None]
                else:
                    nodes[tuple(item)] = [count,None,None]
            except TypeError:
                nodes[tuple(item)] = [count,None,None]
            if item[0]>self.active_coords[0][0]:
                try:
                    nodes[tuple(item - [1,0])][2] = count
                except KeyError:
                    pass
            count +=1
        self.nodes = nodes

    def __blocksizes_from_coords(self):
        block_size = 1
        block_sizes = []
        for i in range(1,len(self.active_coords)):
            if self.active_coords[i][0] == self.active_coords[i-1][0]:
                block_size +=1
            else:
                block_sizes.append(block_size)
                block_size = 1
        block_sizes.append(block_size)
        self.block_sizes = block_sizes
