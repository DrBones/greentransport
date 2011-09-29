class World:

    defaultatlas = 'canvas/wire70x30spinorbit.bmp'

    def __init__(self, atlas=defaultatlas):
        self.atlas = atlas
        self.q = 1.6e-19 #Coulomb
        self.hbar = 6.58211928e-16 #eV * s
        self.c = 299792458
        self.m0 = 0.510e6/(self.c**2) #eV*s**2/m**2
        self.eps0 = 8.854e-12 # Vacuum permittivity C/V*m
        self.kb = 8.6173324e-5 # ev /K

        self.__read_geometry()
        #self.__compose_nodes()
        self.__blocksizes_from_coords()

    def __read_geometry(self):
        """reads in bmp and generates contacts and other
        objects of interest"""
        from PIL import Image
        from scipy import where, asarray, array, transpose,logical_or
        from aux import Contact
        img = Image.open(self.atlas)
        arr = asarray(img)
        contacts = []
        shades = [(109,119), (139,149), (169,179), (199,209)]
        contact_index = 0
        leads = []
        for shade in shades:
            a = Contact(where(arr == shade[1]))
            if a.shape[1] == 0: continue
            lead = Contact(logical_or(arr == shade[0],arr ==shade[1]).nonzero())
            lead.index = contact_index
            leads.append(lead)
            a.index = contact_index
            #from pudb import set_trace; set_trace()
            try:
                if any(arr[a[0][0]+1,:]==239):
                    a.SO = True
                else:
                    a.SO = False
            except IndexError:
                if any(arr[a[0][0]-1,:] == 239):
                    a.SO = True
                else:
                    a.SO = False
            contacts.append(a)
            contact_index +=1
        self.raw_coords = where(arr > 0)
        self.active_coords = transpose(self.raw_coords)
        self.wafer = arr
        self.canvas = arr.shape
        self.contacts = contacts
        self.leads = leads

    def __compose_nodes(self):
        from collections import OrderedDict
        from aux import pairwise
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
