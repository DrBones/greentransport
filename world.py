class World:

    defaultatlas = 'wirecentralcontact100x50.bmp'

    def __init__(self, atlas=defaultatlas):
        self.atlas = atlas
        self.q = 1.6e-19
        self.hbar = 1.0545e-34/self.q
        self.m0 = 0.510e6/((3e8)**2)

        self.__read_geometry()
        self.__compose_nodes()

    def __read_geometry(self):
        from PIL import Image
        import scipy
        img = Image.open(self.atlas)
        arr = scipy.asarray(img)
        contact = []
        contact_shades = [149, 179, 209, 239]
        for shade in contact_shades:
            a = scipy.array(scipy.where(arr == shade))
            if a.shape[1] == 0: continue
            contact.append(a)
        conductor = scipy.transpose(scipy.where(arr > 0))

        self.canvas = arr.shape
        self.contact = contact
        self.conductor = conductor

    def __compose_nodes(self):
        from collections import OrderedDict
        import customiterators
        nodes = OrderedDict()
        count = 0
        for item,next_item in customiterators.pairwise(self.conductor):
            try:
                if item[1] == next_item[1]-1:
                    nodes[tuple(item)] = [count,count+1,None]
                else:
                    nodes[tuple(item)] = [count,None,None]
            except TypeError:
                nodes[tuple(item)] = [count,None,None]
            if item[0]>1:
                try:
                    nodes[tuple(item - [1,0])][2] = count
                except KeyError:
                    pass
            count +=1
        self.nodes = nodes
