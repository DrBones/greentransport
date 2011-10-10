class World:

    def __init__(self, atlas=None):
        from parameters import p

        self.p = p
        self.atlas = atlas
        self.__read_geometry()
        #self.__compose_nodes()
        # self.__blocksizes_from_coords()

    def __read_geometry(self):
        """reads in bmp and generates contacts and other
        objects of interest"""
        from PIL import Image
        from scipy import where, asarray, array, transpose,logical_or,logical_and
        from aux import Contact
        img = Image.open(self.atlas)
        arr = asarray(img)
        self.p.canvas = arr
        self.p.raw_coords =logical_and(arr > 0,arr %5 ==0).nonzero()
        self.p.tuple_canvas_coordinates = tuple(zip(*self.p.raw_coords))
        contacts = []
        shades = [(103,115), (133,145), (163,175), (193,205)]
        contact_index = 0
        for shade in shades:
            contact_temp_tuple_coords = tuple(zip(*where(arr == shade[1])))
            interface_temp_tuple_coords = tuple(zip(*logical_or(arr == shade[0],arr ==shade[1]).nonzero()))
            a = Contact(contact_temp_tuple_coords, interface_temp_tuple_coords)
            if a.length == 0: continue
            a.index = contact_index
            # try:
            #     if any(arr[a.interface_raw_coordinates[0][0]+1,:]==239):
            #         a.SO = True
            #     else:
            #         a.SO = False
            # except IndexError:
            #     if any(arr[a.interface_raw_coordinates[0][0]-1,:] == 239):
            #         a.SO = True
            #     else:
            #         a.SO = False
            contacts.append(a)
            contact_index +=1
        if len(contacts) == 0:
            print '--------------------------------------'
            print 'No contacts found! This will not work!'
            print '--------------------------------------'
        # from pudb import set_trace; set_trace()
        self.p.contacts = contacts

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
