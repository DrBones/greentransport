class Contact(object):
    import aux

    def __init__(self, contact_temp_tuple_coords, interface_temp_tuple_coords, SO=None,):
        from parameters import p
        self.p = p
        # self.interface_raw_coordinates = None
        # self.contact_raw_coordinates = sc.asarray(contact_temp_coords)
        # the star unpacks the array into lists first and zip zips them to (x,y) tuples
        self.interface_tuple_coordinates = interface_temp_tuple_coords
        self.contact_tuple_coordinates = contact_temp_tuple_coords
        self.length = len(self.contact_tuple_coordinates)
        if self.length != 0:
            print "Generating graph for interface" 
            self.graph = Contact.aux.digraph_from_tuple_coords(self.interface_tuple_coordinates)
            self.add_nodenames()
            self.transverse_sorted_nodelist_from_graph()

    def add_nodenames(self):
        self.internal_node_names = set()
        self.node_names = set()
        for initial_tuple in self.contact_tuple_coordinates:
            initial_name = self.p.tuple_canvas_coordinates.index(initial_tuple)
            self.node_names.add(initial_name)
            self.internal_node_names.add(self.interface_tuple_coordinates.index(initial_tuple))
            #for partner_tuple in contact_tuple:
            #    partner_name = self.tuple_of_coords.index(partner_tuple)
            #    self.graph.add_edge(initial_name,partner_name)

    def transverse_sorted_nodelist_from_graph(self):
        self.levelstructure = list(Contact.aux.BreadthFirstLevels(self.graph, self.internal_node_names))
        self.nodelist = [item for level in self.levelstructure for item in sorted(level)]
