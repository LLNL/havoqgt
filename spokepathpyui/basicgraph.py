class BasicGraph :
    def __init__(self, edge_list_filename = None) :
        self.vertices = list()
        self.edges = list()
        
        if edge_list_filename is not None :
            self.read_edge_list_from_file(edge_list_filename)
        
    def add_vertex(self, vertex) :
        self.vertices.append(vertex)
    
    def add_edge(self, source_vertex, targe_vertex) :
        self.edges.append((source_vertex, targe_vertex))
        
    def read_edge_list_from_file(self, filename) :
        tmp_vertices = set()
        with open (filename, "r") as read_file :
            for line in read_file :
                line = line.strip()
                tokens = line.split(" ")
                s = int(tokens[0].strip(' \t\n\r'))
                t = int(tokens[1].strip(' \t\n\r'))
                self.add_edge(s, t)
                tmp_vertices.add(s)
                tmp_vertices.add(t)
                
        self.vertices.clear()
        self.vertices.extend(sorted(tmp_vertices))
        
    def read_edge_list(self, edge_list) :
        tmp_vertices = set()    
        for line in edge_list :
            line = line.strip()
            tokens = line.split(" ")
            s = int(tokens[0].strip(' \t\n\r'))
            t = int(tokens[1].strip(' \t\n\r'))
            self.add_edge(s, t)
            tmp_vertices.add(s)
            tmp_vertices.add(t)

        self.vertices.clear()
        self.vertices.extend(sorted(tmp_vertices))
        
    def __str__(self) :
        return "Vertices(" + str(len(self.vertices)) + "): " + str(self.vertices) + "\n" + \
    "Edges(" + str(len(self.edges)) + "): " + str(self.edges)
