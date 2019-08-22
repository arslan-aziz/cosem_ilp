import networkx as nx
import numpy as np
import zarr
import json
from scipy import ndimage
from skimage.future.graph import RAG

class SupervoxelRAG:

    def compute_vertex_cost(self,indices):
        cost_array = []

        for j in range(len(self.ROI_META)):
            cost_array.append(np.sum(self.data[j][indices[0],indices[1],indices[2]]))

        return np.array(cost_array).astype(np.float64)

    def __init__(self,data_source,data_meta_path,graph_path):
        self.ROI = data_source
        with open(data_meta_path, 'r') as f:
            self.ROI_META = json.load(f)

        self.labeled_supervoxels = self.ROI['watershed_tot'][:, :, :]
        self.N_vertices = len(np.unique(self.labeled_supervoxels))

        self.data = []
        for leaf in self.ROI_META:
            self.data.append(self.ROI[leaf][:, :, :])

        print('creating initial graph')
        self.G = RAG(label_image=self.labeled_supervoxels)

        print('adding node attributes')
        x=0
        for i in self.G.nodes:
            indices = np.where(self.labeled_supervoxels == i)
            self.G.nodes[i]['cost'] = self.compute_vertex_cost(indices)
            self.G.nodes[i]['label'] = -1
            #self.G.nodes[i]['neighbors'] = []
            if (x % 100 == 0):
                print(x)
            x+=1

        self.nodes = self.G.nodes
        self.edges = self.G.edges

        print('saving graph')
        nx.write_gpickle(self.G, graph_path)

class SupervoxelGraph:

    def compute_vertex_cost(self,indices):
        cost_array = []

        for j in range(len(self.ROI_META)):
            cost_array.append(np.sum(self.data[j][indices[0],indices[1],indices[2]]))

        return cost_array

    def compute_shared_boundary(self,u,v):
        u_coor = np.array(np.where(self.labeled_supervoxels==u))
        v_coor = np.array(np.where(self.labeled_supervoxels==v))
        u_coor = self.to_origin(u_coor)+1
        v_coor = self.to_origin(v_coor)+1

        vol = np.zeros((np.max(np.concatenate((u_coor[0],v_coor[0])))+2,
                       np.max(np.concatenate((u_coor[1],v_coor[1])))+2,
                       np.max(np.concatenate((u_coor[2],v_coor[2])))+2))
        vol[u_coor[0],u_coor[1],u_coor[2]] = 1
        vol[v_coor[0],v_coor[1],v_coor[2]] = 2

        delta = np.array([-1,0,1])

        if(len(u_coor[0])<len(v_coor[0])):
            for i in delta:
                for j in delta:
                    for k in delta:
                        neighbors = vol[u_coor[0]+i,u_coor[1]+j,u_coor[2]+k]
                        if(np.any(neighbors==2)):
                            return True
        else:
            for i in delta:
                for j in delta:
                    for k in delta:
                        neighbors = vol[v_coor[0]+i,v_coor[1]+j,v_coor[2]+k]
                        if(np.any(neighbors==2)):
                            return True


        return False

        #
        # if(len(u_coor[0])<len(v_coor[0])):
        # #check all 26 neighbors
        #
        #     neighbors = [vol[u_coor[0],u_coor[1],u_coor[2]],
        #                  vol[u_coor[0]+1,u_coor[1],u_coor[2]],
        #                  vol[u_coor[0],u_coor[1]+1,u_coor[2]],
        #                  vol[u_coor[0],u_coor[1],u_coor[2]+1],
        #                  vol[u_coor[0]+1,u_coor[1]+1,u_coor[2]],
        #                  vol[u_coor[0]+1,u_coor[1],u_coor[2]+1],
        #                  vol[u_coor[0],u_coor[1]+1,u_coor[2]+1],
        #                  vol[u_coor[0]+1,u_coor[1]+1,u_coor[2]+1]]
        #
        #     neighbors = np.array(neighbors)
        #
        # else:
        #     pass
        #
        # return (np.any(neighbors==1) & np.any(neighbors==2))



        #num_comp = ndimage.measurements.label(vol)[1]
        #
        # if(num_comp==1):
        #     return True
        # else:
        #     return False

        # sx = ndimage.sobel(vol,axis=0,mode='constant')
        # sy = ndimage.sobel(vol,axis=1,mode='constant')
        # sz = ndimage.sobel(vol,axis=2,mode='constant')
        #
        # edge_map = sx+sy+sz
        # edge_map[edge_map != 0] = 1
        #
        # edge_map = ndimage.morphology.binary_erosion(edge_map)

        #return np.sum(edge_map)


    def to_origin(self,coord):
        coord[0] -= np.min(coord[0])
        coord[1] -= np.min(coord[1])
        coord[2] -= np.min(coord[2])
        return coord

    def compute_edge_cost(self,vertex_1,vertex_2):
        #TODO
        return 0

    def get_num_vertices(self):
        return self.N_vertices

    def get_num_edges(self):
        return self.N_edges

    def get_graph(self):
        return self.G

    def add_edges(self):
        print('adding graph edges')
        # Edge added if adjacent
        for i in range(1, self.N_vertices + 1):
            for j in range(1, self.N_vertices + 1):
                if (i != j):
                    s_bound = self.compute_shared_boundary(i, j)
                    if (s_bound == True):
                        # self.G.add_edge(i,j,shared_boundary = s_bound)
                        self.G.add_edge(i, j)
                        self.G.nodes[i]['neighbors'].append(j)
                        self.G.nodes[j]['neighbors'].append(i)


    def __init__(self,data_source,data_meta_path,graph_path,command=None,mode='new'):
        if(mode=='new'):
            self.G = nx.Graph()

            self.ROI = data_source
            with open(data_meta_path,'r') as f:
                self.ROI_META = json.load(f)

            self.labeled_supervoxels = self.ROI['watershed_tot'][:,:,:]
            self.N_vertices = len(np.unique(self.labeled_supervoxels))


            self.data = []
            for leaf in self.ROI_META:
                self.data.append(self.ROI[leaf][:,:,:])



            print('adding graph nodes')
            for i in range(1,self.N_vertices+1):
                indices = np.where(self.labeled_supervoxels==i)
                #assert (len(indices[0]) != 0), ('iteration' + str(i))
                #TODO: why do some watershed components have zero length?
                if(len(indices[0])!=0):
                    self.G.add_node(i,cost = self.compute_vertex_cost(indices), label = -1, neighbors=[])
                if(i%100==0):
                    print(i)


            print('adding graph edges')
            #Edge added if adjacent
            for i in self.G.nodes:
                for j in self.G.nodes:
                    if(i!=j):
                        s_bound = self.compute_shared_boundary(i,j)
                        if(s_bound==True):
                            #self.G.add_edge(i,j,shared_boundary = s_bound)
                            self.G.add_edge(i,j)
                            self.G.nodes[i]['neighbors'].append(j)
                            self.G.nodes[j]['neighbors'].append(i)
                    if((i-1)*j+j%1000==0):
                        print((i-1)*j+j)

            self.N_edges = len(self.G.edges)

            nx.write_gpickle(self.G,graph_path)

        else:
            self.G  = nx.read_gpickle(graph_path)
            self.ROI = data_source
            with open(data_meta_path, 'r') as f:
                self.ROI_META = json.load(f)

            self.labeled_supervoxels = self.ROI['watershed_tot'][:, :, :]
            self.N_vertices = len(np.unique(self.labeled_supervoxels))

            eval('self.'+command+'()')


