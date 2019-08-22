import zarr
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from supervoxel_graph import SupervoxelRAG
from graph_solver import GraphSolver
import os
import networkx as nx

def __binarize__(vol_in,thresh,mode):
    vol = np.copy(vol_in)
    if(mode == 'step'):
        vol[vol<thresh] = 0
        vol[vol>=thresh] = 1
    if(mode == 'delta'):
        vol[vol != thresh] = 0
        vol[vol == thresh] = 1
    return vol

def main():

    SUPERVOXELS = ROI['watershed_tot'][:,:,:]

    # #REMOVE PLASMA MEMBRANE AND NUCLEUS FROM THIS PROCESSING
    # PM_EVIDENCE  = ROI['plasma_membrane'][:,:,:]
    # NUCLEUS_EVIDENCE = ROI['nucleolus'][:,:,:] +\
    #     ROI['nucleus_lumen'][:,:,:] +\
    #     ROI['NE'][:,:,:] +\
    #     ROI['NE_lumen']
    #
    # PM_EVIDENCE /= np.max(PM_EVIDENCE)
    # NUCLEUS_EVIDENCE /= np.max(NUCLEUS_EVIDENCE)
    # PM_EVIDENCE = __binarize__(PM_EVIDENCE,0.5,'step')
    # NUCLEUS_EVIDENCE = __binarize__(NUCLEUS_EVIDENCE,0.5,'step')
    #


    GRAPH_PATH = 'supervoxel_graph.gpickle'

    if (os.path.isfile(GRAPH_PATH) != True):
        graph = SupervoxelRAG(data_source = ROI,
                                data_meta_path = DATA_META_PATH,
                                graph_path = GRAPH_PATH)
    elif (COMMAND==True):
        graph = SupervoxelRAG(data_source=ROI,
                                data_meta_path = DATA_META_PATH,
                                graph_path = GRAPH_PATH,
                                command = command,
                                mode = 'command')
    else:
        graph = nx.read_gpickle(GRAPH_PATH)

    #which indicator variable belongs to which class
    membrane_labels = [2,4,7,9,10,11,17,19,21,23] #classB
    cytosol_labels = [0,25] #classC
    lumen_labels = [1,3,5,6,8,12,13,14,15,16,17,18,20,22,24] #classA

    NUM_SUPERVOXELS = np.max(np.unique(SUPERVOXELS))
    NUM_LEAVES = len(DATA_META)
    NUM_GRAPH_NODES = len(graph.nodes)
    NUM_GRAPH_EDGES = len(graph.edges)
    NUM_SOLVER_COEFF = NUM_GRAPH_NODES*NUM_LEAVES + NUM_GRAPH_EDGES

    GRAPH_NODES = np.array(list(graph.nodes))
    GRAPH_EDGES = np.array(list(graph.edges))

    solver = GraphSolver(graph,lumen_labels,membrane_labels,cytosol_labels)

    edge_nodes = solver.get_edge_nodes()
    obj_edge_id = solver.get_obj_edge_id()

    errors = 100
    num_path_constraints = 0

    iteration = 0
    while (errors > 5):

        print('Iteration: ' + str(iteration))
        print('Errors Remaining: ' + str(errors))
        print('Constraints Added: ' + str(num_path_constraints))

        print('solving')

        results = solver.solve()

        print('done solving')

        assert(len(results) == NUM_SOLVER_COEFF)

        #graph vertex id matches supervoxel id
        vertex_array = []
        x=0
        for i in range(NUM_GRAPH_NODES):
            leaf_array = []
            for j in range(NUM_LEAVES):
                leaf_array.append(results[x])
                x+=1
            vertex_array.append(leaf_array)

        vertex_array = np.array(vertex_array)

        edge_array = []
        for j in range(x,x+NUM_GRAPH_EDGES):
            edge_array.append(results[j])

        assert(np.sum(vertex_array)==NUM_GRAPH_NODES)


        cytosol_supervoxel = []
        membrane_supervoxel = []

        for i in range(len(vertex_array)):
            max_leaf = np.argmax(vertex_array[i])
            graph.nodes[GRAPH_NODES[i]]['label'] = max_leaf
            if(max_leaf == NUM_LEAVES-1):
                cytosol_supervoxel.append(GRAPH_NODES[i])
            if(max_leaf in membrane_labels):
                membrane_supervoxel.append(GRAPH_NODES[i])

        if(SEPARATE_COMPONENTS==True):

            ##CREATE MERGE GRAPH
            nbunch = np.copy(GRAPH_NODES)
            delete = np.array(cytosol_supervoxel).astype(np.int)
            mask = np.isin(nbunch,delete)
            nbunch = np.logical_not(mask)*nbunch
            nbunch = np.delete(nbunch,np.argwhere(nbunch==0))
            # merge_graph = all lumen and membrane supervoxels, all adjacent edges connected
            #constructor required to "unfreeze" subgraph
            merge_graph = nx.Graph(graph.subgraph(nbunch))

            #TODO: remove edges NOT nodes
            ##REMOVE nodes corresponding to "NO MERGE" edges from merge_graph
            # no_merge_edges = np.argwhere(edge_array==1)
            # no_merge_nodes = GRAPH_EDGES[no_merge_edges]
            # merge_graph.remove_nodes_from(np.unique(no_merge_nodes).tolist())
            if(iteration != 0):
                no_merge_edges = [edge_nodes[i] for i in np.where(np.array(edge_array)==1)[0].tolist()]
                merge_graph.remove_edges_from(no_merge_edges)


            ##GENERATE connected components on merge_graph
            super_components = nx.connected_components(merge_graph)

            #On each component, remove membrane supervoxels and check if more than one component
            errors=0
            for comp in super_components:
                nbunch = np.array(list(comp))
                delete = np.array(membrane_supervoxel).astype(np.int)
                mask = np.isin(nbunch, delete)
                nbunch = np.logical_not(mask)*nbunch
                nbunch = np.delete(nbunch,np.argwhere(nbunch==0))
                comp_graph_no_membrane = nx.Graph(merge_graph.subgraph(nbunch))

                sub_components = [c for c in sorted(nx.connected_components(comp_graph_no_membrane),key=len,reverse=True)]

                #TODO: repeat for pairwise components?
                if(len(sub_components)>1):
                    errors += len(sub_components)
                    #add source and sink nodes
                    source_node = np.max(GRAPH_NODES)+1
                    merge_graph.add_node(source_node)
                    for i in list(sub_components[0]):
                        merge_graph.add_edge(source_node,i)
                    sink_node = np.max(GRAPH_NODES)+2
                    merge_graph.add_node(sink_node)
                    for i in list(sub_components[1]):
                        merge_graph.add_edge(sink_node,i)

                    #find shortest path from source to sink
                    short_path = nx.shortest_path(merge_graph,source=source_node,target=sink_node)
                    #TODO: are source and sink indices 0 and last on path list?
                    short_path.remove(source_node)
                    short_path.remove(sink_node)
                    merge_graph.remove_node(source_node)
                    merge_graph.remove_node(sink_node)
                    solver.add_path_constraint(short_path)
                    num_path_constraints += 1

                # if(len(sub_components)>1):
                #     errors += len(sub_components)
                #     # select random node in each component
                #     rnode1 = np.array(list(sub_components[0]))
                #     rnode1 = rnode1[np.random.randint(0,len(rnode1))]
                #     rnode2 = np.array(list(sub_components[1]))
                #     rnode2 = rnode2[np.random.randint(0,len(rnode2))]
                #     #rnode1 = sub_components[0][np.random.randint(0,len(sub_components[0]))]
                #     #rnode2 = sub_components[1][np.random.randint(0,len(sub_components[1]))]
                #     # get shortest path in merge_graph and add as constraint
                #     solver.add_path_constraint(nx.shortest_path(merge_graph,source=rnode1,target=rnode2))
                #     num_path_constraints+=1
        iteration+=1


    print('done solving after ' + str(iteration) + ' iterations')

    LABELED_LEAF_OUTPUT = np.zeros(SUPERVOXELS.shape)
    for i in graph.nodes:
        indices = np.where(SUPERVOXELS==i)
        LABELED_LEAF_OUTPUT[indices[0],indices[1],indices[2]] =\
            graph.nodes[i]['label']

    MITO_EXAMPLE = np.zeros(SUPERVOXELS.shape)
    final_graph = graph.copy()
    no_merge_edges = np.argwhere(edge_array==1)
    no_merge_edges = edge_nodes[no_merge_edges]
    final_graph.remove_nodes_from(no_merge_edges)
    components = nx.connected_components(merge_graph)
    x=1
    for c in components:
        ids = list(c)
        for i in ids:
            MITO_EXAMPLE[SUPERVOXELS==i] = x
        x+=1

    OUTPUT_PATH = '/groups/funke/home/aziza/Documents/cosem_scripts/final.zarr'
    f = zarr.open(OUTPUT_PATH, mode='w')
    f['processed_components'] = LABELED_LEAF_OUTPUT.astype(np.uint64)
    f['processed_components'].attrs['resolution'] = (4,4,4)
    f['watershed'] = SUPERVOXELS.astype(np.uint64)
    f['watershed'].attrs['resolution'] = (4,4,4)
    f['separated'] = MITO_EXAMPLE.astype(np.uint64)
    f['separated'].attrs['resolution'] = (4,4,4)


if __name__ == "__main__":



    # with open('meta/cosem_metadata.json', 'r') as f:
    #     META_DATA = json.load(f)

    print('dataset opened')

    # ROI_NAME = META_DATA['ROI_NAME']
    # ORG_KEYS = META_DATA['ORG_KEYS']
    # COMB_YHAT = META_DATA['COMB_YHAT']

    DATA_SOURCE_PATH = '/groups/funke/home/aziza/Documents/cosem_scripts/watershed.zarr'
    DATA_META_PATH = 'watershed_metadata.json'

    ROI = zarr.open(DATA_SOURCE_PATH, mode='r')
    with open(DATA_META_PATH,'r') as fm:
        DATA_META = json.load(fm)

    COMMAND = False
    command = 'add_edges'

    SEPARATE_COMPONENTS = True

    main()