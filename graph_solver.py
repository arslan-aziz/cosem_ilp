import pylp
import networkx as nx
import numpy as np

class GraphSolver:

    def __init__(self,graph_in,classA,classB,classC):

        self.graph = graph_in

        #n_vertices = graph.get_num_vertices()
        self.backend = pylp.create_linear_solver(pylp.Preference.Scip)
        #self.graph = graph
        #vertices = graph.vertices
        NUM_LEAVES = len(self.graph.nodes[0]['cost'])
        NUM_GRAPH_NODES = len(self.graph.nodes)
        NUM_GRAPH_EDGES = len(self.graph.edges)
        NUM_SOLVER_COEFF = NUM_GRAPH_NODES * NUM_LEAVES + NUM_GRAPH_EDGES

        GRAPH_NODES = list(self.graph.nodes)
        GRAPH_EDGES = list(self.graph.edges)
        self.backend.initialize(NUM_SOLVER_COEFF,pylp.VariableType.Binary)

        # backend.set_num_threads(1)
        objective = pylp.LinearObjective(NUM_SOLVER_COEFF)

        print('adding objective node coefficients')
        #add variables to objective per vertex
        id = 0
        self.vertex_id = []
        self.vertex_class_id = []
        for v in GRAPH_NODES:
            leaf_id = []
            class_id = []
            for leaf in range(NUM_LEAVES):
                objective.set_coefficient(id,-1*self.graph.nodes[v]['cost'][leaf])
                leaf_id.append(id)
                id+=1
            # add indicator variables for class A,B,C as well
            for c in range(3):
                objective.set_coefficient(id,1)
                class_id.append(id)
                id+=1
            self.vertex_id.append(leaf_id)
            self.vertex_class_id.append(class_id)

        assert (id == NUM_GRAPH_NODES*NUM_LEAVES)

        print('adding objective edge coefficients')
        self.edge_nodes = []
        self.obj_edge_id = []
        for e in GRAPH_EDGES:
            objective.set_coefficient(id,1)
            self.edge_nodes.append(e)
            self.obj_edge_id.append(id)
            id+=1

        assert(id == NUM_SOLVER_COEFF)

        # #add variables to objective per vertex
        # binary_id = 0
        # vertex_to_binary = {}
        # binary_to_vertex = {}
        # for v in g1.get_vertex_iterator():
        #     objective.set_coefficient(binary_id,graph[v])
        #
        #     vertex_to_binary[v] = binary_id
        #     binary_to_vertex[binary_id] = v
        #     binary_id += 1
        #
        # assert (binary_id == n_vertices)

        # #add variables to objective per edge
        # edge_to_binary = {}
        # binary_to_edge = {}
        # for e in g1.get_edge_iterator():
        #     objective.set_coefficient(binary_id,
        #                                 edge_cost[e])
        #     edge_to_binary[e] = binary_id
        #     binary_to_edge[binary_id] = e
        #     binary_id += 1
        # assert(binary_id == n_vertices + n_edges)
        #
        # self.backend.set_objective(objective)

        print('adding constraints')

        self.constraints = pylp.LinearConstraints()

        #add node constraints: binary, one leaf per node
        for v in self.vertex_id:
            constraintV = pylp.LinearConstraint()
            for l in v:
                constraint1 = pylp.LinearConstraint()
                constraint2 = pylp.LinearConstraint()

                constraint1.set_coefficient(l,1)
                constraint1.set_relation(pylp.Relation.GreaterEqual)
                constraint1.set_value(0)

                constraint2.set_coefficient(l,1)
                constraint2.set_relation(pylp.Relation.LessEqual)
                constraint2.set_value(1)

                constraintV.set_coefficient(l,1)

                self.constraints.add(constraint1)
                self.constraints.add(constraint2)

            constraintV.set_relation(pylp.Relation.Equal)
            constraintV.set_value(1)
            self.constraints.add(constraintV)

        #add node constraints: binary, one class per node (redundant?)
        for v in self.vertex_class_id:
            constraintV = pylp.LinearConstraint()
            for l in v:
                constraint1 = pylp.LinearConstraint()
                constraint2 = pylp.LinearConstraint()

                constraint1.set_coefficient(l, 1)
                constraint1.set_relation(pylp.Relation.GreaterEqual)
                constraint1.set_value(0)

                constraint2.set_coefficient(l, 1)
                constraint2.set_relation(pylp.Relation.LessEqual)
                constraint2.set_value(1)

                constraintV.set_coefficient(l,1)

                self.constraints.add(constraint1)
                self.constraints.add(constraint2)

            constraintV.set_relation(pylp.Relation.Equal)
            constraintV.set_value(1)
            self.constraints.add(constraintV)

        #add node constraints: node leaf sets node class
        for v in self.vertex_id:
            constraintA = pylp.LinearConstraint()
            for a in classA:
                constraintA.set_coefficient(v[a],1)
            constraintA.set_coefficient(self.vertex_class_id[0],-1)
            constraintA.set_relation(pylp.Relation.Equal)
            constraintA.set_value(0)

            constraintB = pylp.LinearConstraint()
            for b in classB:
                constraintB.set_coefficient(v[b],1)
            constraintB.set_coefficient(self.vertex_class_id[1],-1)
            constraintB.set_relation(pylp.Relation.Equal)
            constraintB.set_value(0)

            constraintC = pylp.LinearConstraint()
            for c in classC:
                constraintC.set_coefficient(v[c],1)
            constraintC.set_coefficient(self.vertex_class_id[2],-1)
            constraintC.set_relation(pylp.Relation.Equal)
            constraintC.set_value(0)

            self.constraints.add(constraintA)
            self.constraints.add(constraintB)
            self.constraints.add(constraintC)

        #add edge constraints: edge can't span class A(lumen) to class C(cytosol)
        #TODO: double check
        for e in GRAPH_EDGES:
            constraint_orderF = pylp.LinearConstraint()
            constraint_orderR = pylp.LinearConstraint()
            start_idx = GRAPH_NODES.index(e[0])
            end_idx = GRAPH_NODES.index(e[1])

            constraint_orderF.add_coefficient(self.vertex_class_id[start_idx][0],1)
            constraint_orderF.add_coefficient(self.vertex_class_id[end_idx][2],1)
            constraint_orderF.add_relation(pylp.Relation.Equal)
            constraint_orderF.add_value(0)
            self.constraints.add(constraint_orderF)

            constraint_orderR.add_coefficient(self.vertex_class_id[start_idx][2],1)
            constraint_orderR.add_coefficient(self.vertex_class_id[end_idx][0],1)
            constraint_orderR.add_relation(pylp.Relation.Equal)
            constraint_orderR.add_value(0)
            self.constraints.add(constraint_orderR)


        #TODO: add constraint edges with different classes are cut




        # #Edge selection implies vertex selection
        # for e in g1.get_edge_iterator():
        #     v0 = e.source()
        #     v1 = e.target()
        #
        #     #TODO: XOR logic?
        #     constraint = pylp.LinearConstraint()
        #     constraint.set_coefficient(edge_to_binary[e], 2)
        #     constraint.set_coefficient(vertex_to_binary[v0], -1)
        #     constraint.set_coefficient(vertex_to_binary[v1], -1)
        #     constraint.set_relation(pylp.Relation.LessEqual)
        #     constraint.set_value(0)

        self.backend.set_objective(objective)
        self.backend.set_constraints(self.constraints)

    def solve(self):
        solution,msg = self.backend.solve()
        return solution

    def get_vertex_ids(self):
        return self.vertex_id

    def get_edge_nodes(self):
        return self.edge_nodes

    def get_obj_edge_id(self):
        return self.obj_edge_id

    def add_path_constraint(self,node_path):
        path_constraint = pylp.LinearConstraint()
        for i in range(len(node_path)-1):
            try:
                idx = self.edge_nodes.index((node_path[i],node_path[i+1]))
            except:
                idx = self.edge_nodes.index((node_path[i+1],node_path[i]))
            #edge_idx = np.argwhere(self.edge_nodes == (node_path[i],node_path[i+1]))
            #edge_id = self.edge_id[edge_idx]
            edge_id = self.obj_edge_id[idx]
            path_constraint.set_coefficient(edge_id,1)
        path_constraint.set_relation(pylp.Relation.NotEqual)
        path_constraint.set_value(2)
        self.constraints.add(path_constraint)
        self.backend.set_constraints(self.constraints)


