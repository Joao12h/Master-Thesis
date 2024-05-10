import networkx as nx
from networkx.utils import py_random_state
import random
import numpy as np
from itertools import combinations, islice,product
import copy


def havel_hakimi2(degree_distribution):
    nb_links=sum(np.array(degree_distribution)[:,1])//2
    simple_graph=[]
    sorted_degree_distribution=copy.deepcopy(degree_distribution)

    while sorted_degree_distribution!=[]:
        random.shuffle(sorted_degree_distribution)
        sorted_degree_distribution=sorted(sorted_degree_distribution, reverse=True, key=lambda x: x[1])
        ego=sorted_degree_distribution.pop(0)
        index=0
        sub_simple_graph=[]

        while ego[1]>0:
            sub_simple_graph+=[[ego[0],sorted_degree_distribution[index][0]]]
            ego[1]=ego[1]-1
            sorted_degree_distribution[index][1]=sorted_degree_distribution[index][1]-1
            index=index+1
        simple_graph+=sub_simple_graph
    return (simple_graph, len(simple_graph)==nb_links)

def double_edge_swap(initial_graph, nswap=1, max_tries=100, seed=None):

    G=initial_graph.copy()
    if G.is_directed():
        raise nx.NetworkXError(
            "double_edge_swap() not defined for directed graphs. Use directed_edge_swap instead."
        )
    if nswap > max_tries:
        raise nx.NetworkXError("Number of swaps > number of tries allowed.")
    if len(G) < 4:
        raise nx.NetworkXError("Graph has fewer than four nodes.")
    if len(G.edges) < 2:
        raise nx.NetworkXError("Graph has fewer than 2 edges")
    # Instead of choosing uniformly at random from a generated edge list,
    # this algorithm chooses nonuniformly from the set of nodes with
    # probability weighted by degree.
    candidate_edges=[]
    edges_proposed=[]
    n = 0
    swapcount = 0
    keys, degrees = zip(*G.degree())  # keys, degree
    cdf = nx.utils.cumulative_distribution(degrees)  # cdf of degree
    discrete_sequence = nx.utils.discrete_sequence
    while swapcount < nswap:
        #print('n:',n)
        #        if random.random() < 0.5: continue # trick to avoid periodicities?
        # pick two random edges without creating edge list
        # choose source node indices from discrete distribution
        (ui, xi) = discrete_sequence(2, cdistribution=cdf, seed=seed)
        if ui == xi:
            continue  # same source, skip
        u = keys[ui]  # convert index to label
        #print('u:',u)
        x = keys[xi]
        #print('x:',x)
        # choose target uniformly from neighbors
        v = random.choice(list(G[u]))
        y = random.choice(list(G[x]))

        if v == y:
            continue  # same target, skip
        if (x not in G[u]) and (y not in G[v]):  # don't create parallel edges
            G.add_edge(u, x)
            G.add_edge(v, y)
            G.remove_edge(u, v)
            G.remove_edge(x, y)
            swapcount += 1
        if n >= max_tries:
            e = (
                f"Maximum number of swap attempts ({n}) exceeded "
                f"before desired swaps achieved ({nswap})."
            )
            raise nx.NetworkXAlgorithmError(e)
        n += 1
        candidate_edges=[(u,v),(x,y)]
        edges_proposed=[(u,x),(v,y)]

    return G,candidate_edges,edges_proposed

#comb_lists_different_neighbours gets edge_lists of the neighbours of two different nodes that are not connected to each other. It will be useful to increase the clustering coefficient and maximize number of triangles formed per iteration
def comb_lists_different_neighbours(neighbourhood1, neighbourhood2):
    set_neighbourhood1 = set(neighbourhood1)
    set_neighbourhood2 = set(neighbourhood2)

    return list(product(set_neighbourhood1, set_neighbourhood2))


def markov_chain(initial_graph, nb_iterations=2,seed=None):
    G=initial_graph.copy()
    if G.is_directed():
        raise nx.NetworkXError(
            "double_edge_swap() not defined for directed graphs. Use directed_edge_swap instead."
        )
    if len(G) < 4:
        raise nx.NetworkXError("Graph has fewer than four nodes.")
    if len(G.edges) < 2:
        raise nx.NetworkXError("Graph has fewer than 2 edges")
    # Instead of choosing uniformly at random from a generated edge list,
    # this algorithm chooses nonuniformly from the set of nodes with
    # probability weighted by degree.
    swapcount = 0
    keys, degrees = zip(*G.degree())  # keys, degree
    cdf = nx.utils.cumulative_distribution(degrees)  # cdf of degree
    discrete_sequence = nx.utils.discrete_sequence
    current_transitivity=nx.transitivity(G)
    largest_cc = max(nx.connected_components(G), key=len)
    S = G.subgraph(largest_cc).copy()
    current_diameter=nx.diameter(S)
    current_assortativity=nx.degree_assortativity_coefficient(G)
    list_of_transitivity=[current_transitivity]
    list_of_diameter=[current_diameter]
    list_of_assortativity=[current_assortativity]
    list_of_sample_mean_diameter=[current_diameter]
    iteration=1

    while iteration <= nb_iterations:
        print('iteration:',iteration)
        #        if random.random() < 0.5: continue # trick to avoid periodicities?
        # pick two random edges without creating edge list
        # choose source node indices from discrete distribution
        # Generate a random number between 0 and 1
        random_number = random.random()
        if random_number<0.5:
            (ui, xi) = discrete_sequence(2, cdistribution=cdf, seed=seed)
            if ui != xi:
                u = keys[ui]  
                x = keys[xi]
                v = random.choice(list(G[u]))
                y = random.choice(list(G[x]))
                if v != y:
                    if (x not in G[u]) and (y not in G[v]):  # don't create parallel edges
                        G.add_edge(u, x)
                        G.add_edge(v, y)
                        G.remove_edge(u, v)
                        G.remove_edge(x, y)
                        swapcount += 1
        if iteration%3000==0:
            current_transitivity=nx.transitivity(G)
            largest_cc = max(nx.connected_components(G), key=len)
            S = G.subgraph(largest_cc).copy()
            current_diameter=nx.diameter(S)
            current_assortativity=nx.degree_assortativity_coefficient(G) 

            list_of_transitivity.append(current_transitivity)
            list_of_diameter.append(current_diameter)
            list_of_sample_mean_diameter.append(np.mean(list_of_diameter))
            list_of_assortativity.append(current_assortativity)
        iteration+=1
        print('iteration:',iteration)
    
    return G,list_of_transitivity,list_of_diameter,list_of_assortativity,list_of_sample_mean_diameter,swapcount



def clustering_maybe_fast(initial_graph, nb_iterations=2,seed=None):
    G=initial_graph.copy()
    if G.is_directed():
        raise nx.NetworkXError(
            "double_edge_swap() not defined for directed graphs. Use directed_edge_swap instead."
        )
    if len(G) < 4:
        raise nx.NetworkXError("Graph has fewer than four nodes.")
    if len(G.edges) < 2:
        raise nx.NetworkXError("Graph has fewer than 2 edges")
   
    #save_graphs=[G]
    current_nb_triangles=sum(nx.triangles(G).values())/3
    number_of_triangles_per_iteration=[current_nb_triangles]
    swapcount = 0
    keys, degrees = zip(*G.degree())  # keys, degree
    cdf = nx.utils.cumulative_distribution(degrees)  # cdf of degree
    discrete_sequence = nx.utils.discrete_sequence
    iteration=1

    while iteration <= nb_iterations:
        print('iteration:',iteration)
        #        if random.random() < 0.5: continue # trick to avoid periodicities?
        # pick two random edges without creating edge list
        # choose source node indices from discrete distribution
       
        # Generate a random number between 0 and 1
        random_number = random.random()
        #print('random_number:',random_number)
        
        (ui, xi) = discrete_sequence(2, cdistribution=cdf, seed=seed)
        if ui != xi:
            u = keys[ui]
            x = keys[xi]
            comb_lists=comb_lists_different_neighbours(G[u],G[x])
            number_common_neighbours_current_per_edge_of_neighbourhood=[len(set(nx.common_neighbors(G,edge[0], edge[1]))) for edge in comb_lists]
        
            all_common_neighbours=sum(number_common_neighbours_current_per_edge_of_neighbourhood)

            if all_common_neighbours==0:
                weights=None
            else:
                weights=[number_common_neighbours_current_per_edge_of_neighbourhood[i]/sum(number_common_neighbours_current_per_edge_of_neighbourhood) for i in range(0,len(number_common_neighbours_current_per_edge_of_neighbourhood))]

            random_sublist = random.choices(comb_lists, weights=weights, k=1)
            #print('random_sublist:', random_sublist)

            v = random_sublist[0][0]
            #print('v:',v)
            y = random_sublist[0][1]
            #print('y:',y)
                    #print('y:',y)
            if v != y:
                if (x not in G[u]) and (y not in G[v]):
                    candidate_edges=[(u,v),(x,y)]
                    edges_proposed=[(u,x),(v,y)]

                    number_common_neighbours_current_per_edge=[len(set(nx.common_neighbors(G,edge[0], edge[1]))) for edge in candidate_edges]

                    number_common_neighbours_proposed_per_edge=[len(set(nx.common_neighbors(G,edge[0], edge[1]))) for edge in edges_proposed]
                
                    soma_number_common_neighbours_current=sum(number_common_neighbours_current_per_edge)

                    soma_number_common_neighbours_proposed=sum(number_common_neighbours_proposed_per_edge)

                    if soma_number_common_neighbours_proposed>soma_number_common_neighbours_current:
                        gain_nb_triangles=soma_number_common_neighbours_proposed-soma_number_common_neighbours_current
                        current_nb_triangles+=gain_nb_triangles
                        G.add_edge(u, x)
                        G.add_edge(v, y)
                        G.remove_edge(u, v)
                        G.remove_edge(x, y)
                        swapcount += 1            

        if iteration%200==0: 
            number_of_triangles_per_iteration.append(current_nb_triangles)
        iteration+=1
        

    return G,number_of_triangles_per_iteration,swapcount


def increase_assortativity_degree(initial_graph,nb_iterations=2):
    iteration=1
    current_graph=initial_graph.copy()
    save_assortativity_coefficient=[]
    assortativity_coefficient=nx.degree_assortativity_coefficient(current_graph)
    save_graphs=[]

    while iteration<=nb_iterations:
        proposed_graph=double_edge_swap(current_graph, nswap=1, max_tries=100, seed=None)

        if   current_graph.degree(proposed_graph[2][0][0])*current_graph.degree(proposed_graph[2][0][1])+\
                    current_graph.degree(proposed_graph[2][1][0])*current_graph.degree(proposed_graph[2][1][1])-\
                    (current_graph.degree(proposed_graph[1][0][0])*current_graph.degree(proposed_graph[1][0][1])+\
                    current_graph.degree(proposed_graph[1][1][0])*current_graph.degree(proposed_graph[1][1][1]))>0:
                    current_graph=proposed_graph[0]
                    assortativity_coefficient=nx.degree_assortativity_coefficient(current_graph)

        if iteration%200==0:
                save_assortativity_coefficient.append(assortativity_coefficient)
        if iteration%40000==0:
                save_graphs.append()

        iteration+=1
            
    return save_graphs,save_assortativity_coefficient

def increase_assortativity_race(initial_graph,atribute,nb_iterations=2):
    iteration=1
    current_graph=initial_graph.copy()
    save_assortativity_coefficient=[]
    assortativity_coefficient=nx.attribute_assortativity_coefficient(current_graph,atribute)
    save_graphs=[]

    while iteration<=nb_iterations:

        proposed_graph=double_edge_swap(current_graph, nswap=1, max_tries=100, seed=None)

        if nx.attribute_assortativity_coefficient(proposed_graph[0],atribute)>nx.attribute_assortativity_coefficient(current_graph,atribute):
                    current_graph=proposed_graph[0]
                    assortativity_coefficient=nx.attribute_assortativity_coefficient(current_graph,atribute)
        
        if iteration%200==0:
                save_assortativity_coefficient.append(assortativity_coefficient)
        if iteration%40000==0:
                save_graphs.append()

        iteration+=1
    return save_graphs,save_assortativity_coefficient

def increase_assortativity_age(initial_graph,atribute,nb_iterations=2):
    iteration=1
    current_graph=initial_graph.copy()
    save_assortativity_coefficient=[]
    assortativity_coefficient=nx.numeric_assortativity_coefficient(current_graph,atribute)
    save_graphs=[]

    while iteration<=nb_iterations:
        proposed_graph=double_edge_swap(current_graph, nswap=1, max_tries=100, seed=None)

        if nx.numeric_assortativity_coefficient(proposed_graph[0],atribute)>nx.numeric_assortativity_coefficient(current_graph,atribute):
                current_graph=proposed_graph[0]
                assortativity_coefficient=nx.numeric_assortativity_coefficient(current_graph,atribute)
        
        if iteration%200==0:
                save_assortativity_coefficient.append(assortativity_coefficient)
        if iteration%40000==0:
                save_graphs.append()

        iteration+=1
    return save_graphs,save_assortativity_coefficient


