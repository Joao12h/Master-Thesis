#libraries
import os, pickle
import sys
import networkx as nx
import numpy as np
import pandas as pd
import community
from pathlib import Path

property_under_study=sys.argv[1]

#auxiliary functions

def flatten_comprehension(matrix):
        return [item for row in matrix for item in row]

def graph_with_atributes(graph,agents):
    for nodes in graph.nodes():
            #print(nodes)
            index=agents.loc[agents['Participant']==nodes].index[0]
            attrs = {nodes: {"age": agents["age"][index],"Number of partners": agents["number_of_partners"][index], "race": agents["race"][index] , "HIV_status": agents["HIV_status"][index],"Sexual_position": 'versatile',"PrEP adherence": 0,"Testing Frequency": agents["TEST2YRS"][index]/2,"Condom usage": 0}}
            nx.set_node_attributes(graph, attrs)
    return graph

if property_under_study=='Clustering Coefficient':

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'agents_complete_final_use_in_analysis.pkl'), 'rb') #rb = read bytes because we are reading the file
    agents = pickle.load(file)

    #Clustering

    list_of_clustered_networks_15=[]
    id_networks=list(range(0,7+1,1))
    directory_path_15 = "/home/jbrazia/metrics/Network_15/Clustering"

    for id in id_networks:
        file_name_15 = "clustering_network_500000_iterations_15_latest_stochastic_ID_"+str(id)+".pkl"
        file_path_15 = os.path.join(directory_path_15, file_name_15)

        with open(file_path_15, "rb") as open_file_15:
            stationary_networks_15 = pickle.load(open_file_15)
            list_of_clustered_networks_15.append(stationary_networks_15)

    list_of_clustered_networks_16=[]
    id_networks=list(range(0,7+1,1))
    directory_path_16 = "/home/jbrazia/metrics/Network_16/Clustering"

    for id in id_networks:
        file_name_16 = "clustering_network_500000_iterations_16_latest_stochastic_ID_"+str(id)+".pkl"
        file_path_16 = os.path.join(directory_path_16, file_name_16)

        with open(file_path_16, "rb") as open_file_16:
            stationary_networks_16 = pickle.load(open_file_16)
            list_of_clustered_networks_16.append(stationary_networks_16)

    list_of_clustered_networks_17=[]
    id_networks=list(range(0,7+1,1))
    directory_path_17 = "/home/jbrazia/metrics/Network_17/Clustering"

    for id in id_networks:
        file_name_17 = "clustering_network_500000_iterations_17_latest_stochastic_ID_"+str(id)+".pkl"
        file_path_17 = os.path.join(directory_path_17, file_name_17)

        with open(file_path_17, "rb") as open_file_17:
            stationary_networks_17 = pickle.load(open_file_17)
            list_of_clustered_networks_17.append(stationary_networks_17)

    list_of_all_clustered_networks=[list_of_clustered_networks_15[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_clustered_networks_15))]+[list_of_clustered_networks_16[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_clustered_networks_16))]+[list_of_clustered_networks_17[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_clustered_networks_17))]
    
    flattened_list_of_all_clustered_networks = flatten_comprehension(list_of_all_clustered_networks)


    
    for graph in flattened_list_of_all_clustered_networks:
        graph=graph_with_atributes(graph,agents)
    
    #group lists according to timesteps
    
    grouped_lists=[]
    for timestep in range(0,len(list_of_clustered_networks_16[0][0])):
        graph_ID=timestep
        group_per_time_step=[]
        while graph_ID<len(flattened_list_of_all_clustered_networks):
            group_per_time_step.append(flattened_list_of_all_clustered_networks[graph_ID])
            graph_ID+=6
        grouped_lists.append(group_per_time_step)

    #global properties per generated graph, we compute and compare the properties of the 
    counter=0
    clustering_global=[]
    assortativity_degree_global=[]
    assortativity_age_global=[]
    assortativity_race_global=[]
    average_diameter_global=[]
    average_shortest_path_length_global=[]
    largest_component_size_global=[]
    number_of_subgraphs_global=[]
    modularity_global=[]
    betweenness_global=[]



    #assortativity by HIV status and sexual role checked in the largest component size graph
    for timestep in range(0,len(list_of_clustered_networks_16[0][0])):
        clustering_per_time_step=[]
        assortativity_degree_per_time_step=[]
        assortativity_age_per_time_step=[]
        assortativity_race_per_time_step=[]
        average_diameter_per_time_step=[]
        average_shortest_path_length_per_time_step=[]
        number_of_subgraphs_per_time_step=[]
        largest_component_per_time_step=[]
        modularity_per_time_step=[]
        betweenness_per_time_step=[]


        for graph in grouped_lists[timestep]:
            
            largest_cc = max(nx.connected_components(graph), key=len)
            S = graph.subgraph(largest_cc).copy()
            #connectivity.append(approx.node_connectivity(nx.Graph(all_adjacency_matrices[i])))

            clustering_coefficient=nx.transitivity(graph)
            assortativity_coefficient_degree=nx.degree_assortativity_coefficient(graph)
            assortativity_coefficient_race=nx.attribute_assortativity_coefficient(graph,'race')
            assortativity_coefficient_age=nx.numeric_assortativity_coefficient(graph,'age')
            
            list_sizes_connected_components=[len(c) for c in sorted(nx.connected_components(graph), key=len, reverse=True)]
            #print('list_sizes_connected_components:', list_sizes_connected_components)
            largest_component_size=list_sizes_connected_components[0]
            
            number_of_subgraphs=len(list_sizes_connected_components)
            
            average_diameter=nx.diameter(S)
            average_shortest_path_length=nx.average_shortest_path_length(S)

            partition = community.best_partition(S)

            modularity=community.modularity(partition,S)

            # Compute degree for each node
            degrees = dict(S.degree())

            # Sort nodes based on degree in descending order
            sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)

            # Compute the number of nodes to consider (top 10%)
            top_10_percent = int(len(sorted_nodes) * 0.1)

            # Select the top 10% nodes
            top_nodes = sorted_nodes[:top_10_percent]

            # Compute betweenness centrality for selected nodes
            betweenness = nx.betweenness_centrality_subset(S, sources=top_nodes, targets=top_nodes)

            # Compute average betweenness centrality
            average_betweenness = sum(betweenness.values()) / len(betweenness)
            

            clustering_per_time_step.append(clustering_coefficient)
            assortativity_degree_per_time_step.append(assortativity_coefficient_degree)
            assortativity_age_per_time_step.append(assortativity_coefficient_age)
            assortativity_race_per_time_step.append(assortativity_coefficient_race)
            average_diameter_per_time_step.append(average_diameter)
            average_shortest_path_length_per_time_step.append(average_shortest_path_length)
            largest_component_per_time_step.append(largest_component_size)
            number_of_subgraphs_per_time_step.append(number_of_subgraphs)
            modularity_per_time_step.append(modularity)
            betweenness_per_time_step.append(average_betweenness)


            counter+=1

        clustering_global.append(clustering_per_time_step)
        assortativity_degree_global.append(assortativity_degree_per_time_step)
        assortativity_age_global.append(assortativity_age_per_time_step)
        assortativity_race_global.append(assortativity_race_per_time_step)
        average_diameter_global.append(average_diameter_per_time_step)
        average_shortest_path_length_global.append(average_shortest_path_length_per_time_step)
        largest_component_size_global.append(largest_component_per_time_step)
        number_of_subgraphs_global.append(number_of_subgraphs_per_time_step)
        modularity_global.append(modularity_per_time_step)
        betweenness_global.append(betweenness_per_time_step)
        
       

        #get mean for time step and standard deviation

    clustering_global_summary=[]
    assortativity_degree_global_summary=[]
    assortativity_age_global_summary=[]
    assortativity_race_global_summary=[]
    average_diameter_global_summary=[]
    average_shortest_path_length_global_summary=[]
    largest_component_summary=[]
    number_of_subgraphs_global_summary=[]
    modularity_global_summary=[]
    betweenness_global_summary=[]
    


    for timestep in range(0,len(list_of_clustered_networks_16[0][0])):
            clustering_global_summary.append([np.mean(clustering_global[timestep]),np.std(clustering_global[timestep])])
            assortativity_degree_global_summary.append([np.mean(assortativity_degree_global[timestep]),np.std(assortativity_degree_global[timestep])])
            assortativity_age_global_summary.append([np.mean(assortativity_age_global[timestep]),np.std(assortativity_age_global[timestep])])
            assortativity_race_global_summary.append([np.mean(assortativity_race_global[timestep]),np.std(assortativity_race_global[timestep])])
            average_diameter_global_summary.append([np.mean(average_diameter_global[timestep]),np.std(average_diameter_global[timestep])])
            average_shortest_path_length_global_summary.append([np.mean(average_shortest_path_length_global[timestep]),np.std(average_shortest_path_length_global[timestep])])
            largest_component_summary.append([np.mean(largest_component_size_global[timestep]),np.std(largest_component_size_global[timestep])])
            number_of_subgraphs_global_summary.append([np.mean(number_of_subgraphs_global[timestep]),np.std(number_of_subgraphs_global[timestep])])
            modularity_global_summary.append([np.mean(modularity_global[timestep]),np.std(modularity_global[timestep])])
            betweenness_global_summary.append([np.mean(betweenness_global[timestep]),np.std(betweenness_global[timestep])])

    # Create multi-level column index for multiple columns
    columns = pd.MultiIndex.from_tuples([
        ('Clustering ', 'Mean'), ('Clustering ', 'Std'),
        ('Assortativity Degree Coefficient', 'Mean'), ('Assortativity Degree Coefficient', 'Std'),
        ('Assortativity Age Coefficient', 'Mean'), ('Assortativity Age Coefficient', 'Std'),
        ('Assortativity Race Coefficient', 'Mean'), ('Assortativity Race Coefficient', 'Std'),
        ('Average Diameter', 'Mean'), ('Average Diameter', 'Std'),
        ('Average Shortest Path Length', 'Mean'), ('Average Shortest Path Length', 'Std'),
        ('Largest Component Size', 'Mean'), ('Largest Component Size', 'Std'),
        ('Number of Subgraphs', 'Mean'), ('Number of Subgraphs', 'Std'),
        ('Modularity', 'Mean'), ('Modularity', 'Std'),
        ('Betweenness', 'Mean'), ('Betweenness', 'Std')
    ])

    time_step_0=flatten_comprehension([clustering_global_summary[0]+assortativity_degree_global_summary[0]+assortativity_age_global_summary[0]+assortativity_race_global_summary[0]\
                 +average_diameter_global_summary[0]+average_shortest_path_length_global_summary[0]+largest_component_summary[0]+number_of_subgraphs_global_summary[0]+\
                    modularity_global_summary[0]+betweenness_global_summary[0]])
    
    time_step_1=flatten_comprehension([clustering_global_summary[1]+assortativity_degree_global_summary[1]+assortativity_age_global_summary[1]+assortativity_race_global_summary[1]\
                 +average_diameter_global_summary[1]+average_shortest_path_length_global_summary[1]+largest_component_summary[1]+number_of_subgraphs_global_summary[1]+\
                    modularity_global_summary[1]+betweenness_global_summary[1]])
    
    time_step_2=flatten_comprehension([clustering_global_summary[2]+assortativity_degree_global_summary[2]+assortativity_age_global_summary[2]+assortativity_race_global_summary[2]\
                 +average_diameter_global_summary[2]+average_shortest_path_length_global_summary[2]+largest_component_summary[2]+number_of_subgraphs_global_summary[2]+\
                    modularity_global_summary[2]+betweenness_global_summary[2]])
    
    time_step_3=flatten_comprehension([clustering_global_summary[3]+assortativity_degree_global_summary[3]+assortativity_age_global_summary[3]+assortativity_race_global_summary[3]\
                 +average_diameter_global_summary[3]+average_shortest_path_length_global_summary[3]+largest_component_summary[3]+number_of_subgraphs_global_summary[3]+\
                    modularity_global_summary[3]+betweenness_global_summary[3]])
    
    time_step_4=flatten_comprehension([clustering_global_summary[4]+assortativity_degree_global_summary[4]+assortativity_age_global_summary[4]+assortativity_race_global_summary[4]\
                 +average_diameter_global_summary[4]+average_shortest_path_length_global_summary[4]+largest_component_summary[4]+number_of_subgraphs_global_summary[4]+\
                    modularity_global_summary[4]+betweenness_global_summary[4]])
    
    time_step_5=flatten_comprehension([clustering_global_summary[5]+assortativity_degree_global_summary[5]+assortativity_age_global_summary[5]+assortativity_race_global_summary[5]\
                 +average_diameter_global_summary[5]+average_shortest_path_length_global_summary[5]+largest_component_summary[5]+number_of_subgraphs_global_summary[5]+\
                    modularity_global_summary[5]+betweenness_global_summary[5]])
    
    # Sample data
    data = [time_step_0,time_step_1,time_step_2,time_step_3,time_step_4,time_step_5]

    # Create a DataFrame with multi-level columns
    dataframe_metrics_clustering_effect = pd.DataFrame(data, columns=columns)

    file_name = "dataframe_metrics_clustering_effect_new_with_modularity_and_betweenness.pkl"

    open_file = open(file_name, "wb")
    pickle.dump(dataframe_metrics_clustering_effect, open_file)
    open_file.close()


if property_under_study=='Assortativity By Degree':

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'agents_complete_final_use_in_analysis.pkl'), 'rb') #rb = read bytes because we are reading the file
    agents = pickle.load(file)

    #Assortativity By Degree
        
    list_of_assortativity_by_degree_15=[]
    id_networks=list(range(0,7+1,1))
    directory_path_assortativity_by_degree_15 = "/home/jbrazia/metrics/Network_15/Assortativity/Degree"

    for id in id_networks:
        file_name_assortativity_by_degree_15 = "assortativity_degree_500000_iterations_15_stochastic_ID_"+str(id)+".pkl"
        file_path_assortativity_by_degree_15 = os.path.join(directory_path_assortativity_by_degree_15, file_name_assortativity_by_degree_15)

        with open(file_path_assortativity_by_degree_15, "rb") as open_file_assortativity_by_degree_15:
            stationary_networks_assortativity_by_degree_15 = pickle.load(open_file_assortativity_by_degree_15)
            list_of_assortativity_by_degree_15.append(stationary_networks_assortativity_by_degree_15)

    list_of_assortativity_by_degree_16=[]
    id_networks=list(range(0,7+1,1))
    directory_path_assortativity_by_degree_16 = "/home/jbrazia/metrics/Network_16/Assortativity/Degree"

    for id in id_networks:
        file_name_assortativity_by_degree_16 = "assortativity_degree_500000_iterations_16_stochastic_ID_"+str(id)+".pkl"
        file_path_assortativity_by_degree_16 = os.path.join(directory_path_assortativity_by_degree_16, file_name_assortativity_by_degree_16)

        with open(file_path_assortativity_by_degree_16, "rb") as open_file_assortativity_by_degree_16:
            stationary_networks_assortativity_by_degree_16 = pickle.load(open_file_assortativity_by_degree_16)
            list_of_assortativity_by_degree_16.append(stationary_networks_assortativity_by_degree_16)

    list_of_assortativity_by_degree_17=[]
    id_networks=list(range(0,7+1,1))
    directory_path_assortativity_by_degree_17 = "/home/jbrazia/metrics/Network_17/Assortativity/Degree"

    for id in id_networks:
        file_name_assortativity_by_degree_17 = "assortativity_degree_500000_iterations_17_stochastic_ID_"+str(id)+".pkl"
        file_path_assortativity_by_degree_17 = os.path.join(directory_path_assortativity_by_degree_17, file_name_assortativity_by_degree_17)

        with open(file_path_assortativity_by_degree_17, "rb") as open_file_assortativity_by_degree_17:
            stationary_networks_assortativity_by_degree_17 = pickle.load(open_file_assortativity_by_degree_17)
            list_of_assortativity_by_degree_17.append(stationary_networks_assortativity_by_degree_17)

    list_of_all_assortativity_by_degree_networks=[list_of_assortativity_by_degree_15[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_assortativity_by_degree_15))]+[list_of_assortativity_by_degree_16[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_assortativity_by_degree_16))]+[list_of_assortativity_by_degree_17[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_assortativity_by_degree_17))]
    
    flattened_list_of_assortativity_by_degree_networks = flatten_comprehension(list_of_all_assortativity_by_degree_networks)
    

    for graph in flattened_list_of_assortativity_by_degree_networks:
        graph=graph_with_atributes(graph,agents)
    #group lists according to timesteps
    
    grouped_lists=[]
    for timestep in range(0,len(list_of_assortativity_by_degree_15[0][0])):
        graph_ID=timestep
        group_per_time_step=[]
        while graph_ID<len(flattened_list_of_assortativity_by_degree_networks):
            group_per_time_step.append(flattened_list_of_assortativity_by_degree_networks[graph_ID])
            graph_ID+=6
        grouped_lists.append(group_per_time_step)

    #global properties per generated graph, we compute and compare the properties of the 
    counter=0
    clustering_global=[]
    assortativity_degree_global=[]
    assortativity_age_global=[]
    assortativity_race_global=[]
    average_diameter_global=[]
    average_shortest_path_length_global=[]
    largest_component_size_global=[]
    number_of_subgraphs_global=[]
    modularity_global=[]
    betweenness_global=[]

    for timestep in range(0,len(list_of_assortativity_by_degree_15[0][0])):
        clustering_per_time_step=[]
        assortativity_degree_per_time_step=[]
        assortativity_age_per_time_step=[]
        assortativity_race_per_time_step=[]
        average_diameter_per_time_step=[]
        average_shortest_path_length_per_time_step=[]
        number_of_subgraphs_per_time_step=[]
        largest_component_per_time_step=[]
        modularity_per_time_step=[]
        betweenness_per_time_step=[]


        for graph in grouped_lists[timestep]:
            
            largest_cc = max(nx.connected_components(graph), key=len)
            S = graph.subgraph(largest_cc).copy()
            #connectivity.append(approx.node_connectivity(nx.Graph(all_adjacency_matrices[i])))

            
            assortativity_coefficient_degree=nx.degree_assortativity_coefficient(graph)
            clustering_coefficient=nx.transitivity(graph)
            assortativity_coefficient_race=nx.attribute_assortativity_coefficient(graph,'race')
            assortativity_coefficient_age=nx.numeric_assortativity_coefficient(graph,'age')
            
            list_sizes_connected_components=[len(c) for c in sorted(nx.connected_components(graph), key=len, reverse=True)]
            #print('list_sizes_connected_components:', list_sizes_connected_components)
            largest_component_size=list_sizes_connected_components[0]
            
            number_of_subgraphs=len(list_sizes_connected_components)
            
            average_diameter=nx.diameter(S)
            average_shortest_path_length=nx.average_shortest_path_length(S)

            partition = community.best_partition(S)

            modularity=community.modularity(partition,S)

            # Compute degree for each node
            degrees = dict(S.degree())

            # Sort nodes based on degree in descending order
            sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)

            # Compute the number of nodes to consider (top 10%)
            top_10_percent = int(len(sorted_nodes) * 0.1)

            # Select the top 10% nodes
            top_nodes = sorted_nodes[:top_10_percent]

            # Compute betweenness centrality for selected nodes
            betweenness = nx.betweenness_centrality_subset(S, sources=top_nodes, targets=top_nodes)

            # Compute average betweenness centrality
            average_betweenness = sum(betweenness.values()) / len(betweenness)
            
            
            assortativity_degree_per_time_step.append(assortativity_coefficient_degree)
            clustering_per_time_step.append(clustering_coefficient)
            assortativity_age_per_time_step.append(assortativity_coefficient_age)
            assortativity_race_per_time_step.append(assortativity_coefficient_race)
            average_diameter_per_time_step.append(average_diameter)
            average_shortest_path_length_per_time_step.append(average_shortest_path_length)
            largest_component_per_time_step.append(largest_component_size)
            number_of_subgraphs_per_time_step.append(number_of_subgraphs)
            modularity_per_time_step.append(modularity)
            betweenness_per_time_step.append(average_betweenness)
            
            counter+=1
            print('counter:',counter)

        
        assortativity_degree_global.append(assortativity_degree_per_time_step)
        clustering_global.append(clustering_per_time_step)
        assortativity_age_global.append(assortativity_age_per_time_step)
        assortativity_race_global.append(assortativity_race_per_time_step)
        average_diameter_global.append(average_diameter_per_time_step)
        average_shortest_path_length_global.append(average_shortest_path_length_per_time_step)
        largest_component_size_global.append(largest_component_per_time_step)
        number_of_subgraphs_global.append(number_of_subgraphs_per_time_step)
        modularity_global.append(modularity_per_time_step)
        betweenness_global.append(betweenness_per_time_step)
        

    clustering_global_summary=[]
    assortativity_degree_global_summary=[]
    assortativity_age_global_summary=[]
    assortativity_race_global_summary=[]
    average_diameter_global_summary=[]
    average_shortest_path_length_global_summary=[]
    largest_component_summary=[]
    number_of_subgraphs_global_summary=[]
    modularity_global_summary=[]
    betweenness_global_summary=[]
    

    for timestep in range(0,len(list_of_assortativity_by_degree_15[0][0])):
            assortativity_degree_global_summary.append([np.mean(assortativity_degree_global[timestep]),np.std(assortativity_degree_global[timestep])])
            clustering_global_summary.append([np.mean(clustering_global[timestep]),np.std(clustering_global[timestep])])
            assortativity_age_global_summary.append([np.mean(assortativity_age_global[timestep]),np.std(assortativity_age_global[timestep])])
            assortativity_race_global_summary.append([np.mean(assortativity_race_global[timestep]),np.std(assortativity_race_global[timestep])])
            average_diameter_global_summary.append([np.mean(average_diameter_global[timestep]),np.std(average_diameter_global[timestep])])
            average_shortest_path_length_global_summary.append([np.mean(average_shortest_path_length_global[timestep]),np.std(average_shortest_path_length_global[timestep])])
            largest_component_summary.append([np.mean(largest_component_size_global[timestep]),np.std(largest_component_size_global[timestep])])
            number_of_subgraphs_global_summary.append([np.mean(number_of_subgraphs_global[timestep]),np.std(number_of_subgraphs_global[timestep])])
            modularity_global_summary.append([np.mean(modularity_global[timestep]),np.std(modularity_global[timestep])])
            betweenness_global_summary.append([np.mean(betweenness_global[timestep]),np.std(betweenness_global[timestep])])
            
    
    # Create multi-level column index for multiple columns
    columns = pd.MultiIndex.from_tuples([
        ('Assortativity Degree Coefficient', 'Mean'), ('Assortativity Degree Coefficient', 'Std'),
        ('Clustering ', 'Mean'), ('Clustering ', 'Std'),
        ('Assortativity Age Coefficient', 'Mean'), ('Assortativity Age Coefficient', 'Std'),
        ('Assortativity Race Coefficient', 'Mean'), ('Assortativity Race Coefficient', 'Std'),
        ('Average Diameter', 'Mean'), ('Average Diameter', 'Std'),
        ('Average Shortest Path Length', 'Mean'), ('Average Shortest Path Length', 'Std'),
        ('Largest Component Size', 'Mean'), ('Largest Component Size', 'Std'),
        ('Number of Subgraphs', 'Mean'), ('Number of Subgraphs', 'Std'),
        ('Modularity', 'Mean'), ('Modularity', 'Std'),
        ('Betweenness', 'Mean'), ('Betweenness', 'Std')
    ])

    time_step_0=flatten_comprehension([assortativity_degree_global_summary[0]+clustering_global_summary[0]+assortativity_age_global_summary[0]+assortativity_race_global_summary[0]\
                 +average_diameter_global_summary[0]+average_shortest_path_length_global_summary[0]+largest_component_summary[0]+number_of_subgraphs_global_summary[0]\
                    +modularity_global_summary[0]+betweenness_global_summary[0]])
    
    time_step_1=flatten_comprehension([assortativity_degree_global_summary[1]+clustering_global_summary[1]+assortativity_age_global_summary[1]+assortativity_race_global_summary[1]\
                 +average_diameter_global_summary[1]+average_shortest_path_length_global_summary[1]+largest_component_summary[1]+number_of_subgraphs_global_summary[1]\
                    +modularity_global_summary[1]+betweenness_global_summary[1]])
    
    time_step_2=flatten_comprehension([assortativity_degree_global_summary[2]+clustering_global_summary[2]+assortativity_age_global_summary[2]+assortativity_race_global_summary[2]\
                 +average_diameter_global_summary[2]+average_shortest_path_length_global_summary[2]+largest_component_summary[2]+number_of_subgraphs_global_summary[2]\
                    +modularity_global_summary[2]+betweenness_global_summary[2]])
    
    time_step_3=flatten_comprehension([assortativity_degree_global_summary[3]+clustering_global_summary[3]+assortativity_age_global_summary[3]+assortativity_race_global_summary[3]\
                 +average_diameter_global_summary[3]+average_shortest_path_length_global_summary[3]+largest_component_summary[3]+number_of_subgraphs_global_summary[3]\
                    +modularity_global_summary[3]+betweenness_global_summary[3]])
    
    time_step_4=flatten_comprehension([assortativity_degree_global_summary[4]+clustering_global_summary[4]+assortativity_age_global_summary[4]+assortativity_race_global_summary[4]\
                 +average_diameter_global_summary[4]+average_shortest_path_length_global_summary[4]+largest_component_summary[4]+number_of_subgraphs_global_summary[4]\
                    +modularity_global_summary[4]+betweenness_global_summary[4]])
    
    time_step_5=flatten_comprehension([assortativity_degree_global_summary[5]+clustering_global_summary[5]+assortativity_age_global_summary[5]+assortativity_race_global_summary[5]\
                 +average_diameter_global_summary[5]+average_shortest_path_length_global_summary[5]+largest_component_summary[5]+number_of_subgraphs_global_summary[5]\
                    +modularity_global_summary[5]+betweenness_global_summary[5]])

    # Sample data
    data = [time_step_0,time_step_1,time_step_2,time_step_3,time_step_4,time_step_5]

    # Create a DataFrame with multi-level columns
    dataframe_metrics_assortativity_degree_effect = pd.DataFrame(data, columns=columns)

    file_name = "dataframe_metrics_assortativity_degree_effect_new_with_modularity_and_betweenness.pkl"

    open_file = open(file_name, "wb")
    pickle.dump(dataframe_metrics_assortativity_degree_effect, open_file)
    open_file.close()


if property_under_study=='Assortativity By Age':

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'agents_complete_final_use_in_analysis.pkl'), 'rb') #rb = read bytes because we are reading the file
    agents = pickle.load(file)

    #Assortativity By Age
        
    list_of_assortativity_by_age_15=[]
    id_networks=list(range(0,7+1,1))
    directory_path_assortativity_by_age_15 = "/home/jbrazia/metrics/Network_15/Assortativity/Age"

    for id in id_networks:
        file_name_assortativity_by_age_15 = "assortativity_age_network_500000_iterations_15_stochastic_ID_"+str(id)+".pkl"
        file_path_assortativity_by_age_15 = os.path.join(directory_path_assortativity_by_age_15, file_name_assortativity_by_age_15)

        with open(file_path_assortativity_by_age_15, "rb") as open_file_assortativity_by_age_15:
            stationary_networks_assortativity_by_age_15 = pickle.load(open_file_assortativity_by_age_15)
            list_of_assortativity_by_age_15.append(stationary_networks_assortativity_by_age_15)

    list_of_assortativity_by_age_16=[]
    id_networks=list(range(0,7+1,1))
    directory_path_assortativity_by_age_16 = "/home/jbrazia/metrics/Network_16/Assortativity/Age"

    for id in id_networks:
        file_name_assortativity_by_age_16 = "assortativity_age_network_500000_iterations_16_stochastic_ID_"+str(id)+".pkl"
        file_path_assortativity_by_age_16 = os.path.join(directory_path_assortativity_by_age_16, file_name_assortativity_by_age_16)

        with open(file_path_assortativity_by_age_16, "rb") as open_file_assortativity_by_age_16:
            stationary_networks_assortativity_by_age_16 = pickle.load(open_file_assortativity_by_age_16)
            list_of_assortativity_by_age_16.append(stationary_networks_assortativity_by_age_16)

    list_of_assortativity_by_age_17=[]
    id_networks=list(range(0,7+1,1))
    directory_path_assortativity_by_age_17 = "/home/jbrazia/metrics/Network_17/Assortativity/Age"

    for id in id_networks:
        file_name_assortativity_by_age_17 = "assortativity_age_network_500000_iterations_17_stochastic_ID_"+str(id)+".pkl"
        file_path_assortativity_by_age_17 = os.path.join(directory_path_assortativity_by_age_17, file_name_assortativity_by_age_17)

        with open(file_path_assortativity_by_age_17, "rb") as open_file_assortativity_by_age_17:
            stationary_networks_assortativity_by_age_17 = pickle.load(open_file_assortativity_by_age_17)
            list_of_assortativity_by_age_17.append(stationary_networks_assortativity_by_age_17)

    list_of_all_assortativity_by_age_networks=[list_of_assortativity_by_age_15[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_assortativity_by_age_15))]+[list_of_assortativity_by_age_16[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_assortativity_by_age_16))]+[list_of_assortativity_by_age_17[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_assortativity_by_age_17))]
    
    flattened_list_of_assortativity_by_age_networks = flatten_comprehension(list_of_all_assortativity_by_age_networks)

    for graph in flattened_list_of_assortativity_by_age_networks:
        graph=graph_with_atributes(graph,agents)

    #group lists according to timesteps
    
    grouped_lists=[]
    for timestep in range(0,len(list_of_assortativity_by_age_15[0][0])):
        graph_ID=timestep
        group_per_time_step=[]
        while graph_ID<len(flattened_list_of_assortativity_by_age_networks):
            group_per_time_step.append(flattened_list_of_assortativity_by_age_networks[graph_ID])
            graph_ID+=6
        grouped_lists.append(group_per_time_step)

    #global properties per generated graph, we compute and compare the properties of the 
    counter=0
    clustering_global=[]
    assortativity_degree_global=[]
    assortativity_age_global=[]
    assortativity_race_global=[]
    average_diameter_global=[]
    average_shortest_path_length_global=[]
    largest_component_size_global=[]
    number_of_subgraphs_global=[]
    modularity_global=[]
    betweenness_global=[]


    for timestep in range(0,len(list_of_assortativity_by_age_15[0][0])):
        clustering_per_time_step=[]
        assortativity_degree_per_time_step=[]
        assortativity_age_per_time_step=[]
        assortativity_race_per_time_step=[]
        average_diameter_per_time_step=[]
        average_shortest_path_length_per_time_step=[]
        number_of_subgraphs_per_time_step=[]
        largest_component_per_time_step=[]
        modularity_per_time_step=[]
        betweenness_per_time_step=[]
        

        for graph in grouped_lists[timestep]:
            
            largest_cc = max(nx.connected_components(graph), key=len)
            S = graph.subgraph(largest_cc).copy()
            #connectivity.append(approx.node_connectivity(nx.Graph(all_adjacency_matrices[i])))

            
            assortativity_coefficient_degree=nx.degree_assortativity_coefficient(graph)
            clustering_coefficient=nx.transitivity(graph)
            assortativity_coefficient_race=nx.attribute_assortativity_coefficient(graph,'race')
            assortativity_coefficient_age=nx.numeric_assortativity_coefficient(graph,'age')
            
            list_sizes_connected_components=[len(c) for c in sorted(nx.connected_components(graph), key=len, reverse=True)]
            #print('list_sizes_connected_components:', list_sizes_connected_components)
            largest_component_size=list_sizes_connected_components[0]
            
            number_of_subgraphs=len(list_sizes_connected_components)
            
            average_diameter=nx.diameter(S)
            average_shortest_path_length=nx.average_shortest_path_length(S)

            

            partition = community.best_partition(S)

            modularity=community.modularity(partition,S)

            # Compute degree for each node
            degrees = dict(S.degree())

            # Sort nodes based on degree in descending order
            sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)

            # Compute the number of nodes to consider (top 10%)
            top_10_percent = int(len(sorted_nodes) * 0.1)

            # Select the top 10% nodes
            top_nodes = sorted_nodes[:top_10_percent]

            # Compute betweenness centrality for selected nodes
            betweenness = nx.betweenness_centrality_subset(S, sources=top_nodes, targets=top_nodes)

            # Compute average betweenness centrality
            average_betweenness = sum(betweenness.values()) / len(betweenness)
            
            
            assortativity_age_per_time_step.append(assortativity_coefficient_age)
            assortativity_degree_per_time_step.append(assortativity_coefficient_degree)
            clustering_per_time_step.append(clustering_coefficient)
            assortativity_race_per_time_step.append(assortativity_coefficient_race)
            average_diameter_per_time_step.append(average_diameter)
            average_shortest_path_length_per_time_step.append(average_shortest_path_length)
            largest_component_per_time_step.append(largest_component_size)
            number_of_subgraphs_per_time_step.append(number_of_subgraphs)
            modularity_per_time_step.append(modularity)
            betweenness_per_time_step.append(average_betweenness)

            
            counter+=1
            print('counter:',counter)

        assortativity_age_global.append(assortativity_age_per_time_step)
        assortativity_degree_global.append(assortativity_degree_per_time_step)
        clustering_global.append(clustering_per_time_step)
        assortativity_race_global.append(assortativity_race_per_time_step)
        average_diameter_global.append(average_diameter_per_time_step)
        average_shortest_path_length_global.append(average_shortest_path_length_per_time_step)
        largest_component_size_global.append(largest_component_per_time_step)
        number_of_subgraphs_global.append(number_of_subgraphs_per_time_step)
        modularity_global.append(modularity_per_time_step)
        betweenness_global.append(betweenness_per_time_step)


    clustering_global_summary=[]
    assortativity_degree_global_summary=[]
    assortativity_age_global_summary=[]
    assortativity_race_global_summary=[]
    average_diameter_global_summary=[]
    average_shortest_path_length_global_summary=[]
    largest_component_summary=[]
    number_of_subgraphs_global_summary=[]
    modularity_global_summary=[]
    betweenness_global_summary=[]
    

    for timestep in range(0,len(list_of_assortativity_by_age_15[0][0])):
            assortativity_degree_global_summary.append([np.mean(assortativity_degree_global[timestep]),np.std(assortativity_degree_global[timestep])])
            clustering_global_summary.append([np.mean(clustering_global[timestep]),np.std(clustering_global[timestep])])
            assortativity_age_global_summary.append([np.mean(assortativity_age_global[timestep]),np.std(assortativity_age_global[timestep])])
            assortativity_race_global_summary.append([np.mean(assortativity_race_global[timestep]),np.std(assortativity_race_global[timestep])])
            average_diameter_global_summary.append([np.mean(average_diameter_global[timestep]),np.std(average_diameter_global[timestep])])
            average_shortest_path_length_global_summary.append([np.mean(average_shortest_path_length_global[timestep]),np.std(average_shortest_path_length_global[timestep])])
            largest_component_summary.append([np.mean(largest_component_size_global[timestep]),np.std(largest_component_size_global[timestep])])
            number_of_subgraphs_global_summary.append([np.mean(number_of_subgraphs_global[timestep]),np.std(number_of_subgraphs_global[timestep])])
            modularity_global_summary.append([np.mean(modularity_global[timestep]),np.std(modularity_global[timestep])])
            betweenness_global_summary.append([np.mean(betweenness_global[timestep]),np.std(betweenness_global[timestep])])
            
    
    # Create multi-level column index for multiple columns
    columns = pd.MultiIndex.from_tuples([
        ('Assortativity Age Coefficient', 'Mean'), ('Assortativity Age Coefficient', 'Std'),
        ('Assortativity Degree Coefficient', 'Mean'), ('Assortativity Degree Coefficient', 'Std'),
        ('Clustering ', 'Mean'), ('Clustering ', 'Std'),
        ('Assortativity Race Coefficient', 'Mean'), ('Assortativity Race Coefficient', 'Std'),
        ('Average Diameter', 'Mean'), ('Average Diameter', 'Std'),
        ('Average Shortest Path Length', 'Mean'), ('Average Shortest Path Length', 'Std'),
        ('Largest Component Size', 'Mean'), ('Largest Component Size', 'Std'),
        ('Number of Subgraphs', 'Mean'), ('Number of Subgraphs', 'Std'),
        ('Modularity', 'Mean'), ('Modularity', 'Std'),
        ('Betweenness', 'Mean'), ('Betweenness', 'Std')
    ])

    time_step_0=flatten_comprehension([assortativity_age_global_summary[0]+assortativity_degree_global_summary[0]+clustering_global_summary[0]+assortativity_race_global_summary[0]\
                 +average_diameter_global_summary[0]+average_shortest_path_length_global_summary[0]+largest_component_summary[0]+number_of_subgraphs_global_summary[0]\
                    +modularity_global_summary[0]+betweenness_global_summary[0]])
    
    time_step_1=flatten_comprehension([assortativity_age_global_summary[1]+assortativity_degree_global_summary[1]+clustering_global_summary[1]+assortativity_race_global_summary[1]\
                 +average_diameter_global_summary[1]+average_shortest_path_length_global_summary[1]+largest_component_summary[1]+number_of_subgraphs_global_summary[1]\
                    +modularity_global_summary[1]+betweenness_global_summary[1]])
    
    time_step_2=flatten_comprehension([assortativity_age_global_summary[2]+assortativity_degree_global_summary[2]+clustering_global_summary[2]+assortativity_race_global_summary[2]\
                 +average_diameter_global_summary[2]+average_shortest_path_length_global_summary[2]+largest_component_summary[2]+number_of_subgraphs_global_summary[2]\
                    +modularity_global_summary[2]+betweenness_global_summary[2]])
    
    time_step_3=flatten_comprehension([assortativity_age_global_summary[3]+assortativity_degree_global_summary[3]+clustering_global_summary[3]+assortativity_race_global_summary[3]\
                 +average_diameter_global_summary[3]+average_shortest_path_length_global_summary[3]+largest_component_summary[3]+number_of_subgraphs_global_summary[3]\
                    +modularity_global_summary[3]+betweenness_global_summary[3]])
    
    time_step_4=flatten_comprehension([assortativity_age_global_summary[4]+assortativity_degree_global_summary[4]+clustering_global_summary[4]+assortativity_race_global_summary[4]\
                 +average_diameter_global_summary[4]+average_shortest_path_length_global_summary[4]+largest_component_summary[4]+number_of_subgraphs_global_summary[4]\
                    +modularity_global_summary[4]+betweenness_global_summary[4]])
    
    time_step_5=flatten_comprehension([assortativity_age_global_summary[5]+assortativity_degree_global_summary[5]+clustering_global_summary[5]+assortativity_race_global_summary[5]\
                 +average_diameter_global_summary[5]+average_shortest_path_length_global_summary[5]+largest_component_summary[5]+number_of_subgraphs_global_summary[5]\
                    +modularity_global_summary[5]+betweenness_global_summary[5]])

    # Sample data
    data = [time_step_0,time_step_1,time_step_2,time_step_3,time_step_4,time_step_5]

    # Create a DataFrame with multi-level columns
    dataframe_metrics_assortativity_age_effect = pd.DataFrame(data, columns=columns)

    file_name = "dataframe_metrics_assortativity_age_effect_new_with_modularity_and_betweenness.pkl"

    open_file = open(file_name, "wb")
    pickle.dump(dataframe_metrics_assortativity_age_effect, open_file)
    open_file.close()
    


if property_under_study=='Assortativity By Race':

    #Assortativity By Race

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'agents_complete_final_use_in_analysis.pkl'), 'rb') #rb = read bytes because we are reading the file
    agents = pickle.load(file)


    list_of_assortativity_by_race_15=[]
    id_networks=list(range(0,7+1,1))
    directory_path_assortativity_by_race_15 = "/home/jbrazia/metrics/Network_15/Assortativity/Race"

    for id in id_networks:
        file_name_assortativity_by_race_15 = "assortativity_race_network_500000_iterations_15_stochastic_ID_"+str(id)+".pkl"
        file_path_assortativity_by_race_15 = os.path.join(directory_path_assortativity_by_race_15, file_name_assortativity_by_race_15)

        with open(file_path_assortativity_by_race_15, "rb") as open_file_assortativity_by_race_15:
            stationary_networks_assortativity_by_race_15 = pickle.load(open_file_assortativity_by_race_15)
            list_of_assortativity_by_race_15.append(stationary_networks_assortativity_by_race_15)

    list_of_assortativity_by_race_16=[]
    id_networks=list(range(0,7+1,1))
    directory_path_assortativity_by_race_16 = "/home/jbrazia/metrics/Network_16/Assortativity/Race"

    for id in id_networks:
        file_name_assortativity_by_race_16 = "assortativity_race_network_500000_iterations_16_stochastic_ID_"+str(id)+".pkl"
        file_path_assortativity_by_race_16 = os.path.join(directory_path_assortativity_by_race_16, file_name_assortativity_by_race_16)

        with open(file_path_assortativity_by_race_16, "rb") as open_file_assortativity_by_race_16:
            stationary_networks_assortativity_by_race_16 = pickle.load(open_file_assortativity_by_race_16)
            list_of_assortativity_by_race_16.append(stationary_networks_assortativity_by_race_16)

    list_of_assortativity_by_race_17=[]
    id_networks=list(range(0,7+1,1))
    directory_path_assortativity_by_race_17 = "/home/jbrazia/metrics/Network_17/Assortativity/Race"

    for id in id_networks:
        file_name_assortativity_by_race_17 = "assortativity_race_network_500000_iterations_17_stochastic_ID_"+str(id)+".pkl"
        file_path_assortativity_by_race_17 = os.path.join(directory_path_assortativity_by_race_17, file_name_assortativity_by_race_17)

        with open(file_path_assortativity_by_race_17, "rb") as open_file_assortativity_by_race_17:
            stationary_networks_assortativity_by_race_17 = pickle.load(open_file_assortativity_by_race_17)
            list_of_assortativity_by_race_17.append(stationary_networks_assortativity_by_race_17)
    
    list_of_all_assortativity_by_race_networks=[list_of_assortativity_by_race_15[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_assortativity_by_race_15))]+[list_of_assortativity_by_race_16[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_assortativity_by_race_16))]+[list_of_assortativity_by_race_17[stochastic_network][0] for stochastic_network in \
                        range(0,len(list_of_assortativity_by_race_17))]
    
    flattened_list_of_assortativity_by_race_networks = flatten_comprehension(list_of_all_assortativity_by_race_networks)

    for graph in flattened_list_of_assortativity_by_race_networks:
        graph=graph_with_atributes(graph,agents)

    #group lists according to timesteps
    
    grouped_lists=[]
    for timestep in range(0,len(list_of_assortativity_by_race_15[0][0])):
        graph_ID=timestep
        group_per_time_step=[]
        while graph_ID<len(flattened_list_of_assortativity_by_race_networks):
            group_per_time_step.append(flattened_list_of_assortativity_by_race_networks[graph_ID])
            graph_ID+=6
        grouped_lists.append(group_per_time_step)

    #global properties per generated graph, we compute and compare the properties of the 
    counter=0
    clustering_global=[]
    assortativity_degree_global=[]
    assortativity_age_global=[]
    assortativity_race_global=[]
    average_diameter_global=[]
    average_shortest_path_length_global=[]
    largest_component_size_global=[]
    number_of_subgraphs_global=[]
    modularity_global=[]
    betweenness_global=[]

    for timestep in range(0,len(list_of_assortativity_by_race_15[0][0])):
        clustering_per_time_step=[]
        assortativity_degree_per_time_step=[]
        assortativity_age_per_time_step=[]
        assortativity_race_per_time_step=[]
        average_diameter_per_time_step=[]
        average_shortest_path_length_per_time_step=[]
        number_of_subgraphs_per_time_step=[]
        largest_component_per_time_step=[]
        modularity_per_time_step=[]
        betweenness_per_time_step=[]

        for graph in grouped_lists[timestep]:
            
            largest_cc = max(nx.connected_components(graph), key=len)
            S = graph.subgraph(largest_cc).copy()
            #connectivity.append(approx.node_connectivity(nx.Graph(all_adjacency_matrices[i])))

            
            assortativity_coefficient_degree=nx.degree_assortativity_coefficient(graph)
            clustering_coefficient=nx.transitivity(graph)
            assortativity_coefficient_race=nx.attribute_assortativity_coefficient(graph,'race')
            assortativity_coefficient_age=nx.numeric_assortativity_coefficient(graph,'age')
            
            list_sizes_connected_components=[len(c) for c in sorted(nx.connected_components(graph), key=len, reverse=True)]
            #print('list_sizes_connected_components:', list_sizes_connected_components)
            largest_component_size=list_sizes_connected_components[0]
            
            number_of_subgraphs=len(list_sizes_connected_components)
            
            average_diameter=nx.diameter(S)
            average_shortest_path_length=nx.average_shortest_path_length(S)

            partition = community.best_partition(S)

            modularity=community.modularity(partition,S)

            # Compute degree for each node
            degrees = dict(S.degree())

            # Sort nodes based on degree in descending order
            sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)

            # Compute the number of nodes to consider (top 10%)
            top_10_percent = int(len(sorted_nodes) * 0.1)

            # Select the top 10% nodes
            top_nodes = sorted_nodes[:top_10_percent]

            # Compute betweenness centrality for selected nodes
            betweenness = nx.betweenness_centrality_subset(S, sources=top_nodes, targets=top_nodes)

            # Compute average betweenness centrality
            average_betweenness = sum(betweenness.values()) / len(betweenness)
            
            

            assortativity_age_per_time_step.append(assortativity_coefficient_age)
            assortativity_degree_per_time_step.append(assortativity_coefficient_degree)
            clustering_per_time_step.append(clustering_coefficient)
            assortativity_race_per_time_step.append(assortativity_coefficient_race)
            average_diameter_per_time_step.append(average_diameter)
            average_shortest_path_length_per_time_step.append(average_shortest_path_length)
            largest_component_per_time_step.append(largest_component_size)
            number_of_subgraphs_per_time_step.append(number_of_subgraphs)
            modularity_per_time_step.append(modularity)
            betweenness_per_time_step.append(average_betweenness)
            
            counter+=1
            print('counter:',counter)

        assortativity_age_global.append(assortativity_age_per_time_step)
        assortativity_degree_global.append(assortativity_degree_per_time_step)
        clustering_global.append(clustering_per_time_step)
        assortativity_race_global.append(assortativity_race_per_time_step)
        average_diameter_global.append(average_diameter_per_time_step)
        average_shortest_path_length_global.append(average_shortest_path_length_per_time_step)
        largest_component_size_global.append(largest_component_per_time_step)
        number_of_subgraphs_global.append(number_of_subgraphs_per_time_step)
        modularity_global.append(modularity_per_time_step)
        betweenness_global.append(betweenness_per_time_step)

        

    clustering_global_summary=[]
    assortativity_degree_global_summary=[]
    assortativity_age_global_summary=[]
    assortativity_race_global_summary=[]
    average_diameter_global_summary=[]
    average_shortest_path_length_global_summary=[]
    largest_component_summary=[]
    number_of_subgraphs_global_summary=[]
    modularity_global_summary=[]
    betweenness_global_summary=[]

    for timestep in range(0,len(list_of_assortativity_by_race_15[0][0])):
            assortativity_degree_global_summary.append([np.mean(assortativity_degree_global[timestep]),np.std(assortativity_degree_global[timestep])])
            clustering_global_summary.append([np.mean(clustering_global[timestep]),np.std(clustering_global[timestep])])
            assortativity_age_global_summary.append([np.mean(assortativity_age_global[timestep]),np.std(assortativity_age_global[timestep])])
            assortativity_race_global_summary.append([np.mean(assortativity_race_global[timestep]),np.std(assortativity_race_global[timestep])])
            average_diameter_global_summary.append([np.mean(average_diameter_global[timestep]),np.std(average_diameter_global[timestep])])
            average_shortest_path_length_global_summary.append([np.mean(average_shortest_path_length_global[timestep]),np.std(average_shortest_path_length_global[timestep])])
            largest_component_summary.append([np.mean(largest_component_size_global[timestep]),np.std(largest_component_size_global[timestep])])
            number_of_subgraphs_global_summary.append([np.mean(number_of_subgraphs_global[timestep]),np.std(number_of_subgraphs_global[timestep])])
            modularity_global_summary.append([np.mean(modularity_global[timestep]),np.std(modularity_global[timestep])])
            betweenness_global_summary.append([np.mean(betweenness_global[timestep]),np.std(betweenness_global[timestep])])
    
    # Create multi-level column index for multiple columns
    columns = pd.MultiIndex.from_tuples([
        ('Assortativity Race Coefficient', 'Mean'), ('Assortativity Race Coefficient', 'Std'),
        ('Assortativity Age Coefficient', 'Mean'), ('Assortativity Age Coefficient', 'Std'),
        ('Assortativity Degree Coefficient', 'Mean'), ('Assortativity Degree Coefficient', 'Std'),
        ('Clustering ', 'Mean'), ('Clustering ', 'Std'),
        ('Average Diameter', 'Mean'), ('Average Diameter', 'Std'),
        ('Average Shortest Path Length', 'Mean'), ('Average Shortest Path Length', 'Std'),
        ('Largest Component Size', 'Mean'), ('Largest Component Size', 'Std'),
        ('Number of Subgraphs', 'Mean'), ('Number of Subgraphs', 'Std'),
        ('Modularity', 'Mean'), ('Modularity', 'Std'),
        ('Betweenness', 'Mean'), ('Betweenness', 'Std')
        
    ])

    time_step_0=flatten_comprehension([assortativity_race_global_summary[0]+assortativity_age_global_summary[0]+assortativity_degree_global_summary[0]+clustering_global_summary[0]\
                 +average_diameter_global_summary[0]+average_shortest_path_length_global_summary[0]+largest_component_summary[0]+number_of_subgraphs_global_summary[0]\
                    +modularity_global_summary[0]+betweenness_global_summary[0]])
    
    time_step_1=flatten_comprehension([assortativity_race_global_summary[1]+assortativity_age_global_summary[1]+assortativity_degree_global_summary[1]+clustering_global_summary[1]\
                 +average_diameter_global_summary[1]+average_shortest_path_length_global_summary[1]+largest_component_summary[1]+number_of_subgraphs_global_summary[1]\
                    +modularity_global_summary[1]+betweenness_global_summary[1]])
    
    time_step_2=flatten_comprehension([assortativity_race_global_summary[2]+assortativity_age_global_summary[2]+assortativity_degree_global_summary[2]+clustering_global_summary[2]\
                 +average_diameter_global_summary[2]+average_shortest_path_length_global_summary[2]+largest_component_summary[2]+number_of_subgraphs_global_summary[2]\
                    +modularity_global_summary[2]+betweenness_global_summary[2]])
    
    time_step_3=flatten_comprehension([assortativity_race_global_summary[3]+assortativity_age_global_summary[3]+assortativity_degree_global_summary[3]+clustering_global_summary[3]\
                 +average_diameter_global_summary[3]+average_shortest_path_length_global_summary[3]+largest_component_summary[3]+number_of_subgraphs_global_summary[3]\
                    +modularity_global_summary[3]+betweenness_global_summary[3]])
    
    time_step_4=flatten_comprehension([assortativity_race_global_summary[4]+assortativity_age_global_summary[4]+assortativity_degree_global_summary[4]+clustering_global_summary[4]\
                 +average_diameter_global_summary[4]+average_shortest_path_length_global_summary[4]+largest_component_summary[4]+number_of_subgraphs_global_summary[4]\
                    +modularity_global_summary[4]+betweenness_global_summary[4]])
    
    time_step_5=flatten_comprehension([assortativity_race_global_summary[5]+assortativity_age_global_summary[5]+assortativity_degree_global_summary[5]+clustering_global_summary[5]\
                 +average_diameter_global_summary[5]+average_shortest_path_length_global_summary[5]+largest_component_summary[5]+number_of_subgraphs_global_summary[5]\
                    +modularity_global_summary[5]+betweenness_global_summary[5]])

    # Sample data
    data = [time_step_0,time_step_1,time_step_2,time_step_3,time_step_4,time_step_5]

    # Create a DataFrame with multi-level columns
    dataframe_metrics_assortativity_race_effect = pd.DataFrame(data, columns=columns)

    file_name = "dataframe_metrics_assortativity_race_effect_new_with_modularity_and_betweenness.pkl"

    open_file = open(file_name, "wb")
    pickle.dump(dataframe_metrics_assortativity_race_effect, open_file)
    open_file.close()