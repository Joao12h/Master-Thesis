import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt
import math
from copy import deepcopy



#get_infected_neighbours_of_each_susceptible will determine a list of sublists such that each sublist contains the infected neighbours of each susceptible node

def get_infected_neighbours_of_each_susceptible(G,largest_component_graph):
    nodes_attributes=nx.get_node_attributes(largest_component_graph,'HIV_status')
    list_infected_neighbours_per_node=[[] for i in range(0,len(G))]
    ordered_nodes=sorted(largest_component_graph.nodes)
    for node in ordered_nodes:
        for infected_node in largest_component_graph.neighbors(node):
            if nodes_attributes[node]=='Susceptible' and nodes_attributes[infected_node]=='Infected':
                list_infected_neighbours_per_node[node-1].append(infected_node)
    return list_infected_neighbours_per_node

#count_initial_number will receive a graph, atribute's node and label. It will be used in this context to determine the initial number of individuals in each infectious state

def count_initial_number(G,atribute,label):
    nodes_attributes=nx.get_node_attributes(G,atribute)
    counter_individuals_per_atribute=0

    for nodes in G.nodes:
        if nodes_attributes[nodes]==label:
            counter_individuals_per_atribute+=1
        
    return counter_individuals_per_atribute

#initial_state function receives a graph, its largest component and a given node atribute. The output will give a list containing the specified atribute of each node that belongs to the largest component  

def initial_state(graph_with_atributes_test, largest_component, atribute):
    
    initial_state_vector=[0] * len(graph_with_atributes_test)
    attributes_dict = nx.get_node_attributes(largest_component, atribute)

    # Convert the dictionary into a list of tuples
    attributes_list = list(attributes_dict.items())
    #print('attributes_dict:',attributes_list)

    # Sort the list based on the node labels
    attributes_list_sorted = sorted(attributes_list, key=lambda x: x[0])

    # Extract the attributes from the sorted list
    for node_index in range(0,len(initial_state_vector)):
        for node_state in attributes_list_sorted:
            if node_index==node_state[0]-1:
                initial_state_vector[node_index]=node_state[1]

    return initial_state_vector

#get_high_degree_centrality_nodes

def get_high_degree_centrality_nodes(graph,number_of_nodes_to_be_selected,current_state_vector):
    # Calculate degree centrality for each node
    degree_centrality = nx.degree_centrality(graph)
    # Sort nodes based on degree centrality in ascending order
    sorted_nodes = sorted(degree_centrality.items(), key=lambda x: x[1],reverse=True)
    index=0
    selected_nodes=[]
    while len(selected_nodes)!=number_of_nodes_to_be_selected:
        if current_state_vector[sorted_nodes[index][0]-1]=='Susceptible':
            selected_nodes.append(sorted_nodes[index][0])
        index+=1
    return selected_nodes

# Check if there is any node left with degree d
def check(h, d):
    f = 0  # there is no node of deg <= d
    for i in h.nodes():
        if (h.degree(i) <= d):
            f = 1
            break
    return f
 
# Find list of nodes with particular degree
def find_nodes(h, it):
    set1 = []
    for i in h.nodes():
        if (h.degree(i) <= it):
            set1.append(i)
    return set1


def get_k_shell_decomposition(graph):
    # Copy the graph
    h = graph.copy()
    it = 1
    # Bucket being filled currently
    tmp = []
    # list of lists of buckets
    buckets = []
    while (1):
        flag = check(h, it)
        if (flag == 0):
            it += 1
            buckets.append(tmp)
            tmp = []
        if (flag == 1):
            node_set = find_nodes(h, it)
            for each in node_set:
                h.remove_node(each)
                tmp.append(each)
        if (h.number_of_nodes() == 0):
            buckets.append(tmp)
            break
    return buckets

def get_k_inner_shell_nodes(bucket_list,number_of_nodes_to_be_selected,current_state_vector):
    i=len(bucket_list)-1
    inner_nodes=[]
    while i>0:
        if bucket_list[i]!=[]:
            node_in_current_k_shell=0
            while node_in_current_k_shell<len(bucket_list[i]) and len(inner_nodes)!=number_of_nodes_to_be_selected:
                if current_state_vector[bucket_list[i][node_in_current_k_shell]-1]=='Susceptible':
                    inner_nodes.append(bucket_list[i][node_in_current_k_shell])
                node_in_current_k_shell+=1
        i-=1
        #print('i:',i)
    return inner_nodes 

def increase_neighbours_rejection_sampling_with_state_vector(current_infected_neighbours_per_node_udpated, infected_node,state_vector,current_graph):
    current_infected_neighbours_per_node_udpated[infected_node-1].clear()
    neighbours_of_infected_node=list(current_graph[infected_node])
    neighbours_of_infected_node_susceptible=[selected_neighbour for selected_neighbour in neighbours_of_infected_node if(state_vector[selected_neighbour-1]=='Susceptible')]
    if neighbours_of_infected_node_susceptible!=[]:
        for neighbour in neighbours_of_infected_node_susceptible:
            neighbour=neighbour-1
            current_infected_neighbours_per_node_udpated[neighbour].append(infected_node)
    
    return current_infected_neighbours_per_node_udpated

def reduce_infected_neighbours_rejection_sampling(current_infected_neighbours,recovered_node):
    for nodes in range(0,len(current_infected_neighbours)):
        if recovered_node in current_infected_neighbours[nodes]:
            current_infected_neighbours[nodes].remove(recovered_node)
    return current_infected_neighbours


def am_I_a_bottom_or_a_top_prevention_strategies(node,state_vector_prep,condom_effectivness,state_vector_sexual_role,infected_neighbours,transmission_probability_sexual_role_matrix,total_edges_on_condoms):
    #in the article they consider several PrEP reduction risks inversely proportional to the amount of PrEP taken: reduction of 95%, 80%, 23% and 0%. Ln(1) is when there is no risk reduction
    #risk_with_prep=1-state_vector_prep[node-1]
    #print('risk_with_prep:',risk_with_prep)
    
    #convert to log odds
    #log_odds_condom=np.log(risk_with_condom)
    #log_ods_prep=np.log(risk_with_prep)
    probability_of_remaining_susceptible=1
    #neighbours_of_target_node=list(current_graph[node])
            
    reordered_list = [(x[1], x[0]) if x[0] != node else x for x in total_edges_on_condoms if x[0] == node or x[1] == node]
    
    if reordered_list!=[]:
        risk_with_condom=1-condom_effectivness
    else:
        risk_with_condom=1
    
    risk_with_prep=1-state_vector_prep[node-1]

    infected_nodes_with_which_condoms_are_used=[infected_neighbour[1] for infected_neighbour in reordered_list if (infected_neighbour[1] in infected_neighbours)]


    infected_nodes_with_which_condoms_are_not_used=[not_on_condom_infected_neighbour for not_on_condom_infected_neighbour in infected_neighbours if not_on_condom_infected_neighbour not in infected_nodes_with_which_condoms_are_used]

    

    if state_vector_sexual_role[node-1]=='receptive':
    
        transmission_probability_without_condoms=0
        transmission_probability_with_condoms=0
        count_number_of_times_I_am_bottom_being_a_bottom_not_using_condoms=len([sexual_partner for sexual_partner in infected_nodes_with_which_condoms_are_not_used if (state_vector_sexual_role[sexual_partner-1]=='versatile' or state_vector_sexual_role[sexual_partner-1]=='insertive')])
        count_number_of_times_I_am_bottom_being_a_bottom_using_condoms=len([sexual_partner for sexual_partner in infected_nodes_with_which_condoms_are_used if (state_vector_sexual_role[sexual_partner-1]=='versatile' or state_vector_sexual_role[sexual_partner-1]=='insertive')])
        baseline_probability=transmission_probability_sexual_role_matrix[1][2]
        if count_number_of_times_I_am_bottom_being_a_bottom_not_using_condoms!=0:
            log_odds_bot=np.log(baseline_probability/(1-baseline_probability))+np.log(risk_with_prep)
            transmission_probability_without_condoms=math.exp(log_odds_bot)/(1+math.exp(log_odds_bot))
            #print('probability_of_remaining_susceptible if node is a bottom:',probability_of_remaining_susceptible)
        
        elif count_number_of_times_I_am_bottom_being_a_bottom_using_condoms!=0:
            baseline_probability=transmission_probability_sexual_role_matrix[1][2]
            log_odds_bot=np.log(baseline_probability/(1-baseline_probability))+np.log(risk_with_condom)+np.log(risk_with_prep)
            transmission_probability_with_condoms=math.exp(log_odds_bot)/(1+math.exp(log_odds_bot))

        probability_of_remaining_susceptible=(1-transmission_probability_with_condoms)**count_number_of_times_I_am_bottom_being_a_bottom_using_condoms*(1-transmission_probability_without_condoms)**count_number_of_times_I_am_bottom_being_a_bottom_not_using_condoms

    elif state_vector_sexual_role[node-1]=='insertive':
        #print('node is insertive')

        transmission_probability_without_condoms=0
        transmission_probability_with_condoms=0
        count_number_of_times_I_am_top_being_a_top_not_using_condoms=len([sexual_partner for sexual_partner in infected_nodes_with_which_condoms_are_not_used if (state_vector_sexual_role[sexual_partner-1]=='versatile' or state_vector_sexual_role[sexual_partner-1]=='receptive')])
        count_number_of_times_I_am_top_being_a_top_using_condoms=len([sexual_partner for sexual_partner in infected_nodes_with_which_condoms_are_used if (state_vector_sexual_role[sexual_partner-1]=='versatile' or state_vector_sexual_role[sexual_partner-1]=='receptive')])
        baseline_probability=transmission_probability_sexual_role_matrix[1][3]
        if count_number_of_times_I_am_top_being_a_top_not_using_condoms!=0:
            log_odds_top=np.log(baseline_probability/(1-baseline_probability))+np.log(risk_with_prep)
            transmission_probability_without_condoms=math.exp(log_odds_top)/(1+math.exp(log_odds_top))
            
            #print('probability_of_remaining_susceptible if node is a top:',probability_of_remaining_susceptible)
        elif count_number_of_times_I_am_top_being_a_top_using_condoms!=0:
            log_odds_top=np.log(baseline_probability/(1-baseline_probability))+np.log(risk_with_condom)+np.log(risk_with_prep)
            transmission_probability_with_condoms=math.exp(log_odds_top)/(1+math.exp(log_odds_top))
        
        probability_of_remaining_susceptible=(1-transmission_probability_with_condoms)**count_number_of_times_I_am_top_being_a_top_using_condoms*(1-transmission_probability_without_condoms)**count_number_of_times_I_am_top_being_a_top_not_using_condoms

    elif state_vector_sexual_role[node-1]=='versatile':
        
        #print('node is versatile')
        transmission_probability_bottom_without_condoms=0
        transmission_probability_top_without_condoms=0
        transmission_probability_vers_without_condoms=0
        transmission_probability_bottom_with_condoms=0
        transmission_probability_top_with_condoms=0
        transmission_probability_vers_with_condoms=0

        baseline_probability_bottom=transmission_probability_sexual_role_matrix[3][1]
        baseline_probability_top=transmission_probability_sexual_role_matrix[2][1]
        baseline_probability_vers=transmission_probability_sexual_role_matrix[1][1]

        
        count_number_of_times_I_am_bottom_being_vers_without_condoms=len([sexual_partner for sexual_partner in infected_nodes_with_which_condoms_are_not_used if (state_vector_sexual_role[sexual_partner-1]=='insertive')])
        
        count_number_of_times_I_am_top_being_vers_without_condoms=len([sexual_partner for sexual_partner in infected_nodes_with_which_condoms_are_not_used if (state_vector_sexual_role[sexual_partner-1]=='receptive')])
        
        count_number_of_times_I_am_vers_being_vers_without_condoms=len([sexual_partner for sexual_partner in infected_nodes_with_which_condoms_are_not_used if (state_vector_sexual_role[sexual_partner-1]=='versatile')])
        

        count_number_of_times_I_am_bottom_being_vers_with_condoms=len([sexual_partner for sexual_partner in infected_nodes_with_which_condoms_are_used if (state_vector_sexual_role[sexual_partner-1]=='insertive')])
        
        count_number_of_times_I_am_top_being_vers_with_condoms=len([sexual_partner for sexual_partner in infected_nodes_with_which_condoms_are_used if (state_vector_sexual_role[sexual_partner-1]=='receptive')])
        
        count_number_of_times_I_am_vers_being_vers_with_condoms=len([sexual_partner for sexual_partner in infected_nodes_with_which_condoms_are_used if (state_vector_sexual_role[sexual_partner-1]=='versatile')])
        

        if count_number_of_times_I_am_bottom_being_vers_without_condoms!=0:
            log_odds_vers_bottom=np.log(baseline_probability_bottom/(1-baseline_probability_bottom))+np.log(risk_with_prep)
            transmission_probability_bottom_without_condoms=math.exp(log_odds_vers_bottom)/(1+math.exp(log_odds_vers_bottom))

        if count_number_of_times_I_am_top_being_vers_without_condoms!=0:
            log_odds_vers_top=np.log(baseline_probability_top/(1-baseline_probability_top))+np.log(risk_with_prep)
            transmission_probability_top_without_condoms=math.exp(log_odds_vers_top)/(1+math.exp(log_odds_vers_top))
        
        if count_number_of_times_I_am_vers_being_vers_without_condoms!=0:
            log_odds_vers_vers=np.log(baseline_probability_vers/(1-baseline_probability_vers))+np.log(risk_with_prep)
            transmission_probability_vers_without_condoms=math.exp(log_odds_vers_vers)/(1+math.exp(log_odds_vers_vers))
        
        if count_number_of_times_I_am_bottom_being_vers_with_condoms!=0:
            log_odds_vers_bottom=np.log(baseline_probability_bottom/(1-baseline_probability_bottom))+np.log(risk_with_condom)+np.log(risk_with_prep)
            transmission_probability_bottom_with_condoms=math.exp(log_odds_vers_bottom)/(1+math.exp(log_odds_vers_bottom))

        if count_number_of_times_I_am_top_being_vers_with_condoms!=0:
            log_odds_vers_top=np.log(baseline_probability_top/(1-baseline_probability_top))+np.log(risk_with_condom)+np.log(risk_with_prep)
            transmission_probability_top_with_condoms=math.exp(log_odds_vers_top)/(1+math.exp(log_odds_vers_top))
        
        if count_number_of_times_I_am_vers_being_vers_with_condoms!=0:
            log_odds_vers_vers=np.log(baseline_probability_vers/(1-baseline_probability_vers))+np.log(risk_with_condom)+np.log(risk_with_prep)
            transmission_probability_vers_with_condoms=math.exp(log_odds_vers_vers)/(1+math.exp(log_odds_vers_vers))

        probability_of_remaining_susceptible=(1-transmission_probability_bottom_without_condoms)**count_number_of_times_I_am_bottom_being_vers_without_condoms\
            *(1-transmission_probability_top_without_condoms)**count_number_of_times_I_am_top_being_vers_without_condoms*(1-transmission_probability_vers_without_condoms)**count_number_of_times_I_am_vers_being_vers_without_condoms\
            *(1-transmission_probability_bottom_with_condoms)**count_number_of_times_I_am_bottom_being_vers_with_condoms*(1-transmission_probability_top_with_condoms)**count_number_of_times_I_am_top_being_vers_with_condoms\
            *(1-transmission_probability_vers_with_condoms)**count_number_of_times_I_am_vers_being_vers_with_condoms
         
    return probability_of_remaining_susceptible

def synchronous_update_best_prevention_strategies_faster(graph_with_atributes_test,transmission_probability_sexual_role_matrix, prep_strategy,coverage_increase,time_end,N_trials):

    current_graph=graph_with_atributes_test.copy()
     # This parameter is a reference to a 2D vector where the simulation results will be stored. Each row corresponds to a simulation trial, and each column corresponds to the number of infected individuals at different intervals.

    largest_component = max(nx.connected_components(current_graph), key=len)
    largest_component_graph = current_graph.subgraph(largest_component)

    #for nodes in current_graph.nodes():
        #if nodes not in largest_component_graph.nodes():
            #current_graph.nodes[nodes]['HIV_status'] = 0

    delta_t=1

    number_of_intervals=int(time_end/delta_t)

    #prob_beta=beta*delta_t
    #prob_gama=gama*delta_t
    
    
    prevention_strategy_efficacy=[0.23, 0.80, 0.95]
    
    
    condom_effectivness=0.70

    
    #current_number_of_infected=initial_number_infected
    #save_number_infected_per_time_step=[]
    #initial_number_of_infected_neighbours_per_node=get_number_infected_neighbours(graph_with_atributes_test)[0]
    initial_infected_neighbours_per_node=get_infected_neighbours_of_each_susceptible(current_graph,largest_component_graph)
    #print('initial_infected_neighbours_per_node:',initial_infected_neighbours_per_node)
    #list_of_neighbours=[list(graph_with_atributes_test[node]) for node in sorted(graph_with_atributes_test.nodes)]
    all_nodes=sorted(largest_component_graph.nodes())

    #choose here prevention strategy to be analyzed
    

    sublist_coverage_level=[]
    for coverage_level in coverage_increase:
        sublist_efficacy_level=[]
        for  efficacy in prevention_strategy_efficacy:
            #print('initial_number_infected:',initial_number_infected)
            iteration_trial=0
            list_infected_all_trials=[]
            while iteration_trial<=N_trials:
                #print('iteration_trial:',iteration_trial)
                total_edges_on_condoms=[]
                times=[0]
                last_time=delta_t
                #Initial number per state
                current_number_of_susceptible=count_initial_number(largest_component_graph,'HIV_status','Susceptible')
                current_number_of_infected=count_initial_number(largest_component_graph,'HIV_status','Infected')
                current_number_of_treatment=count_initial_number(largest_component_graph,'HIV_status','Treatment')
                current_state_vector=initial_state(current_graph,largest_component_graph,'HIV_status')

                #history of number of people per state
                save_number_susceptible_per_time_step=[current_number_of_susceptible]
                save_number_infected_per_time_step=[current_number_of_infected]
                save_number_recovered_per_time_step=[current_number_of_treatment]
                current_infected_neighbours_per_node=initial_infected_neighbours_per_node

                state_vector_testing=initial_state(current_graph,largest_component_graph,'Testing Frequency')
                state_vector_prep=initial_state(current_graph,largest_component_graph,'PrEP adherence')
                #print('state_vector_prep_before_coverage:',state_vector_prep)
                #state_vector_condoms=initial_state(current_graph,largest_component_graph,'Condom usage')
                state_vector_sexual_position=initial_state(current_graph,largest_component_graph,'Sexual_position')

                
                #eligible nodes are the nodes that practice anal sex and are HIV negative   
                no_anal_intercourse = {node for node, attributes in largest_component_graph.nodes(data=True) if attributes.get('Sexual_position') == 'No Anal Intercourse'}
                eligible_nodes=[node for node in all_nodes if current_state_vector[node-1]=='Susceptible' and node not in no_anal_intercourse]
                number_of_nodes_to_be_selected=int(coverage_level * len(eligible_nodes))

                if prep_strategy=='random':
                    #select nodes to take prep
                    selected_nodes_prep = random.sample(eligible_nodes, k=number_of_nodes_to_be_selected)
                    for node in selected_nodes_prep:
                        state_vector_prep[node-1] = efficacy
                
                elif prep_strategy=='highest_degree_centrality':
                    selected_nodes_prep=get_high_degree_centrality_nodes(largest_component_graph,number_of_nodes_to_be_selected,current_state_vector)
                    for node in selected_nodes_prep:
                        state_vector_prep[node-1] = efficacy
                
                elif prep_strategy=='kshell':

                    bucket_list=get_k_shell_decomposition(largest_component_graph)
                    selected_nodes_prep=get_k_inner_shell_nodes(bucket_list,number_of_nodes_to_be_selected,current_state_vector)
                    for node in selected_nodes_prep:
                        state_vector_prep[node-1] = efficacy

                # Find edges in which there is anal intercourse
                edges_anal_intercourse = [(u, v) for u, v in largest_component_graph.edges() if \
                                     (u not in no_anal_intercourse) if (v not in no_anal_intercourse) if (largest_component_graph.nodes[u]['Sexual_position']!=largest_component_graph.nodes[v]['Sexual_position'] or \
                                                                        largest_component_graph.nodes[u]['Sexual_position']=='versatile' and largest_component_graph.nodes[v]['Sexual_position']=='versatile')]
                
                #print("edge_list:",edge_list)
                total_edges_on_condoms=random.sample(list(edges_anal_intercourse), k=int(0.30 * len(list(edges_anal_intercourse))))

                while last_time<time_end:

                    times.append(last_time)
                    #print('times:',times)
                    updated_state_vector=current_state_vector.copy()
                    #updated_sexual_position=state_vector_sexual_position.copy()
                    current_infected_neighbours_per_node_updated=deepcopy(current_infected_neighbours_per_node)
                    for node in all_nodes:
                        #print('node:',node)
                        #print('current_state_vector:',current_state_vector)
                        #print('current_state_vector[node-1]:',current_state_vector[node-1])
                        #print('random_probability:',random_probability)

                        if current_state_vector[node-1]=='Infected':
                            random_probability=random.random()
                            #the probability an individual will get tested is defined as the frequency of times each one gets tesed per year divided by the number of opportunities it has to get tested in one year
                            """ if random_probability<state_vector_testing[node-1]/365:
                                updated_state_vector[node-1]='Treatment'
                                current_infected_neighbours_per_node_updated=reduce_infected_neighbours_rejection_sampling(current_infected_neighbours_per_node_updated,node)
                                current_number_of_treatment+=1
                                current_number_of_infected-=1 """
                                
                        elif current_state_vector[node-1]=='Susceptible' and current_infected_neighbours_per_node[node-1]!=[]:
                            #based on article Agent-based model projections for reducing HIV infection among MSM: Prevention and care pathways to end the HIV epidemic in Chicago, Illinois
                            #print('node:',node)
                            #base risk will be the probability associated with sexual role and viral load
                                random_probability=random.random()
                                #print('random_probability:',random_probability)
                                probability_of_remaining_susceptible=am_I_a_bottom_or_a_top_prevention_strategies(node,state_vector_prep,condom_effectivness,state_vector_sexual_position,current_infected_neighbours_per_node[node-1],transmission_probability_sexual_role_matrix,total_edges_on_condoms)
                                #print('probability_of_remaining_susceptible:',probability_of_remaining_susceptible)
                                #print('probability_of_remaining_susceptible:',probability_of_remaining_susceptible)
                                if random_probability<1-probability_of_remaining_susceptible:
                                    #print('Node will get infected:',node)
                                    #print('current_infected_neighbours_per_node_before_adding_infected_neighbour:',current_infected_neighbours_per_node)
                                    current_infected_neighbours_per_node_updated=increase_neighbours_rejection_sampling_with_state_vector(current_infected_neighbours_per_node_updated,node,updated_state_vector,current_graph)
                                    updated_state_vector[node-1]='Infected'
                                    #print('updated_state_vector:',updated_state_vector)
                                    #print('current_infected_neighbours_per_node_after_adding_infected_neighbour:',current_infected_neighbours_per_node_updated)
                                    current_number_of_infected+=1
                                    current_number_of_susceptible-=1

                        #print('current_state_vector:',current_state_vector)
                        #print('updated_state_vector:',updated_state_vector)
                        #print('current_infected_neighbours_per_node_updated:',current_infected_neighbours_per_node_updated)
                        #print('current_infected_neighbours_per_node:',current_infected_neighbours_per_node)
                        #print('current_infected_neighbours_per_node_updated:',current_infected_neighbours_per_node_updated)

                    current_state_vector=updated_state_vector
                    #print('current_state_vector:',current_state_vector)
                    current_infected_neighbours_per_node=current_infected_neighbours_per_node_updated

                    #print('counter_recovered:',counter_recovered)
                    #print('counter_infected:',counter_infected)
                    #print('current_number_of_infected:',current_number_of_infected)
                    
                    save_number_susceptible_per_time_step.append(current_number_of_susceptible)
                    save_number_infected_per_time_step.append(current_number_of_infected)
                    save_number_recovered_per_time_step.append(current_number_of_treatment)
                    last_time=times[-1]+delta_t
                    print('last time:',last_time)
                
                information_iteration=[times,save_number_susceptible_per_time_step,save_number_infected_per_time_step,save_number_recovered_per_time_step]
                list_infected_all_trials.append(information_iteration)
                iteration_trial+=1
            sublist_efficacy_level.append(list_infected_all_trials)
        sublist_coverage_level.append(sublist_efficacy_level)


def synchronous_update_best_prevention_strategies_faster_with_testing(graph_with_atributes_test,transmission_probability_sexual_role_matrix, prep_strategy, coverage_increase,time_end,N_trials):

    current_graph=graph_with_atributes_test.copy()
     # This parameter is a reference to a 2D vector where the simulation results will be stored. Each row corresponds to a simulation trial, and each column corresponds to the number of infected individuals at different intervals.

    largest_component = max(nx.connected_components(current_graph), key=len)
    largest_component_graph = current_graph.subgraph(largest_component)

    #for nodes in current_graph.nodes():
        #if nodes not in largest_component_graph.nodes():
            #current_graph.nodes[nodes]['HIV_status'] = 0

    delta_t=1

    number_of_intervals=int(time_end/delta_t)

    #prob_beta=beta*delta_t
    #prob_gama=gama*delta_t
    
    
    prevention_strategy_efficacy=[0.23, 0.80, 0.95]
    
    
    condom_effectivness=0.80

    
    #current_number_of_infected=initial_number_infected
    #save_number_infected_per_time_step=[]
    #initial_number_of_infected_neighbours_per_node=get_number_infected_neighbours(graph_with_atributes_test)[0]
    initial_infected_neighbours_per_node=get_infected_neighbours_of_each_susceptible(current_graph,largest_component_graph)
    #print('initial_infected_neighbours_per_node:',initial_infected_neighbours_per_node)
    #list_of_neighbours=[list(graph_with_atributes_test[node]) for node in sorted(graph_with_atributes_test.nodes)]
    all_nodes=sorted(largest_component_graph.nodes())

    #choose here prevention strategy to be analyzed
    

    sublist_coverage_level=[]
    for coverage_level in coverage_increase:
        sublist_efficacy_level=[]
        for  efficacy in prevention_strategy_efficacy:
            #print('initial_number_infected:',initial_number_infected)
            iteration_trial=0
            list_infected_all_trials=[]
            while iteration_trial<=N_trials:
                #print('iteration_trial:',iteration_trial)
                total_edges_on_condoms=[]
                times=[0]
                last_time=delta_t
                #Initial number per state
                current_number_of_susceptible=count_initial_number(largest_component_graph,'HIV_status','Susceptible')
                current_number_of_infected=count_initial_number(largest_component_graph,'HIV_status','Infected')
                current_number_of_treatment=count_initial_number(largest_component_graph,'HIV_status','Treatment')
                current_state_vector=initial_state(current_graph,largest_component_graph,'HIV_status')

                #history of number of people per state
                save_number_susceptible_per_time_step=[current_number_of_susceptible]
                save_number_infected_per_time_step=[current_number_of_infected]
                save_number_recovered_per_time_step=[current_number_of_treatment]
                current_infected_neighbours_per_node=initial_infected_neighbours_per_node

                state_vector_testing=initial_state(current_graph,largest_component_graph,'Testing Frequency')
                state_vector_prep=initial_state(current_graph,largest_component_graph,'PrEP adherence')
                #print('state_vector_prep_before_coverage:',state_vector_prep)
                #state_vector_condoms=initial_state(current_graph,largest_component_graph,'Condom usage')
                state_vector_sexual_position=initial_state(current_graph,largest_component_graph,'Sexual_position')

                
                

                #susceptible_nodes=nodes_with_certain_atribute(current_state_vector,state_vector_prep, ['Susceptible'])

                #eligible nodes are the nodes that practice anal sex and are HIV negative   
                #no_anal_intercourse = {node for node, attributes in largest_component_graph.nodes(data=True) if attributes.get('Sexual_position') == 'No Anal Intercourse'}
                eligible_nodes=[node for node in all_nodes if current_state_vector[node-1]=='Susceptible']
                number_of_nodes_to_be_selected=int(coverage_level * len(eligible_nodes))

                #define type of strategy to distribute PrEP efficiently

                if prep_strategy=='random':
                    #select nodes to take prep
                    selected_nodes_prep = random.sample(eligible_nodes, k=number_of_nodes_to_be_selected)
                    for node in selected_nodes_prep:
                        state_vector_prep[node-1] = efficacy
                
                elif prep_strategy=='highest_degree_centrality':
                    selected_nodes_prep=get_high_degree_centrality_nodes(largest_component_graph,number_of_nodes_to_be_selected,current_state_vector)
                    print('selected_nodes:',selected_nodes_prep)
                    for node in selected_nodes_prep:
                        state_vector_prep[node-1] = efficacy
                
                elif prep_strategy=='kshell':

                    bucket_list=get_k_shell_decomposition(largest_component_graph)
                    selected_nodes_prep=get_k_inner_shell_nodes(bucket_list,number_of_nodes_to_be_selected,current_state_vector)
                    for node in selected_nodes_prep:
                        state_vector_prep[node-1] = efficacy




                # Find edges in which there is anal intercourse
                #edges_anal_intercourse = [(u, v) for u, v in largest_component_graph.edges() if \
                                     #(u not in no_anal_intercourse) if (v not in no_anal_intercourse) if (largest_component_graph.nodes[u]['Sexual_position']!=largest_component_graph.nodes[v]['Sexual_position'] or \
                                                                        #largest_component_graph.nodes[u]['Sexual_position']=='versatile' and largest_component_graph.nodes[v]['Sexual_position']=='versatile')]
                
                edges_anal_intercourse=[(u, v) for u, v in largest_component_graph.edges()]
                
                #print("edge_list:",edge_list)
                total_edges_on_condoms=random.sample(list(edges_anal_intercourse), k=int(0.30 * len(list(edges_anal_intercourse))))

                while last_time<time_end:

                    times.append(last_time)
                    #print('times:',times)
                    updated_state_vector=current_state_vector.copy()
                    #updated_sexual_position=state_vector_sexual_position.copy()
                    current_infected_neighbours_per_node_updated=deepcopy(current_infected_neighbours_per_node)
                    for node in all_nodes:
                        #print('node:',node)
                        #print('current_state_vector:',current_state_vector)
                        #print('current_state_vector[node-1]:',current_state_vector[node-1])
                        #print('random_probability:',random_probability)

                        if current_state_vector[node-1]=='Infected':
                            random_probability=random.random()
                            #the probability an individual will get tested is defined as the frequency of times each one gets tesed per year divided by the number of opportunities it has to get tested in one year
                            if random_probability<state_vector_testing[node-1]/365:
                                updated_state_vector[node-1]='Treatment'
                                current_infected_neighbours_per_node_updated=reduce_infected_neighbours_rejection_sampling(current_infected_neighbours_per_node_updated,node)
                                current_number_of_treatment+=1
                                current_number_of_infected-=1
                                
                        elif current_state_vector[node-1]=='Susceptible' and current_infected_neighbours_per_node[node-1]!=[]:
                            #based on article Agent-based model projections for reducing HIV infection among MSM: Prevention and care pathways to end the HIV epidemic in Chicago, Illinois
                            #print('node:',node)
                            #base risk will be the probability associated with sexual role and viral load
                                random_probability=random.random()
                                #print('random_probability:',random_probability)
                                probability_of_remaining_susceptible=am_I_a_bottom_or_a_top_prevention_strategies(node,state_vector_prep,condom_effectivness,state_vector_sexual_position,current_infected_neighbours_per_node[node-1],transmission_probability_sexual_role_matrix,total_edges_on_condoms)
                                #print('probability_of_remaining_susceptible:',probability_of_remaining_susceptible)
                                #print('probability_of_remaining_susceptible:',probability_of_remaining_susceptible)
                                if random_probability<1-probability_of_remaining_susceptible:
                                    #print('Node will get infected:',node)
                                    #print('current_infected_neighbours_per_node_before_adding_infected_neighbour:',current_infected_neighbours_per_node)
                                    current_infected_neighbours_per_node_updated=increase_neighbours_rejection_sampling_with_state_vector(current_infected_neighbours_per_node_updated,node,updated_state_vector,current_graph)
                                    updated_state_vector[node-1]='Infected'
                                    #print('updated_state_vector:',updated_state_vector)
                                    #print('current_infected_neighbours_per_node_after_adding_infected_neighbour:',current_infected_neighbours_per_node_updated)
                                    current_number_of_infected+=1
                                    current_number_of_susceptible-=1

                        #print('current_state_vector:',current_state_vector)
                        #print('updated_state_vector:',updated_state_vector)
                        #print('current_infected_neighbours_per_node_updated:',current_infected_neighbours_per_node_updated)
                        #print('current_infected_neighbours_per_node:',current_infected_neighbours_per_node)
                        #print('current_infected_neighbours_per_node_updated:',current_infected_neighbours_per_node_updated)

                    current_state_vector=updated_state_vector
                    #print('current_state_vector:',current_state_vector)
                    current_infected_neighbours_per_node=current_infected_neighbours_per_node_updated

                    #print('counter_recovered:',counter_recovered)
                    #print('counter_infected:',counter_infected)
                    #print('current_number_of_infected:',current_number_of_infected)
                    
                    save_number_susceptible_per_time_step.append(current_number_of_susceptible)
                    save_number_infected_per_time_step.append(current_number_of_infected)
                    save_number_recovered_per_time_step.append(current_number_of_treatment)
                    last_time=times[-1]+delta_t
                    print('last time:',last_time)
                
                information_iteration=[times,save_number_susceptible_per_time_step,save_number_infected_per_time_step,save_number_recovered_per_time_step]
                list_infected_all_trials.append(information_iteration)
                iteration_trial+=1
            sublist_efficacy_level.append(list_infected_all_trials)
        sublist_coverage_level.append(sublist_efficacy_level)


    return sublist_coverage_level
