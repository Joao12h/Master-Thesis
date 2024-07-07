from clustering_assortativity_functions import*
from auxiliary_functions_epidemic_simulation import*
import sys
import pickle
import numpy as np
import os
from pathlib import Path



property_under_study=sys.argv[1]



#based on stanford medicine healthcare

theta_from_top_to_bottom=theta_from_versatile_to_bottom=theta_from_top_to_versatile=0.014
theta_from_bottom_to_top=theta_from_bottom_to_versatile=theta_from_versatile_to_top=0.0062
theta_from_versatile_to_versatile=theta_from_top_to_bottom+theta_from_bottom_to_top
theta_from_bottom_to_bottom=theta_from_top_to_top=theta_from_top_to_no_anal_intercourse=\
    theta_from_bottom_to_no_anal_intercourse=theta_from_no_anal_intercourse_to_no_anal_intercourse=\
    theta_from_versatile_to_no_anal_intercourse=theta_from_no_anal_intercourse_to_bottom =theta_from_no_anal_intercourse_to_top=\
            theta_from_no_anal_intercourse_to_versatile=0

transmission_probability_sexual_role_matrix=\
    [[theta_from_no_anal_intercourse_to_no_anal_intercourse,theta_from_versatile_to_no_anal_intercourse,\
      theta_from_bottom_to_no_anal_intercourse,theta_from_top_to_no_anal_intercourse], \
        [theta_from_versatile_to_no_anal_intercourse,theta_from_versatile_to_versatile,theta_from_versatile_to_bottom,\
         theta_from_versatile_to_top], [theta_from_bottom_to_no_anal_intercourse,theta_from_bottom_to_versatile,\
                                        theta_from_bottom_to_bottom,theta_from_bottom_to_top],\
                                            [theta_from_top_to_no_anal_intercourse,theta_from_top_to_versatile,\
                                             theta_from_top_to_bottom,theta_from_top_to_top]]


#get baseline networks

if property_under_study=='Network Structure Effect':

  current_directory = Path(__file__).parent #Get current directory
  file = open(os.path.join(current_directory, 'baseline_networks.pkl'), 'rb') #rb = read bytes because we are reading the file
  baseline_networks = pickle.load(file)
  file.close()


  network_structure=[]
  for type_of_network in range(0,len(baseline_networks)):
      save_summary_graph=[]
      for graph in baseline_networks[type_of_network]:
          
          for node in graph.nodes():
            graph.nodes[node]['Sexual_position'] = 'versatile'

          largest_component = max(nx.connected_components(graph), key=len)
          largest_component_graph = graph.subgraph(largest_component)
          results_epidemics=synchronous_update_to_study_network_structure(graph,transmission_probability_sexual_role_matrix,time_end=50,N_trials=30)
          prevalence_all_stochastic_events=[]

          for stochastic_event in range(0,len(results_epidemics)):
            infected_history_per_stochastic_event=[current_infected/len(largest_component_graph) for current_infected in results_epidemics[stochastic_event][2]]
            #incidence_per_stochastic_event=computeincidence(infected_history_per_stochastic_event,4667)
            prevalence_all_stochastic_events.append(infected_history_per_stochastic_event)


          infected_all_stochastic_events_array= np.array(prevalence_all_stochastic_events)
          summary_statistics=np.percentile(infected_all_stochastic_events_array, [25, 50, 75], axis=0)
          save_summary_graph.append(summary_statistics)

      network_structure.append(save_summary_graph)


  # open a file, where you ant to store the data

  doc_name= 'baseline_incidence_networks_in_terms_of_prevalence_all_vers.pkl'

  #doc_name= 'clustered_networks_30000_iterations_improved_ID_12.pkl'

  file = open(doc_name, 'wb')

  # dump information to that file
  pickle.dump(network_structure, file)

  # close the file
  file.close()

# Calculate median and quartiles for each timestep


if property_under_study=='PrEP + condom + highest_degree':
        current_directory = Path(__file__).parent #Get current directory
        file = open(os.path.join(current_directory, 'networks_prevention_strategies.pkl'), 'rb') #rb = read bytes because we are reading the file
        networks_prevention_strategies = pickle.load(file)
        file.close()

        for graph in networks_prevention_strategies:
            for node in graph.nodes():
                graph.nodes[node]['Sexual_position'] = 'versatile'
                


        network_structure=[]
        for graph in networks_prevention_strategies:
          
          largest_component = max(nx.connected_components(graph), key=len)
          largest_component_graph = graph.subgraph(largest_component)
          save_summary_graph=[]
          results_epidemics_baseline=synchronous_update_best_prevention_strategies_faster(graph,transmission_probability_sexual_role_matrix,prep_strategy=0,coverage_increase=[0],time_end=50,N_trials=30)
          results_epidemics_prep=synchronous_update_best_prevention_strategies_faster(graph,transmission_probability_sexual_role_matrix,prep_strategy='highest_degree_centrality',coverage_increase=[0.05,0.10,0.20,0.40,0.60],time_end=50,N_trials=30)
          results_epidemics=results_epidemics_baseline+results_epidemics_prep
          results_per_coverage=[]
          for coverage in range(0,len(results_epidemics)):
            incidence_per_efficacy=[]
            for efficacy in range(0,len(results_epidemics[0])):
              incidence_all_stochastic_events=[]
              for stochastic_event in range(0,len(results_epidemics[0][0])):
                  infected_history_per_stochastic_event=[current_infected/len(largest_component_graph) for current_infected in results_epidemics[coverage][efficacy][stochastic_event][2]]
                  #incidence_per_stochastic_event=computeincidence(infected_history_per_stochastic_event,4667)
                  incidence_all_stochastic_events.append(infected_history_per_stochastic_event)

              infected_all_stochastic_events_array= np.array(incidence_all_stochastic_events)
              summary_statistics=np.percentile(infected_all_stochastic_events_array, [25, 50, 75], axis=0)
              incidence_per_efficacy.append(summary_statistics)
            results_per_coverage.append(incidence_per_efficacy)
          network_structure.append(results_per_coverage)

          

        # Define a list of names for the pickle files
        file_names = ["clustered_prevention_strategies_prep_all_vers_highest_degree_centrality.pkl", "assorted_by_degree_prevention_strategies_prep_all_vers_highest_degree_centrality.pkl", "assorted_by_age_prevention_strategies_prep_all_vers_highest_degree_centrality.pkl", "assorted_by_race_prevention_strategies_prep_all_vers_highest_degree_centrality.pkl"]  # Add your desired file names here

        # Make sure the number of file names matches the number of sublists in network_structure
        if len(file_names) != len(network_structure):
            raise ValueError("Number of file names should match the number of sublists in network_structure")

        for file_name, sublist in zip(file_names, network_structure):
            file_path = os.path.join(current_directory, file_name)
            
            with open(file_path, 'wb') as file_prevention_strategies:
                pickle.dump(sublist, file_prevention_strategies)

if property_under_study=='PrEP + condom + kshell':
        current_directory = Path(__file__).parent #Get current directory
        file = open(os.path.join(current_directory, 'networks_prevention_strategies.pkl'), 'rb') #rb = read bytes because we are reading the file
        networks_prevention_strategies = pickle.load(file)
        file.close()

        for graph in networks_prevention_strategies:
            for node in graph.nodes():
                graph.nodes[node]['Sexual_position'] = 'versatile'
                


        network_structure=[]
        for graph in networks_prevention_strategies:
          
          largest_component = max(nx.connected_components(graph), key=len)
          largest_component_graph = graph.subgraph(largest_component)
          save_summary_graph=[]
          results_epidemics_baseline=synchronous_update_best_prevention_strategies_faster(graph,transmission_probability_sexual_role_matrix,prep_strategy=0,coverage_increase=[0],time_end=50,N_trials=30)
          results_epidemics_prep=synchronous_update_best_prevention_strategies_faster(graph,transmission_probability_sexual_role_matrix,prep_strategy='kshell', coverage_increase=[0.05,0.10,0.20,0.40,0.60],time_end=50,N_trials=30)
          results_epidemics=results_epidemics_baseline+results_epidemics_prep
          results_per_coverage=[]
          for coverage in range(0,len(results_epidemics)):
            incidence_per_efficacy=[]
            for efficacy in range(0,len(results_epidemics[0])):
              incidence_all_stochastic_events=[]
              for stochastic_event in range(0,len(results_epidemics[0][0])):
                  infected_history_per_stochastic_event=[current_infected/len(largest_component_graph) for current_infected in results_epidemics[coverage][efficacy][stochastic_event][2]]
                  #incidence_per_stochastic_event=computeincidence(infected_history_per_stochastic_event,4667)
                  incidence_all_stochastic_events.append(infected_history_per_stochastic_event)

              infected_all_stochastic_events_array= np.array(incidence_all_stochastic_events)
              summary_statistics=np.percentile(infected_all_stochastic_events_array, [25, 50, 75], axis=0)
              incidence_per_efficacy.append(summary_statistics)
            results_per_coverage.append(incidence_per_efficacy)
          network_structure.append(results_per_coverage)

          

        # Define a list of names for the pickle files
        file_names = ["clustered_prevention_strategies_prep_all_vers_kshell.pkl", "assorted_by_degree_prevention_strategies_prep_all_vers_kshell.pkl", "assorted_by_age_prevention_strategies_prep_all_vers_kshell.pkl", "assorted_by_race_prevention_strategies_prep_all_vers_kshell.pkl"]  # Add your desired file names here

        # Make sure the number of file names matches the number of sublists in network_structure
        if len(file_names) != len(network_structure):
            raise ValueError("Number of file names should match the number of sublists in network_structure")

        for file_name, sublist in zip(file_names, network_structure):
            file_path = os.path.join(current_directory, file_name)
            
            with open(file_path, 'wb') as file_prevention_strategies:
                pickle.dump(sublist, file_prevention_strategies)


if property_under_study=='PrEP + condom + testing + random':
        

        current_directory = Path(__file__).parent #Get current directory
        file = open(os.path.join(current_directory, 'networks_prevention_strategies.pkl'), 'rb') #rb = read bytes because we are reading the file
        networks_prevention_strategies = pickle.load(file)
        file.close()

        file = open(os.path.join(current_directory, 'agents_complete_final_use_in_analysis.pkl'), 'rb') #rb = read bytes because we are reading the file
        agents_complete = pickle.load(file)
        file.close()

        for graph in networks_prevention_strategies:
            graph=graph_with_atributes(graph,agents_complete)


        network_structure=[]
        for graph in networks_prevention_strategies:
          save_summary_graph=[]

          largest_component = max(nx.connected_components(graph), key=len)
          largest_component_graph = graph.subgraph(largest_component)

          results_epidemics_baseline=synchronous_update_best_prevention_strategies_faster_with_testing(graph,transmission_probability_sexual_role_matrix,prep_strategy=0,coverage_increase=[0],time_end=50,N_trials=30)
          results_epidemics_prep=synchronous_update_best_prevention_strategies_faster_with_testing(graph,transmission_probability_sexual_role_matrix,prep_strategy='random',coverage_increase=[0.05,0.10,0.20,0.40,0.60],time_end=50,N_trials=30)
          results_epidemics=results_epidemics_baseline+results_epidemics_prep
          results_per_coverage=[]
          for coverage in range(0,len(results_epidemics)):
            incidence_per_efficacy=[]
            for efficacy in range(0,len(results_epidemics[0])):
              incidence_all_stochastic_events=[]
              for stochastic_event in range(0,len(results_epidemics[0][0])):
                  infected_history_per_stochastic_event=[current_infected/len(largest_component_graph) for current_infected in results_epidemics[coverage][efficacy][stochastic_event][2]]
                  #incidence_per_stochastic_event=computeincidence(infected_history_per_stochastic_event,4667)
                  incidence_all_stochastic_events.append(infected_history_per_stochastic_event)

              infected_all_stochastic_events_array= np.array(incidence_all_stochastic_events)
              summary_statistics=np.percentile(infected_all_stochastic_events_array, [25, 50, 75], axis=0)
              incidence_per_efficacy.append(summary_statistics)
            results_per_coverage.append(incidence_per_efficacy)
          network_structure.append(results_per_coverage)

          

        # Define a list of names for the pickle files
        file_names = ["clustered_prevention_strategies_prep_testing_all_vers_random.pkl", "assorted_by_degree_prevention_strategies_prep_testing_all_vers_random.pkl", "assorted_by_age_prevention_strategies_prep_testing_all_vers_random.pkl", "assorted_by_race_prevention_strategies_prep_testing_all_vers_random.pkl"]  # Add your desired file names here

        # Make sure the number of file names matches the number of sublists in network_structure
        if len(file_names) != len(network_structure):
            raise ValueError("Number of file names should match the number of sublists in network_structure")

        for file_name, sublist in zip(file_names, network_structure):
            file_path = os.path.join(current_directory, file_name)
            
            with open(file_path, 'wb') as file_prevention_strategies:
                pickle.dump(sublist, file_prevention_strategies)


if property_under_study=='PrEP + condom + testing + highest_degree':
        

        current_directory = Path(__file__).parent #Get current directory
        file = open(os.path.join(current_directory, 'networks_prevention_strategies.pkl'), 'rb') #rb = read bytes because we are reading the file
        networks_prevention_strategies = pickle.load(file)
        file.close()

        file = open(os.path.join(current_directory, 'agents_complete_final_use_in_analysis.pkl'), 'rb') #rb = read bytes because we are reading the file
        agents_complete = pickle.load(file)
        file.close()

        for graph in networks_prevention_strategies:
            graph=graph_with_atributes(graph,agents_complete)


        network_structure=[]
        for graph in networks_prevention_strategies:
          save_summary_graph=[]

          largest_component = max(nx.connected_components(graph), key=len)
          largest_component_graph = graph.subgraph(largest_component)

          results_epidemics_baseline=synchronous_update_best_prevention_strategies_faster_with_testing(graph,transmission_probability_sexual_role_matrix,prep_strategy=0,coverage_increase=[0],time_end=50,N_trials=30)
          results_epidemics_prep=synchronous_update_best_prevention_strategies_faster_with_testing(graph,transmission_probability_sexual_role_matrix,prep_strategy='highest_degree_centrality',coverage_increase=[0.05,0.10,0.20,0.40,0.60],time_end=50,N_trials=30)
          results_epidemics=results_epidemics_baseline+results_epidemics_prep
          results_per_coverage=[]
          for coverage in range(0,len(results_epidemics)):
            incidence_per_efficacy=[]
            for efficacy in range(0,len(results_epidemics[0])):
              incidence_all_stochastic_events=[]
              for stochastic_event in range(0,len(results_epidemics[0][0])):
                  infected_history_per_stochastic_event=[current_infected/len(largest_component_graph) for current_infected in results_epidemics[coverage][efficacy][stochastic_event][2]]
                  #incidence_per_stochastic_event=computeincidence(infected_history_per_stochastic_event,4667)
                  incidence_all_stochastic_events.append(infected_history_per_stochastic_event)

              infected_all_stochastic_events_array= np.array(incidence_all_stochastic_events)
              summary_statistics=np.percentile(infected_all_stochastic_events_array, [25, 50, 75], axis=0)
              incidence_per_efficacy.append(summary_statistics)
            results_per_coverage.append(incidence_per_efficacy)
          network_structure.append(results_per_coverage)

          

        # Define a list of names for the pickle files
        file_names = ["clustered_prevention_strategies_prep_testing_all_vers_highest_degree.pkl", "assorted_by_degree_prevention_strategies_prep_testing_all_vers_highest_degree.pkl", "assorted_by_age_prevention_strategies_prep_testing_all_vers_highest_degree.pkl", "assorted_by_race_prevention_strategies_prep_testing_all_vers_highest_degree.pkl"]  # Add your desired file names here

        # Make sure the number of file names matches the number of sublists in network_structure
        if len(file_names) != len(network_structure):
            raise ValueError("Number of file names should match the number of sublists in network_structure")

        for file_name, sublist in zip(file_names, network_structure):
            file_path = os.path.join(current_directory, file_name)
            
            with open(file_path, 'wb') as file_prevention_strategies:
                pickle.dump(sublist, file_prevention_strategies)


if property_under_study=='PrEP + condom + testing + kshell':
        

        current_directory = Path(__file__).parent #Get current directory
        file = open(os.path.join(current_directory, 'networks_prevention_strategies.pkl'), 'rb') #rb = read bytes because we are reading the file
        networks_prevention_strategies = pickle.load(file)
        file.close()

        file = open(os.path.join(current_directory, 'agents_complete_final_use_in_analysis.pkl'), 'rb') #rb = read bytes because we are reading the file
        agents_complete = pickle.load(file)
        file.close()

        for graph in networks_prevention_strategies:
            graph=graph_with_atributes(graph,agents_complete)


        network_structure=[]
        for graph in networks_prevention_strategies:
          save_summary_graph=[]

          largest_component = max(nx.connected_components(graph), key=len)
          largest_component_graph = graph.subgraph(largest_component)

          results_epidemics_baseline=synchronous_update_best_prevention_strategies_faster_with_testing(graph,transmission_probability_sexual_role_matrix,prep_strategy=0,coverage_increase=[0],time_end=50,N_trials=30)
          results_epidemics_prep=synchronous_update_best_prevention_strategies_faster_with_testing(graph,transmission_probability_sexual_role_matrix,prep_strategy='kshell',coverage_increase=[0.05,0.10,0.20,0.40,0.60],time_end=50,N_trials=30)
          results_epidemics=results_epidemics_baseline+results_epidemics_prep
          results_per_coverage=[]
          for coverage in range(0,len(results_epidemics)):
            incidence_per_efficacy=[]
            for efficacy in range(0,len(results_epidemics[0])):
              incidence_all_stochastic_events=[]
              for stochastic_event in range(0,len(results_epidemics[0][0])):
                  infected_history_per_stochastic_event=[current_infected/len(largest_component_graph) for current_infected in results_epidemics[coverage][efficacy][stochastic_event][2]]
                  #incidence_per_stochastic_event=computeincidence(infected_history_per_stochastic_event,4667)
                  incidence_all_stochastic_events.append(infected_history_per_stochastic_event)

              infected_all_stochastic_events_array= np.array(incidence_all_stochastic_events)
              summary_statistics=np.percentile(infected_all_stochastic_events_array, [25, 50, 75], axis=0)
              incidence_per_efficacy.append(summary_statistics)
            results_per_coverage.append(incidence_per_efficacy)
          network_structure.append(results_per_coverage)

          

        # Define a list of names for the pickle files
        file_names = ["clustered_prevention_strategies_prep_testing_all_vers_kshell.pkl", "assorted_by_degree_prevention_strategies_prep_testing_all_vers_kshell.pkl", "assorted_by_age_prevention_strategies_prep_testing_all_vers_kshell.pkl", "assorted_by_race_prevention_strategies_prep_testing_all_vers_kshell.pkl"]  # Add your desired file names here

        # Make sure the number of file names matches the number of sublists in network_structure
        if len(file_names) != len(network_structure):
            raise ValueError("Number of file names should match the number of sublists in network_structure")

        for file_name, sublist in zip(file_names, network_structure):
            file_path = os.path.join(current_directory, file_name)
            
            with open(file_path, 'wb') as file_prevention_strategies:
                pickle.dump(sublist, file_prevention_strategies)




if property_under_study=='Network Structure Effect Plots':

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'baseline_incidence_networks_in_terms_of_prevalence_all_vers.pkl'), 'rb') #rb = read bytes because we are reading the file
    baseline_incidence_networks = pickle.load(file)
    file.close()


    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 2, figsize=(20, 14))

    # Plot on each subplot

    #Clustering
    axs[0,0].plot(time_steps, baseline_incidence_networks[0][0][1],label=r'$C=0.0350$',color='blue')
    axs[0,0].fill_between(time_steps, baseline_incidence_networks[0][0][0], baseline_incidence_networks[0][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, baseline_incidence_networks[0][1][1],label=r'$C=0.302$',color='orange')
    axs[0,0].fill_between(time_steps, baseline_incidence_networks[0][1][0], baseline_incidence_networks[0][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, baseline_incidence_networks[0][2][1],label=r'$C=0.355$',color='red')
    axs[0,0].fill_between(time_steps, baseline_incidence_networks[0][2][0], baseline_incidence_networks[0][2][2], color='red', alpha=0.3)

    axs[0,0].plot(time_steps, baseline_incidence_networks[0][3][1],label=r'$C=0.382$',color='green')
    axs[0,0].fill_between(time_steps, baseline_incidence_networks[0][3][0], baseline_incidence_networks[0][3][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, baseline_incidence_networks[0][4][1],label=r'$C=0.401$',color='brown')
    axs[0,0].fill_between(time_steps, baseline_incidence_networks[0][4][0], baseline_incidence_networks[0][4][2], color='brown', alpha=0.3)

    axs[0,0].plot(time_steps, baseline_incidence_networks[0][5][1],label=r'$C=0.417$',color='black')
    axs[0,0].fill_between(time_steps, baseline_incidence_networks[0][5][0], baseline_incidence_networks[0][5][2], color='black', alpha=0.3)

    axs[0,0].set_title('Clustering Effect Epidemics')
    axs[0,0].set_xlabel('Timestep')
    axs[0,0].set_ylabel('HIV incidence')
    axs[0,0].legend()
    axs[0,0].grid(True)

    #Assortativity by Degree

    axs[0,1].plot(time_steps, baseline_incidence_networks[1][0][1],label=r'$A_d=-0.0540$',color='b')
    axs[0,1].fill_between(time_steps, baseline_incidence_networks[1][0][0], baseline_incidence_networks[1][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, baseline_incidence_networks[1][1][1],label=r'$A_d=0.376$',color='orange')
    axs[0,1].fill_between(time_steps, baseline_incidence_networks[1][1][0], baseline_incidence_networks[1][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, baseline_incidence_networks[1][2][1],label=r'$A_d=0.442$',color='red')
    axs[0,1].fill_between(time_steps, baseline_incidence_networks[1][2][0], baseline_incidence_networks[1][2][2], color='red', alpha=0.3)

    axs[0,1].plot(time_steps, baseline_incidence_networks[1][3][1],label=r'$A_d=0.466$',color='green')
    axs[0,1].fill_between(time_steps, baseline_incidence_networks[1][3][0], baseline_incidence_networks[1][3][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, baseline_incidence_networks[1][4][1],label=r'$A_d=0.479$',color='brown')
    axs[0,1].fill_between(time_steps, baseline_incidence_networks[1][4][0], baseline_incidence_networks[1][4][2], color='brown', alpha=0.3)

    axs[0,1].plot(time_steps, baseline_incidence_networks[1][5][1],label=r'$A_d=0.486$',color='black')
    axs[0,1].fill_between(time_steps, baseline_incidence_networks[1][5][0], baseline_incidence_networks[1][5][2], color='black', alpha=0.3)

    axs[0,1].set_title('Assortativity by Degree Epidemics')
    axs[0,1].set_xlabel('Timestep')
    axs[0,1].set_ylabel('HIV incidence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Assortativity by Age

    axs[1,0].plot(time_steps, baseline_incidence_networks[2][0][1],label=r'$A_a=-0.00672$',color='b')
    axs[1,0].fill_between(time_steps, baseline_incidence_networks[2][0][0], baseline_incidence_networks[2][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, baseline_incidence_networks[2][1][1],label=r'$A_a=0.805$',color='orange')
    axs[1,0].fill_between(time_steps, baseline_incidence_networks[2][1][0], baseline_incidence_networks[2][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, baseline_incidence_networks[2][2][1],label=r'$A_a=0.915$',color='red')
    axs[1,0].fill_between(time_steps, baseline_incidence_networks[2][2][0], baseline_incidence_networks[2][2][2], color='red', alpha=0.3)

    axs[1,0].plot(time_steps, baseline_incidence_networks[2][3][1],label=r'$A_a=0.950$',color='green')
    axs[1,0].fill_between(time_steps, baseline_incidence_networks[2][3][0], baseline_incidence_networks[2][3][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, baseline_incidence_networks[2][4][1],label=r'$A_a=0.966$',color='brown')
    axs[1,0].fill_between(time_steps, baseline_incidence_networks[2][4][0], baseline_incidence_networks[2][4][2], color='brown', alpha=0.3)

    axs[1,0].plot(time_steps, baseline_incidence_networks[2][5][1],label=r'$A_a=0.974$',color='black')
    axs[1,0].fill_between(time_steps, baseline_incidence_networks[2][5][0], baseline_incidence_networks[2][5][2], color='black', alpha=0.3)
    axs[1,0].set_title('Assortativity by Age Effect Epidemics')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,0].legend()
    axs[1,0].grid(True)


    #Assortativity by Race

    axs[1,1].plot(time_steps, baseline_incidence_networks[3][0][1],label=r'$A_r=-8.98 \times 10^{-4}$',color='b')
    axs[1,1].fill_between(time_steps, baseline_incidence_networks[3][0][0], baseline_incidence_networks[3][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, baseline_incidence_networks[3][1][1],label=r'$A_r=0.480$',color='orange')
    axs[1,1].fill_between(time_steps, baseline_incidence_networks[3][1][0], baseline_incidence_networks[3][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, baseline_incidence_networks[3][2][1],label=r'$A_r=0.632$',color='red')
    axs[1,1].fill_between(time_steps, baseline_incidence_networks[3][2][0], baseline_incidence_networks[3][2][2], color='red', alpha=0.3)

    axs[1,1].plot(time_steps, baseline_incidence_networks[3][3][1],label=r'$A_r=0.710$',color='green')
    axs[1,1].fill_between(time_steps, baseline_incidence_networks[3][3][0], baseline_incidence_networks[3][3][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, baseline_incidence_networks[3][4][1],label=r'$A_r=0.759$',color='brown')
    axs[1,1].fill_between(time_steps, baseline_incidence_networks[3][4][0], baseline_incidence_networks[3][4][2], color='brown', alpha=0.3)

    axs[1,1].plot(time_steps, baseline_incidence_networks[3][5][1],label=r'$A_r=0.793$',color='black')
    axs[1,1].fill_between(time_steps, baseline_incidence_networks[3][5][0], baseline_incidence_networks[3][5][2], color='black', alpha=0.3)

    axs[1,1].set_title('Assortativity by Race Effect Epidemics')
    #axs[3].set_xlabel('Timestep')
    #axs[3].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)


    for ax_row in axs:
        for ax in ax_row:
                ax.set_xlabel('Timestep', fontsize=30)  # Adjust fontsize as needed
                ax.set_ylabel('HIV Prevalence', fontsize=30)  # Adjust fontsize as needed
                ax.tick_params(axis='both', which='major', labelsize=30)  # Adjust fontsize as needed
                ax.tick_params(axis='both', which='minor', labelsize=30)  # Adjust fontsize as needed

    # Adjust legend and title size
    for ax_row in axs:
        for ax in ax_row:
            ax.legend(fontsize=24)  # Adjust fontsize as needed
            ax.set_title(ax.get_title(), fontsize=30)  # Adjust fontsize as needed

    # Adjust layout for better spacing
    plt.tight_layout()

    plt.savefig('network_epidemics_effect_with_quartiles_prevalence_all_vers_1.png', dpi=300)
                

    #plt.show()




if property_under_study=='Plots PrEP + condoms':
    
    #Clustering
    
    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'clustered_prevention_strategies_prep_all_vers.pkl'), 'rb') #rb = read bytes because we are reading the file
    clustered_prevention_strategies_prep = pickle.load(file)


    time_steps=list(range(0,50))

    fig, axs = plt.subplots(1, 3, figsize=(30, 10))

    # Plot on each subplot

    #Coverage=0.20

    axs[0].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0].plot(time_steps, clustered_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0].fill_between(time_steps, clustered_prevention_strategies_prep[1][0][0], clustered_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0].plot(time_steps, clustered_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0].fill_between(time_steps, clustered_prevention_strategies_prep[1][1][0], clustered_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0].plot(time_steps, clustered_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0].fill_between(time_steps, clustered_prevention_strategies_prep[1][2][0], clustered_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0].set_title('Coverage=0.20')
    #axs[0].set_xlabel('Timestep')
    #axs[0].set_ylabel('HIV Prevalence')
    axs[0].legend()
    axs[0].grid(True)

    #Coverage=0.40

    axs[1].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1].plot(time_steps, clustered_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1].fill_between(time_steps, clustered_prevention_strategies_prep[2][0][0], clustered_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[1].plot(time_steps, clustered_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1].fill_between(time_steps, clustered_prevention_strategies_prep[2][1][0], clustered_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[1].plot(time_steps, clustered_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1].fill_between(time_steps, clustered_prevention_strategies_prep[2][2][0], clustered_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[1].set_title('Coverage=0.40')
    #axs[1].set_xlabel('Timestep')
    #axs[1].set_ylabel('HIV Prevalence')
    axs[1].legend()
    axs[1].grid(True)

    ##Coverage=0.60

    axs[2].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[2].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[2].plot(time_steps, clustered_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[2].fill_between(time_steps, clustered_prevention_strategies_prep[3][0][0], clustered_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[2].plot(time_steps, clustered_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[2].fill_between(time_steps, clustered_prevention_strategies_prep[3][1][0], clustered_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[2].plot(time_steps, clustered_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[2].fill_between(time_steps, clustered_prevention_strategies_prep[3][2][0], clustered_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[2].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[2].legend()
    axs[2].grid(True)

    # Add a global title
    fig.suptitle('Clustering', fontsize='xx-large')

    # Customize individual subplot properties
    for ax in axs:
        ax.set_xlabel('Timestep', fontsize='xx-large')
        ax.set_ylabel('HIV Prevalence', fontsize='xx-large')
        ax.tick_params(axis='both', which='major', labelsize='large')
        ax.tick_params(axis='both', which='minor', labelsize='medium')
        ax.legend(fontsize='x-large')
        ax.set_title(ax.get_title(), fontsize='xx-large')


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    plt.savefig('clustered_prevention_strategies_prep_all_vers.png', dpi=300)    
    #plt.show()

    #Assortativity By Degree

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_degree_prevention_strategies_prep_all_vers.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_degree_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(1, 3, figsize=(30, 10))

    #Plot on each subplot

    #Coverage=0.20

    axs[0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][0][0], assorted_by_degree_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][1][0], assorted_by_degree_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][2][0], assorted_by_degree_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0].set_title('Coverage=0.20')
    axs[0].set_xlabel('Timestep')
    axs[0].set_ylabel('HIV Prevalence')
    axs[0].legend()
    axs[0].grid(True)

    #Coverage=0.40

    axs[1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][0][0], assorted_by_degree_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][1][0], assorted_by_degree_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][2][0], assorted_by_degree_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[1].set_title('Coverage=0.40')
    axs[1].set_xlabel('Timestep')
    axs[1].set_ylabel('HIV Prevalence')
    axs[1].legend()
    axs[1].grid(True)

    ##Coverage=0.60

    axs[2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][0][0], assorted_by_degree_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][1][0], assorted_by_degree_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][2][0], assorted_by_degree_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[2].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[2].legend()
    axs[2].grid(True)

    # Add a global title
    fig.suptitle('Assortativity By Degree', fontsize='xx-large')

    # Customize individual subplot properties
    for ax in axs:
        ax.set_xlabel('Timestep', fontsize='xx-large')
        ax.set_ylabel('HIV Prevalence', fontsize='xx-large')
        ax.tick_params(axis='both', which='major', labelsize='large')
        ax.tick_params(axis='both', which='minor', labelsize='medium')
        ax.legend(fontsize='x-large')
        ax.set_title(ax.get_title(), fontsize='xx-large')

    # Adjust layout for better spacing
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    # Adjust layout for better spacing
    plt.tight_layout()

    plt.savefig('assorted_by_degree_prevention_strategies_prep_all_vers.png', dpi=300)
                
    #plt.show()

    #Assortativity By Age

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_age_prevention_strategies_prep_all_vers.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_age_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(1, 3, figsize=(30, 10))

    # Plot on each subplot

    #Coverage=0.20

    axs[0].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][0][0], assorted_by_age_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][1][0], assorted_by_age_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][2][0], assorted_by_age_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0].set_title('Coverage=0.20')
    axs[0].set_xlabel('Timestep')
    axs[0].set_ylabel('HIV Prevalence')
    axs[0].legend()
    axs[0].grid(True)

    #Coverage=0.40

    axs[1].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][0][0], assorted_by_age_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][1][0], assorted_by_age_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][2][0], assorted_by_age_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[1].set_title('Coverage=0.40')
    axs[1].set_xlabel('Timestep')
    axs[1].set_ylabel('HIV Prevalence')
    axs[1].legend()
    axs[1].grid(True)

    ##Coverage=0.60

    axs[2].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][0][0], assorted_by_age_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][1][0], assorted_by_age_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][2][0], assorted_by_age_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[2].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[2].legend()
    axs[2].grid(True)


    # Add a global title
    fig.suptitle('Assortativity By Age', fontsize='xx-large')

    # Customize individual subplot properties
    for ax in axs:
        ax.set_xlabel('Timestep', fontsize='xx-large')
        ax.set_ylabel('HIV Prevalence', fontsize='xx-large')
        ax.tick_params(axis='both', which='major', labelsize='large')
        ax.tick_params(axis='both', which='minor', labelsize='medium')
        ax.legend(fontsize='x-large')
        ax.set_title(ax.get_title(), fontsize='xx-large')


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    # Adjust layout for better spacing
    #plt.tight_layout()

    plt.savefig('assorted_by_age_prevention_strategies_prep_all_vers.png', dpi=300)
                
    #plt.show()

    #Assortativity By Race

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_race_prevention_strategies_prep_all_vers.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_race_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(1, 3, figsize=(30, 10))

    # Plot on each subplot

    # Coverage=0.20

    axs[0].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][0][0], assorted_by_race_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][1][0], assorted_by_race_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][2][0], assorted_by_race_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0].set_title('Coverage=0.20')
    axs[0].set_xlabel('Timestep')
    axs[0].set_ylabel('HIV Prevalence')
    axs[0].legend()
    axs[0].grid(True)

    #Coverage=0.40

    axs[1].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][0][0], assorted_by_race_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][1][0], assorted_by_race_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][2][0], assorted_by_race_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[1].set_title('Coverage=0.40')
    axs[1].set_xlabel('Timestep')
    axs[1].set_ylabel('HIV Prevalence')
    axs[1].legend()
    axs[1].grid(True)

    ##Coverage=0.60

    axs[2].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][0][0], assorted_by_race_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][1][0], assorted_by_race_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][2][0], assorted_by_race_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[2].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[2].legend()
    axs[2].grid(True)

    # Add a global title
    fig.suptitle('Assortativity By Race', fontsize='xx-large')

    # Customize individual subplot properties
    for ax in axs:
        ax.set_xlabel('Timestep', fontsize='xx-large')
        ax.set_ylabel('HIV Prevalence', fontsize='xx-large')
        ax.tick_params(axis='both', which='major', labelsize='large')
        ax.tick_params(axis='both', which='minor', labelsize='medium')
        ax.legend(fontsize='x-large')
        ax.set_title(ax.get_title(), fontsize='xx-large')


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95


    plt.savefig('assorted_by_race_prevention_strategies_prep_all_vers.png', dpi=300)
                
    #plt.show()


if property_under_study=='Plots PrEP + condoms + testing':
    
    #Random

    #Clustering
    
    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'clustered_prevention_strategies_prep_testing_all_vers_random_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    clustered_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    # Coverage=0.05

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[1][0][0], clustered_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[1][1][0], clustered_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[1][2][0], clustered_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    #axs[0].set_xlabel('Timestep')
    #axs[0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)


    #Coverage=0.10

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[2][0][0], clustered_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[2][1][0], clustered_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[2][2][0], clustered_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    #axs[0].set_xlabel('Timestep')
    #axs[0].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[3][0][0], clustered_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[3][1][0], clustered_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[3][2][0], clustered_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    #axs[0].set_xlabel('Timestep')
    #axs[0].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[4][0][0], clustered_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[4][1][0], clustered_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[4][2][0], clustered_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    #axs[1].set_xlabel('Timestep')
    #axs[1].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    #Coverage=0.60

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[5][0][0], clustered_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[5][1][0], clustered_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[5][2][0], clustered_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    #Coverage=0.80

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[6][0][0], clustered_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[6][1][0], clustered_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[6][2][0], clustered_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)

    

    # Add a global title
    fig.suptitle('Clustering + Testing + Random', fontsize=28)

    # Customize individual subplot properties
    """ for ax in axs:
        ax.set_xlabel('Timestep', fontsize='xx-large')
        ax.set_ylabel('HIV Prevalence', fontsize='xx-large')
        ax.tick_params(axis='both', which='major', labelsize='large')
        ax.tick_params(axis='both', which='minor', labelsize='medium')
        ax.legend(fontsize='x-large')
        ax.set_title(ax.get_title(), fontsize='xx-large') """
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95


    plt.savefig('clustered_prevention_strategies_prep_testing_all_vers_random_new_new.png', dpi=300)    
    #plt.show()

    #Assortativity By Degree

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_degree_prevention_strategies_prep_testing_all_vers_random_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_degree_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][0][0], assorted_by_degree_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][1][0], assorted_by_degree_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][2][0], assorted_by_degree_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    axs[0,0].set_xlabel('Timestep')
    axs[0,0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)

    #Coverage=0.10

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][0][0], assorted_by_degree_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][1][0], assorted_by_degree_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][2][0], assorted_by_degree_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    axs[0,1].set_xlabel('Timestep')
    axs[0,1].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][0][0], assorted_by_degree_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][1][0], assorted_by_degree_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][2][0], assorted_by_degree_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    axs[0,2].set_xlabel('Timestep')
    axs[0,2].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[4][0][0], assorted_by_degree_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[4][1][0], assorted_by_degree_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[4][2][0], assorted_by_degree_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    axs[1,0].set_xlabel('Timestep')
    axs[1,0].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    ##Coverage=0.60

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[5][0][0], assorted_by_degree_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[5][1][0], assorted_by_degree_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[5][2][0], assorted_by_degree_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    ##Coverage=0.80

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[6][0][0], assorted_by_degree_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[6][1][0], assorted_by_degree_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[6][2][0], assorted_by_degree_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)

    # Add a global title
    fig.suptitle('Assortativity By Degree + Testing + Random', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    plt.savefig('assorted_by_degree_prevention_strategies_prep_testing_all_vers_random_new_new.png', dpi=300)
                
    #plt.show()

    #Assortativity By Age

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_age_prevention_strategies_prep_testing_all_vers_random_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_age_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][0][0], assorted_by_age_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][1][0], assorted_by_age_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][2][0], assorted_by_age_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    axs[0,0].set_xlabel('Timestep')
    axs[0,0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)


    #Coverage=0.10

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][0][0], assorted_by_age_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][1][0], assorted_by_age_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][2][0], assorted_by_age_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    axs[0,1].set_xlabel('Timestep')
    axs[0,1].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][0][0], assorted_by_age_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][1][0], assorted_by_age_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][2][0], assorted_by_age_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    axs[0,2].set_xlabel('Timestep')
    axs[0,2].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[4][0][0], assorted_by_age_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[4][1][0], assorted_by_age_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[4][2][0], assorted_by_age_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    axs[1,0].set_xlabel('Timestep')
    axs[1,0].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    ##Coverage=0.60

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[5][0][0], assorted_by_age_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[5][1][0], assorted_by_age_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[5][2][0], assorted_by_age_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    ##Coverage=0.60

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[6][0][0], assorted_by_age_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[6][1][0], assorted_by_age_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[6][2][0], assorted_by_age_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)


    # Add a global title
    fig.suptitle('Assortativity By Age + Testing + Random', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    plt.savefig('assorted_by_age_prevention_strategies_prep_testing_all_vers_random_new_new.png', dpi=300)
                
    #plt.show()

    #Assortativity By Race

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_race_prevention_strategies_prep_testing_all_vers_random_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_race_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][0][0], assorted_by_race_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][1][0], assorted_by_race_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][2][0], assorted_by_race_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    axs[0,0].set_xlabel('Timestep')
    axs[0,0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)

    #Coverage=0.10

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][0][0], assorted_by_race_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][1][0], assorted_by_race_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][2][0], assorted_by_race_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    axs[0,1].set_xlabel('Timestep')
    axs[0,1].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][0][0], assorted_by_race_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][1][0], assorted_by_race_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][2][0], assorted_by_race_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    axs[0,2].set_xlabel('Timestep')
    axs[0,2].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[4][0][0], assorted_by_race_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[4][1][0], assorted_by_race_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[4][2][0], assorted_by_race_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    axs[1,0].set_xlabel('Timestep')
    axs[1,0].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    ##Coverage=0.60

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[5][0][0], assorted_by_race_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[5][1][0], assorted_by_race_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[5][2][0], assorted_by_race_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    ##Coverage=0.80

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[6][0][0], assorted_by_race_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[6][1][0], assorted_by_race_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[6][2][0], assorted_by_race_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)


    # Add a global title
    fig.suptitle('Assortativity By Race + Testing + Random', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    plt.savefig('assorted_by_race_prevention_strategies_prep_testing_all_vers_random_new_new.png', dpi=300)
                
    #plt.show()


    #Highest Degree

    #Clustering
    
    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'clustered_prevention_strategies_prep_testing_all_vers_highest_degree_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    clustered_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[1][0][0], clustered_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[1][1][0], clustered_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[1][2][0], clustered_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    #axs[0].set_xlabel('Timestep')
    #axs[0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)


    #Coverage=0.10

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[2][0][0], clustered_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[2][1][0], clustered_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[2][2][0], clustered_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    #axs[0].set_xlabel('Timestep')
    #axs[0].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[3][0][0], clustered_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[3][1][0], clustered_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[3][2][0], clustered_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    #axs[0].set_xlabel('Timestep')
    #axs[0].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[4][0][0], clustered_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[4][1][0], clustered_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[4][2][0], clustered_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    #axs[1].set_xlabel('Timestep')
    #axs[1].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    #Coverage=0.60

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[5][0][0], clustered_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[5][1][0], clustered_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[5][2][0], clustered_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    #Coverage=0.80

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[6][0][0], clustered_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[6][1][0], clustered_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[6][2][0], clustered_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)

    

    # Add a global title
    fig.suptitle('Clustering + Testing + Highest Degree', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95


    plt.savefig('clustered_prevention_strategies_prep_testing_all_vers_highest_degree_new_new.png', dpi=300)    
    #plt.show()

    #Assortativity By Degree

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_degree_prevention_strategies_prep_testing_all_vers_highest_degree_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_degree_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][0][0], assorted_by_degree_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][1][0], assorted_by_degree_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][2][0], assorted_by_degree_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    axs[0,0].set_xlabel('Timestep')
    axs[0,0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)

    #Coverage=0.10

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][0][0], assorted_by_degree_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][1][0], assorted_by_degree_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][2][0], assorted_by_degree_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    axs[0,1].set_xlabel('Timestep')
    axs[0,1].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][0][0], assorted_by_degree_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][1][0], assorted_by_degree_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][2][0], assorted_by_degree_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    axs[0,2].set_xlabel('Timestep')
    axs[0,2].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[4][0][0], assorted_by_degree_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[4][1][0], assorted_by_degree_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[4][2][0], assorted_by_degree_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    axs[1,0].set_xlabel('Timestep')
    axs[1,0].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    ##Coverage=0.60

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[5][0][0], assorted_by_degree_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[5][1][0], assorted_by_degree_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[5][2][0], assorted_by_degree_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    ##Coverage=0.80

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[6][0][0], assorted_by_degree_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[6][1][0], assorted_by_degree_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[6][2][0], assorted_by_degree_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)

    # Add a global title
    fig.suptitle('Assortativity By Degree + Testing + Highest Degree', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    plt.savefig('assorted_by_degree_prevention_strategies_prep_testing_all_vers_highest_degree_new_new.png', dpi=300)
                
    #plt.show()

    #Assortativity By Age

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_age_prevention_strategies_prep_testing_all_vers_highest_degree_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_age_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][0][0], assorted_by_age_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][1][0], assorted_by_age_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][2][0], assorted_by_age_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    axs[0,0].set_xlabel('Timestep')
    axs[0,0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)


    #Coverage=0.10

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][0][0], assorted_by_age_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][1][0], assorted_by_age_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][2][0], assorted_by_age_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    axs[0,1].set_xlabel('Timestep')
    axs[0,1].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][0][0], assorted_by_age_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][1][0], assorted_by_age_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][2][0], assorted_by_age_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    axs[0,2].set_xlabel('Timestep')
    axs[0,2].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[4][0][0], assorted_by_age_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[4][1][0], assorted_by_age_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[4][2][0], assorted_by_age_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    axs[1,0].set_xlabel('Timestep')
    axs[1,0].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    ##Coverage=0.60

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[5][0][0], assorted_by_age_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[5][1][0], assorted_by_age_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[5][2][0], assorted_by_age_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    ##Coverage=0.80

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[6][0][0], assorted_by_age_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[6][1][0], assorted_by_age_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[6][2][0], assorted_by_age_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)


    # Add a global title
    fig.suptitle('Assortativity By Age + Testing + Highest Degree', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    plt.savefig('assorted_by_age_prevention_strategies_prep_testing_all_vers_highest_degree_new_new.png', dpi=300)
                
    #plt.show()

    #Assortativity By Race

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_race_prevention_strategies_prep_testing_all_vers_highest_degree_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_race_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][0][0], assorted_by_race_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][1][0], assorted_by_race_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][2][0], assorted_by_race_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    axs[0,0].set_xlabel('Timestep')
    axs[0,0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)

    #Coverage=0.10

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][0][0], assorted_by_race_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][1][0], assorted_by_race_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][2][0], assorted_by_race_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    axs[0,1].set_xlabel('Timestep')
    axs[0,1].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][0][0], assorted_by_race_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][1][0], assorted_by_race_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][2][0], assorted_by_race_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    axs[0,2].set_xlabel('Timestep')
    axs[0,2].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[4][0][0], assorted_by_race_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[4][1][0], assorted_by_race_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[4][2][0], assorted_by_race_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    axs[1,0].set_xlabel('Timestep')
    axs[1,0].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    ##Coverage=0.60

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[5][0][0], assorted_by_race_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[5][1][0], assorted_by_race_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[5][2][0], assorted_by_race_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    ##Coverage=0.80

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[6][0][0], assorted_by_race_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[6][1][0], assorted_by_race_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[6][2][0], assorted_by_race_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)


    # Add a global title
    fig.suptitle('Assortativity By Race + Testing + Highest Degree', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    plt.savefig('assorted_by_race_prevention_strategies_prep_testing_all_vers_highest_degree_new_new.png', dpi=300)


    #Kshell

    #Clustering
    
    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'clustered_prevention_strategies_prep_testing_all_vers_kshell_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    clustered_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[1][0][0], clustered_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[1][1][0], clustered_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, clustered_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, clustered_prevention_strategies_prep[1][2][0], clustered_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    #axs[0].set_xlabel('Timestep')
    #axs[0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)


    #Coverage=0.10

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[2][0][0], clustered_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[2][1][0], clustered_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, clustered_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, clustered_prevention_strategies_prep[2][2][0], clustered_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    #axs[0].set_xlabel('Timestep')
    #axs[0].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[3][0][0], clustered_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[3][1][0], clustered_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, clustered_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, clustered_prevention_strategies_prep[3][2][0], clustered_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    #axs[0].set_xlabel('Timestep')
    #axs[0].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[4][0][0], clustered_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[4][1][0], clustered_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, clustered_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, clustered_prevention_strategies_prep[4][2][0], clustered_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    #axs[1].set_xlabel('Timestep')
    #axs[1].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    #Coverage=0.60

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[5][0][0], clustered_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[5][1][0], clustered_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, clustered_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, clustered_prevention_strategies_prep[5][2][0], clustered_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    #Coverage=0.80

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[0][0][0], clustered_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[6][0][0], clustered_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[6][1][0], clustered_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, clustered_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, clustered_prevention_strategies_prep[6][2][0], clustered_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)

    

    # Add a global title
    fig.suptitle('Clustering + Testing + kshell', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95


    plt.savefig('clustered_prevention_strategies_prep_testing_all_vers_kshell_new_new.png', dpi=300)    
    #plt.show()

    #Assortativity By Degree

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_degree_prevention_strategies_prep_testing_all_vers_kshell_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_degree_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][0][0], assorted_by_degree_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][1][0], assorted_by_degree_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[1][2][0], assorted_by_degree_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    axs[0,0].set_xlabel('Timestep')
    axs[0,0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)

    #Coverage=0.10

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][0][0], assorted_by_degree_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][1][0], assorted_by_degree_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[2][2][0], assorted_by_degree_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    axs[0,1].set_xlabel('Timestep')
    axs[0,1].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][0][0], assorted_by_degree_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][1][0], assorted_by_degree_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[3][2][0], assorted_by_degree_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    axs[0,2].set_xlabel('Timestep')
    axs[0,2].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[4][0][0], assorted_by_degree_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[4][1][0], assorted_by_degree_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_degree_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[4][2][0], assorted_by_degree_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    axs[1,0].set_xlabel('Timestep')
    axs[1,0].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    ##Coverage=0.60

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[5][0][0], assorted_by_degree_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[5][1][0], assorted_by_degree_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_degree_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[5][2][0], assorted_by_degree_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    ##Coverage=0.80

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[0][0][0], assorted_by_degree_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[6][0][0], assorted_by_degree_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[6][1][0], assorted_by_degree_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_degree_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, assorted_by_degree_prevention_strategies_prep[6][2][0], assorted_by_degree_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)

    # Add a global title
    fig.suptitle('Assortativity By Degree + Testing + kshell', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    plt.savefig('assorted_by_degree_prevention_strategies_prep_testing_all_vers_kshell_new_new.png', dpi=300)
                
    #plt.show()

    #Assortativity By Age

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_age_prevention_strategies_prep_testing_all_vers_kshell_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_age_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][0][0], assorted_by_age_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][1][0], assorted_by_age_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[1][2][0], assorted_by_age_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    axs[0,0].set_xlabel('Timestep')
    axs[0,0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)


    #Coverage=0.10

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][0][0], assorted_by_age_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][1][0], assorted_by_age_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[2][2][0], assorted_by_age_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    axs[0,1].set_xlabel('Timestep')
    axs[0,1].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][0][0], assorted_by_age_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][1][0], assorted_by_age_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[3][2][0], assorted_by_age_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    axs[0,2].set_xlabel('Timestep')
    axs[0,2].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[4][0][0], assorted_by_age_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[4][1][0], assorted_by_age_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_age_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[4][2][0], assorted_by_age_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    axs[1,0].set_xlabel('Timestep')
    axs[1,0].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    ##Coverage=0.60

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[5][0][0], assorted_by_age_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[5][1][0], assorted_by_age_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_age_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[5][2][0], assorted_by_age_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    ##Coverage=0.80

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[0][0][0], assorted_by_age_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[6][0][0], assorted_by_age_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[6][1][0], assorted_by_age_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_age_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, assorted_by_age_prevention_strategies_prep[6][2][0], assorted_by_age_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)


    # Add a global title
    fig.suptitle('Assortativity By Age + Testing + kshell', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    plt.savefig('assorted_by_age_prevention_strategies_prep_testing_all_vers_kshell_new_new.png', dpi=300)
                
    #plt.show()

    #Assortativity By Race

    current_directory = Path(__file__).parent #Get current directory
    file = open(os.path.join(current_directory, 'assorted_by_race_prevention_strategies_prep_testing_all_vers_kshell_new_new.pkl'), 'rb') #rb = read bytes because we are reading the file
    assorted_by_race_prevention_strategies_prep = pickle.load(file)

    time_steps=list(range(0,50))

    fig, axs = plt.subplots(2, 3, figsize=(20, 10))

    # Plot on each subplot

    #Coverage=0.05

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][0][0], assorted_by_race_prevention_strategies_prep[1][0][2], color='blue', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][1][0], assorted_by_race_prevention_strategies_prep[1][1][2], color='orange', alpha=0.3)

    axs[0,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[1][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[1][2][0], assorted_by_race_prevention_strategies_prep[1][2][2], color='purple', alpha=0.3)

    axs[0,0].set_title('Coverage=0.05')
    axs[0,0].set_xlabel('Timestep')
    axs[0,0].set_ylabel('HIV Prevalence')
    axs[0,0].legend()
    axs[0,0].grid(True)

    #Coverage=0.10

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][0][0], assorted_by_race_prevention_strategies_prep[2][0][2], color='blue', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][1][0], assorted_by_race_prevention_strategies_prep[2][1][2], color='orange', alpha=0.3)

    axs[0,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[2][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[2][2][0], assorted_by_race_prevention_strategies_prep[2][2][2], color='purple', alpha=0.3)

    axs[0,1].set_title('Coverage=0.10')
    axs[0,1].set_xlabel('Timestep')
    axs[0,1].set_ylabel('HIV Prevalence')
    axs[0,1].legend()
    axs[0,1].grid(True)

    #Coverage=0.20

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][0][1],label='PrEP Efficacy=0.23',color='blue')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][0][0], assorted_by_race_prevention_strategies_prep[3][0][2], color='blue', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][1][0], assorted_by_race_prevention_strategies_prep[3][1][2], color='orange', alpha=0.3)

    axs[0,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[3][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[0,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[3][2][0], assorted_by_race_prevention_strategies_prep[3][2][2], color='purple', alpha=0.3)

    axs[0,2].set_title('Coverage=0.20')
    axs[0,2].set_xlabel('Timestep')
    axs[0,2].set_ylabel('HIV Prevalence')
    axs[0,2].legend()
    axs[0,2].grid(True)

    #Coverage=0.40

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[4][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[4][0][0], assorted_by_race_prevention_strategies_prep[4][0][2], color='blue', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[4][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[4][1][0], assorted_by_race_prevention_strategies_prep[4][1][2], color='orange', alpha=0.3)

    axs[1,0].plot(time_steps, assorted_by_race_prevention_strategies_prep[4][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,0].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[4][2][0], assorted_by_race_prevention_strategies_prep[4][2][2], color='purple', alpha=0.3)

    axs[1,0].set_title('Coverage=0.40')
    axs[1,0].set_xlabel('Timestep')
    axs[1,0].set_ylabel('HIV Prevalence')
    axs[1,0].legend()
    axs[1,0].grid(True)

    ##Coverage=0.60

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[5][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[5][0][0], assorted_by_race_prevention_strategies_prep[5][0][2], color='blue', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[5][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[5][1][0], assorted_by_race_prevention_strategies_prep[5][1][2], color='orange', alpha=0.3)

    axs[1,1].plot(time_steps, assorted_by_race_prevention_strategies_prep[5][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,1].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[5][2][0], assorted_by_race_prevention_strategies_prep[5][2][2], color='purple', alpha=0.3)

    axs[1,1].set_title('Coverage=0.60')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,1].legend()
    axs[1,1].grid(True)

    ##Coverage=0.80

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[0][0][1],label='No PrEP',color='green')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[0][0][0], assorted_by_race_prevention_strategies_prep[0][0][2], color='green', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[6][0][1],label='PrEP Efficacy=0.23',color='b')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[6][0][0], assorted_by_race_prevention_strategies_prep[6][0][2], color='blue', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[6][1][1],label='PrEP Efficacy=0.80',color='orange')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[6][1][0], assorted_by_race_prevention_strategies_prep[6][1][2], color='orange', alpha=0.3)

    axs[1,2].plot(time_steps, assorted_by_race_prevention_strategies_prep[6][2][1],label='PrEP Efficacy=0.95',color='purple')
    axs[1,2].fill_between(time_steps, assorted_by_race_prevention_strategies_prep[6][2][0], assorted_by_race_prevention_strategies_prep[6][2][2], color='purple', alpha=0.3)

    axs[1,2].set_title('Coverage=0.80')
    #axs[2].set_xlabel('Timestep')
    #axs[2].set_ylabel('HIV incidence')
    axs[1,2].legend()
    axs[1,2].grid(True)


    # Add a global title
    fig.suptitle('Assortativity By Race + Testing + kshell', fontsize=28)

    # Customize individual subplot properties
    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel('Timestep', fontsize=26)
            ax.set_ylabel('HIV Prevalence', fontsize=26)
            ax.tick_params(axis='both', which='major', labelsize=26)
            ax.tick_params(axis='both', which='minor', labelsize=26)
            ax.legend(fontsize=17,loc='upper left')
            ax.set_title(ax.get_title(), fontsize=26)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Set the top of the figure title to be at 0.95

    plt.savefig('assorted_by_race_prevention_strategies_prep_testing_all_vers_kshell_new_new.png', dpi=300)
    


    