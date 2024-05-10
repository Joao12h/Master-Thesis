
import os, pickle
from pathlib import Path
import networkx as nx
import community as community_louvain
from synchronous_epidemics.py import*
import matplotlib as plt

current_directory = Path(__file__).parent #Get current directory
file = open(os.path.join(current_directory, 'networks_prevention_strategies.pkl'), 'rb') 
#network networks_prevention_strategies corresponds to a list that contains the most clustered, assorted by degree, age and race networks
networks_prevention_strategies = pickle.load(file)

# largest_connected_component will allow to obtain the largest component for each network
def largest_connected_component(G):
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    return G.subgraph(components[0])


largest_component_networks_prevention_strategies=[]
for graph in networks_prevention_strategies:
    largest_component_networks_prevention_strategies.append(largest_connected_component(graph))

#get_proportion will receive a list and will compute the proportion of each element that belongs to that list. In this context, will be used to compute the proportion of nodes and infected nodes for each community detected by the louvain algorithm
def get_proportion(lista):
    proportions=[]
    for number_nodes in lista:
        proportions.append(number_nodes/sum(lista))
    return proportions


networks = largest_component_networks_prevention_strategies

all_communities=[]
all_HIV_prevalence=[]
all_degrees=[]

for i, G in enumerate(networks):
    print(f"Network {i+1}:")
    
    # Detect communities using Louvain method
    partition = community_louvain.best_partition(G)
    number_of_infected=count_initial_number(G,'HIV_status','Infected')
    
    # dictionaries are initialized to store community sizes and number of infected per community infected 
    community_sizes = {}
    infected_counts = {}
    median_degrees = {}

    # Iterate over nodes and update community sizes and infected number of agents
    for node, comm_id in partition.items():
        # Increment community size
        community_sizes[comm_id] = community_sizes.get(comm_id, 0) + 1
        
        # Check if the node is infected and increment infected count for the corresponding community
        if G.nodes[node]['HIV_status'] == 'Infected':
            infected_counts[comm_id] = infected_counts.get(comm_id, 0) + 1
        
        # Calculate the degree of the node and store it in a list for the corresponding community
        degree = G.degree[node]
        median_degrees.setdefault(comm_id, []).append(degree)

    
    for comm_id in community_sizes:
        if comm_id not in infected_counts:
            infected_counts[comm_id] = 0
        else:
            infected_counts[comm_id]=infected_counts[comm_id]
            #/number_of_infected
    
    for comm_id in community_sizes:
        community_sizes[comm_id]=community_sizes[comm_id]
        
    community_sizes=dict(sorted(community_sizes.items()))
    infected_counts=dict(sorted(infected_counts.items()))
    degrees=dict(sorted(median_degrees.items()))
    
    all_communities.append(community_sizes)
    all_HIV_prevalence.append(infected_counts)
    all_degrees.append(degrees)


#Get figures of node distribution, infected node distribution and median degree per community detected for each network type

#Clustered

# Create a new figure and divide the plotting area into a 1x3 grid of subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 5))


axs[0].bar(list(all_communities[0].keys()), get_proportion(list(all_communities[0].values())), color='blue')
axs[0].set_xlabel('Communities',fontsize=20)
axs[0].set_ylabel('Proportion of Nodes',fontsize=20)
axs[0].set_xticks(range(0, len(all_communities[0]), 2))

axs[1].bar(list(all_HIV_prevalence[0].keys()), get_proportion(list(all_HIV_prevalence[0].values())), color='orange')
axs[1].set_xlabel('Communities',fontsize=20)
axs[1].set_ylabel('Proportion of Infected Nodes',fontsize=20)
axs[1].set_xticks(range(0, len(all_HIV_prevalence[0]), 2))
#axs[1].set_title('Assortativity By Degree')

# Plot boxplot in the third subplot

positions=list(range(0, len(all_degrees[0])))

axs[2].boxplot(all_degrees[0].values(), positions=positions, patch_artist=True, boxprops=dict(facecolor='purple'), medianprops=dict(color='grey', linewidth=3), widths=0.7)
axs[2].set_xlabel('Communities',fontsize=20)
axs[2].set_ylabel('Degree Distribution',fontsize=20)
axs[2].set_yscale('log')


tick_positions = range(0, 12+1, 2)
tick_labels = range(0, 12+1, 2)
axs[2].set_xticks(tick_positions)
axs[2].set_xticklabels(tick_labels)


for ax in axs:
    ax.tick_params(axis='both', which='major', labelsize=20)


plt.suptitle('Clustering', fontsize=20)
plt.tight_layout()

plt.savefig('community_characterization_clustering.png')

plt.show()

#Assortativity By Age
# Create a new figure and divide the plotting area into a 1x3 grid of subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Plot data1 in the first subplot
axs[0].bar(list(all_communities[1].keys()), get_proportion(list(all_communities[1].values())), color='blue')
axs[0].set_xlabel('Communities',fontsize=20)
axs[0].set_ylabel('Proportion of Nodes',fontsize=20)
axs[0].set_xticks(range(0, 10+1, 1))
#axs[0].set_title('Clustering')

# Plot data2 in the second subplot
axs[1].bar(list(all_HIV_prevalence[1].keys()), get_proportion(list(all_HIV_prevalence[1].values())), color='orange')
axs[1].set_xlabel('Communities',fontsize=20)
axs[1].set_ylabel('Proportion of Infected Nodes',fontsize=20)
axs[1].set_xticks(range(0, 10+1, 1))
#axs[1].set_title('Assortativity By Degree')

data=all_degrees[1].values()

# Plot boxplot in the third subplot
positions = list(range(0, len(all_degrees[1])))
axs[2].boxplot(all_degrees[1].values(), positions=positions, patch_artist=True, boxprops=dict(facecolor='purple'), medianprops=dict(color='grey', linewidth=3), widths=0.7)
#axs[2].boxplot(data,patch_artist=True, boxprops=dict(facecolor='purple'),medianprops=dict(color='grey'),widths=0.7)
axs[2].set_xlabel('Communities',fontsize=20)
axs[2].set_ylabel('Degree Distribution',fontsize=20)
axs[2].set_yscale('log')  # Set y-axis to log scale
#axs[2].set_title('Boxplot of Degrees')

#tick_positions = range(0, 17, 2)
#tick_labels = range(0, 17, 2)
#axs[2].set_xticks(tick_positions)
#axs[2].set_xticklabels(tick_labels)

for ax in axs:
    ax.tick_params(axis='both', which='major', labelsize=20)

# Add main title to the entire figure
plt.suptitle('Assortativity By Degree',fontsize=20)
# Adjust layout to prevent overlap
plt.tight_layout()

plt.savefig('community_characterization_assortativity_by_degree.png')

# Show the plot
plt.show()

#Assortativity By Degree
# Create a new figure and divide the plotting area into a 1x3 grid of subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Plot data1 in the first subplot
axs[0].bar(list(all_communities[2].keys()), get_proportion(list(all_communities[2].values())), color='blue')
axs[0].set_xlabel('Communities',fontsize=20)
axs[0].set_ylabel('Proportion of Nodes',fontsize=20)
axs[0].set_xticks(range(0, 7+1, 1))
#axs[0].set_title('Clustering')

# Plot data2 in the second subplot
axs[1].bar(list(all_HIV_prevalence[2].keys()), get_proportion(list(all_HIV_prevalence[2].values())), color='orange')
axs[1].set_xlabel('Communities',fontsize=20)
axs[1].set_ylabel('Proportion of Infected Nodes',fontsize=20)
axs[1].set_xticks(range(0, 7+1, 1))

#axs[1].set_title('Assortativity By Degree')

# Plot boxplot in the third subplot
positions = list(range(0, len(all_degrees[2])))
axs[2].boxplot(all_degrees[2].values(), positions=positions, patch_artist=True, boxprops=dict(facecolor='purple'), medianprops=dict(color='grey', linewidth=3), widths=0.7)
#axs[2].boxplot(all_degrees[2].values(),patch_artist=True, boxprops=dict(facecolor='purple'),medianprops=dict(color='grey'),widths=0.7)
axs[2].set_xlabel('Communities',fontsize=20)
axs[2].set_ylabel('Degree Distribution',fontsize=20)
axs[2].set_yscale('log')  # Set y-axis to log scale
# Set tick positions and labels correctly
#tick_positions = range(0, 7, 1)
#tick_labels = range(0, 7, 1)
#axs[2].set_xticks(tick_positions)
#axs[2].set_xticklabels(tick_labels)

for ax in axs:
    ax.tick_params(axis='both', which='major', labelsize=20)

# Add main title to the entire figure
plt.suptitle('Assortativity By Age', fontsize=20)
# Adjust layout to prevent overlap
plt.tight_layout()

plt.savefig('community_characterization_assortativity_by_age.png')

# Show the plot
plt.show()

#Assortativity By Race

fig, axs = plt.subplots(1, 3, figsize=(15, 5))


axs[0].bar(list(all_communities[3].keys()), get_proportion(list(all_communities[3].values())), color='blue')
axs[0].set_xlabel('Communities',fontsize=20)
axs[0].set_ylabel('Proportion of Nodes',fontsize=20)
axs[0].set_xticks(range(0, 20+1, 5))


axs[1].bar(list(all_HIV_prevalence[3].keys()), get_proportion(list(all_HIV_prevalence[3].values())), color='orange')
axs[1].set_xlabel('Communities',fontsize=20)
axs[1].set_ylabel('Proportion of Infected Nodes',fontsize=20)
axs[1].set_xticks(range(0, 20+1, 5))


# Plot boxplot in the third subplot
positions = list(range(0, len(all_degrees[3])))
axs[2].boxplot(all_degrees[3].values(), positions=positions, patch_artist=True, boxprops=dict(facecolor='purple'), medianprops=dict(color='grey', linewidth=3), widths=0.7)

axs[2].set_xlabel('Communities',fontsize=20)
axs[2].set_ylabel('Degree Distribution',fontsize=20)
axs[2].set_yscale('log')  

tick_positions = range(0, 20+1, 5)
tick_labels = range(0, 20+1, 5)
axs[2].set_xticks(tick_positions)
axs[2].set_xticklabels(tick_labels)

for ax in axs:
    ax.tick_params(axis='both', which='major', labelsize=20)

plt.suptitle('Assortativity By Race', fontsize=20)

plt.tight_layout()

plt.savefig('community_characterization_assortativity_by_race.png')

# Show the plot
plt.show()


    

    

