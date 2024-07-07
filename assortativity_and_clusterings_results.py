#To import
from clustering_assortativity_functions import *
import sys

#sys arguments

#Load the data
network_ID=sys.argv[1]
uniform_sample_ID=sys.argv[2]
current_directory = Path(__file__).parent #Get current directory
file = open(os.path.join(current_directory, 'baseline.pkl'), 'rb') #rb = read bytes because we are reading the file
list_of_graphs = pickle.load(file)
file.close()
#Select and generate arbitrary network
#Select and generate arbitrary network

all_adjacency_matrices=get_adjacency_matrix(list_of_graphs)

teste=nx.Graph(all_adjacency_matrices[int(network_ID)])

# #Load data from already runned simulations

# #network_ID=sys.argv[1]
# current_directory = Path(__file__).parent #Get current directory
# file = open(os.path.join(current_directory, 'clustered_networks_10000_improved_ID_12.pkl'), 'rb') #rb = read bytes because we are reading the file
# clustered_networks_10000_improved_ID_12 = pickle.load(file)
# file.close()

# teste=clustered_networks_10000_improved_ID_12[0][-1]


# Record the start time
start_time = time.time()

#run
clustering_teste=markov_chain(teste,nb_iterations=160000,seed=None)

# Record the end time
end_time = time.time()

# Calculate the elapsed time
elapsed_time = end_time - start_time

print(f"Elapsed Time: {elapsed_time} seconds")

# open a file, where you ant to store the data

doc_name= 'clustered_networks_150000_iterations_improved_with_new_clustering_algorithm_ID_'+network_ID+'_stochastic_'+uniform_sample_ID+'.pkl'

#doc_name= 'clustered_networks_30000_iterations_improved_ID_12.pkl'

file = open(doc_name, 'wb')

# dump information to that file
pickle.dump(clustering_teste, file)

# close the file
file.close()

