###########################
# Run options
###########################  

levels = { 'county' }
iteration_options = { 10000 }

###########################
# Imports
###########################  

import time
import json
import os

from gerrychain import (GeographicPartition, Graph, MarkovChain, updaters, constraints, accept)
from gerrychain.tree import recursive_tree_part
from gerrychain.proposals import recom
from functools import partial

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

#import oapackage
#import numpy as np

###########################
# Hard-coded inputs
###########################  
'''
state_codes = {
    'WA': '53', 'DE': '10', 'WI': '55', 'WV': '54', 'HI': '15',
    'FL': '12', 'WY': '56', 'NJ': '34', 'NM': '35', 'TX': '48',
    'LA': '22', 'NC': '37', 'ND': '38', 'NE': '31', 'TN': '47', 'NY': '36',
    'PA': '42', 'AK': '02', 'NV': '32', 'NH': '33', 'VA': '51', 'CO': '08',
    'CA': '06', 'AL': '01', 'AR': '05', 'VT': '50', 'IL': '17', 'GA': '13',
    'IN': '18', 'IA': '19', 'MA': '25', 'AZ': '04', 'ID': '16', 'CT': '09',
    'ME': '23', 'MD': '24', 'OK': '40', 'OH': '39', 'UT': '49', 'MO': '29',
    'MN': '27', 'MI': '26', 'RI': '44', 'KS': '20', 'MT': '30', 'MS': '28',
    'SC': '45', 'KY': '21', 'OR': '41', 'SD': '46'
}
'''
state_codes = {
    'MS': '28'
}

congressional_districts = {
    'WA': 10, 'DE': 1, 'WI': 8, 'WV': 3, 'HI': 2,
    'FL': 27, 'WY': 1, 'NJ': 12, 'NM': 3, 'TX': 36,
    'LA': 6, 'NC': 13, 'ND': 1, 'NE': 3, 'TN': 9, 'NY': 27,
    'PA': 18, 'AK': 1, 'NV': 4, 'NH': 2, 'VA': 11, 'CO': 7,
    'CA': 53, 'AL': 7, 'AR': 4, 'VT': 1, 'IL': 18, 'GA': 14,
    'IN': 9, 'IA': 4, 'MA': 9, 'AZ': 9, 'ID': 2, 'CT': 5,
    'ME': 2, 'MD': 8, 'OK': 5, 'OH': 16, 'UT': 4, 'MO': 8,
    'MN': 8, 'MI': 14, 'RI': 2, 'KS': 4, 'MT': 1, 'MS': 4,
    'SC': 7, 'KY': 6, 'OR': 5, 'SD': 1
}

skips = {
    ('WA','tract'), ('WA','county'), ('DE','tract'), ('DE','county'), ('WI','tract'), 
    ('WI','county'), ('HI','tract'), ('HI','county'), ('FL','tract'), ('FL','county'), 
    ('WY','tract'), ('WY','county'), ('NJ','tract'), ('NJ','county'), ('TX','tract'), 
    ('TX','county'), ('LA','tract'), ('LA','county'), ('NC','tract'), ('NC','county'), 
    ('ND','tract'), ('ND','county'), ('TN','tract'), ('TN','county'), ('NY','tract'), 
    ('NY','county'), ('PA','tract'), ('PA','county'), ('AK','tract'), ('AK','county'), 
    ('NV','county'), ('NH','county'), ('VA','tract'), ('VA','county'), ('CO','tract'), 
    ('CO','county'), ('CA','tract'), ('CA','county'), ('AL','tract'), ('VT','tract'), 
    ('VT','county'), ('IL','tract'), ('IL','county'), ('GA','tract'), ('GA','county'),
    ('IN','tract'), ('IN','county'), ('IA','tract'), ('MA','tract'), ('MA','county'), 
    ('AZ','tract'), ('AZ','county'), ('CT','tract'), ('CT','county'), ('ME','county'), 
    ('MD','tract'), ('MD','county'), ('OK','tract'), ('OH','tract'), ('OH','county'), 
    ('UT','county'), ('MO','tract'), ('MO','county'), ('MN','tract'), ('MN','county'), 
    ('MI','tract'), ('MI','county'), ('RI','tract'), ('RI','county'), ('KS','tract'), 
    ('MT','tract'), ('MT','county'), ('SC','tract'), ('SC','county'), ('KY','tract'), 
    ('KY','county'), ('OR','tract'), ('OR','county'), ('SD','tract'), ('SD','county')
}

################################################
# Draws districts and saves to png file
################################################
def draw_pareto(datapoints, cut_minor, png_fn_pareto):
    for data in datapoints:
        front = False
        #plt.plot(data[0], data[1], 'bo')
        for pareto in cut_minor:
            if data == pareto:
                front = True
        #plt.plot(data[0], data[1], 'ro')
        if front:
            plt.plot(data[0], data[1], 'ro')
        else:
            plt.plot(data[0], data[1], 'bo')
                         
        # plt.plot(data[0], data[1], 'o')
        plt.xlabel('cut_edges')
        plt.ylabel('minority_districts')

    plt.savefig(png_fn_pareto)

def export_to_png(G, threshold, df, districts, filename1, filename2):
    
    
    #print("districts: ", districts)
    
    assignment = [ -1 for u in G.nodes ]
    
    minority_district = [0 for u in G.nodes]
    
    for j in range(len(districts)):
        total_population = 0
        minority_population = 0
        for i in districts[j]:
            geoID = G.nodes[i]["GEOID10"]
            total_population += G.node[i]["TOTPOP"]
            minority_population += G.node[i]["BVAP"]
            minority_population += G.node[i]["HVAP"]
            for u in G.nodes:
                if geoID == df['GEOID10'][u]:
                    assignment[u] = j
                    
                    
                    
        if minority_population > threshold*total_population:
            for v in districts[j]:
                geoID = G.nodes[v]["GEOID10"]
                for u in G.nodes:
                    if geoID == df['GEOID10'][u]:
                        minority_district[u] = 1
                        
    
    if min(assignment[v] for v in G.nodes) < 0:
        print("Error: did not assign all nodes in district map png.")
    else:
        df['assignment'] = assignment
        my_fig = df.plot(column='assignment').get_figure()
        RESIZE_FACTOR = 3
        my_fig.set_size_inches(my_fig.get_size_inches()*RESIZE_FACTOR)
        plt.axis('off')
        my_fig.savefig(filename1)
        
        df['minority'] = minority_district
        cmap = LinearSegmentedColormap.from_list('minority', [(0, 'white'), (1, 'lightgray')])
        splot = df.plot(cmap=cmap, column='minority',figsize=(10, 10), linewidth=1, edgecolor='0.25').get_figure()  # display the S map
        plt.axis('off')
        splot.savefig(filename2)
        
        

####################################
# Function for GerryChain call
####################################                       

def run_GerryChain_heuristic(G,population_deviation,k,threshold,iterations):
    
    my_updaters = {"population": updaters.Tally("TOTPOP", alias="population")}
    start = recursive_tree_part(G,range(k),sum(G.nodes[i]["TOTPOP"] for i in G.nodes())/k,"TOTPOP", population_deviation/2,1)
    initial_partition = GeographicPartition(G, start, updaters = my_updaters)
    
    proposal = partial(recom,
                       pop_col="TOTPOP",
                       pop_target=sum(G.nodes[i]["TOTPOP"] for i in G.nodes())/k,
                       epsilon=population_deviation/2,
                       node_repeats=2
                      )
    
    compactness_bound = constraints.UpperBound(
        lambda p: len(p["cut_edges"]),
        1.5*len(initial_partition["cut_edges"])
    )
    
    pop_constraint = constraints.within_percent_of_ideal_population(initial_partition, population_deviation/2)
    
    my_chain = MarkovChain(
        proposal=proposal,
        constraints=[
            pop_constraint,
            compactness_bound
        ],
        accept=accept.always_accept,
        initial_state=initial_partition,
        total_steps=iterations
    )
    
    min_cut_edges = sum(G[i][j]['edge_length'] for i,j in G.edges)
    print("In GerryChain heuristic, current # of cut edges and minority districts: ",end='')
    print(min_cut_edges,",",sep='',end=' ')
    max_of_minority_districts = -1
    all_maps = []
    pareto_frontier = []
    obj_vals = []
    for partition in my_chain:
        current_cut_edges = sum(G[i][j]['edge_length'] for i,j in partition["cut_edges"])
        #print(current_cut_edges,",",sep='',end=' ')
        number_minority_district = 0
        for district in range(k):
            total_pop_district = 0
            total_pop_minority = 0           
            for node in partition.graph:
                    if partition.assignment[node] == district:
                        total_pop_district += G.node[node]["TOTPOP"]
                        total_pop_minority += G.node[node]["BVAP"]
                        total_pop_minority += G.node[node]["HVAP"]
            if (total_pop_minority > threshold*total_pop_district):
                    number_minority_district += 1
        if number_minority_district > max_of_minority_districts:
            max_of_minority_districts = number_minority_district
            #best_partition = partition
        if current_cut_edges < min_cut_edges:
            #best_partition = partition
            min_cut_edges = current_cut_edges    
        print((current_cut_edges, number_minority_district),",",sep='',end=' ')
        obj_vals.append([current_cut_edges, number_minority_district])
        all_maps.append([partition, current_cut_edges, number_minority_district])
        
    
    print("Best heuristic solution has # cut edges =",min_cut_edges)
    print("Best heuristic solution has # minority districts =",max_of_minority_districts)
    #find Pareto in O(nlogn)
    all_maps.sort(key= lambda x: x[1])
    all_maps.sort(key= lambda x: x[2], reverse = True)
    pareto_frontier.append(all_maps[0])
    #highest_number_of_minority_districts = all_maps[0][1]
    least_number_of_cut_edges = all_maps[0][1]
    for i in range(1,len(all_maps)):
        if all_maps[i][1] < least_number_of_cut_edges:
            pareto_frontier.append(all_maps[i])
            least_number_of_cut_edges = all_maps[i][1]
    
    print("Pareto Frontier: ", pareto_frontier)  

    optimal_maps = []
    optimal_cut_edges = []
    optimal_minority_districts = []

    for plan in pareto_frontier:
        optimal_maps.append([[i for i in G.nodes if plan[0].assignment[i]==j] for j in range(k)])
        optimal_cut_edges.append(plan[1])
        optimal_minority_districts.append(plan[2])
    
    return (optimal_maps, optimal_cut_edges, optimal_minority_districts, obj_vals)


###########################
# Main part of the code
###########################  

# create directories for results
os.mkdir("../heuristic-results")
for iterations in iteration_options:
    os.mkdir("../heuristic-results/"+str(iterations)+"-iterations") 

# run all settings
for state in state_codes.keys():
    
    # parameters            
    k = congressional_districts[state]
    deviation = 0.01
    threshold = 0.3
    code = state_codes[state]
    
    for level in levels:
        
        # skip certain (state,level) pairs that we know:
        #   1. are infeasible, 
        #   2. are outside our scope (because of size), or
        #   3. gerrychain gets stuck on (infinite loop).
        
        #if (state,level) in skips:
         #   continue
        
        # read input graph and shapefile df
        G = Graph.from_json("../data/"+level+"/dual_graphs/"+level+code+".json")
        df = gpd.read_file("../data/"+level+"/shape_files/"+state+"_"+level+".shp")
        
        # give each edge a "length" of one
        for i,j in G.edges:
            G[i][j]['edge_length'] = 1
        
        for iterations in iteration_options:
        
            # run GerryChain 
            #print("Hello Hamid!")
            start = time.time()
            (optimal_maps, optimal_cut_edges, optimal_minority_districts, obj_vals) = run_GerryChain_heuristic(G,deviation,k,threshold,iterations)
            stop = time.time()
            
            cut_minor = []
            for i in range(len(optimal_cut_edges)):
                cut_minor.append([optimal_cut_edges[i], optimal_minority_districts[i]])
            
            
            
            #print("Optimal maps: ", optimal_maps)
            
            # filename for outputs
            fn = "../heuristic-results/"+str(iterations)+"-iterations/heur_"+state+"_"+level
            
            png_fn_pareto = fn + "_" + "pareto.png"
            
            draw_pareto(obj_vals, cut_minor, png_fn_pareto)
            
            # draw the solution on a map
            for index in range(len(optimal_maps)):
                png_fn = fn + "_" + str(index) + ".png"
                png_fn_minority = fn + "_" + str(index) + "_minority.png"
                export_to_png(G, threshold, df, optimal_maps[index], png_fn, png_fn_minority)
                
                # dump the solution info to json file
                json_fn = fn + "_" + str(index) + ".json"
                with open(json_fn, 'w') as outfile:
                    data = {}
                    data['number_cut_edges'] = optimal_cut_edges[index]
                    data['number_minority_districts'] = optimal_minority_districts[index]
                    #data['time'] = '{0:.2f}'.format(stop-start)
                    data['iterations'] = iterations
                    data['nodes'] = list()
            
                    for j in range(k):
                        for i in optimal_maps[index][j]:
                            data['nodes'].append({
                                    'name': G.nodes[i]["NAME10"],
                                    'index': i,
                                    'district': j
                                    })
                    json.dump(data, outfile)
                
