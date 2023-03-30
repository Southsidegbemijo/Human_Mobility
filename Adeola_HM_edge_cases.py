import pycxsimulator
from pylab import *


import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from numpy import random
import kaleido as kld
import matplotlib.animation as animation
from scipy.spatial import distance
import numpy as np
import pathpy as pp
from kaleido.scopes.plotly import PlotlyScope


susceptible = 0
infected = 0
recovered = 0

populationSize = 200
linkProbability = 0.01
linkCuttingProb = 0.2
init_p_infection = 0.3
init_percent_strict= 0.8    #for regon 1 alone 
init_percent_none = 0.2  #for region 1 alone
Re = 1.34 #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4255970/
# Re = Ro *susceptible#effective reproduction rate
# Dying_Prob = np.random.exponential()                        
days = 50 #https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1541-0420.2006.00609.x
# Roe = infectionProb*Ro
#initial state of the system
simulation_length=100
total_infections = 0
infected_counts_reg1= 0
infected_counts_reg2 = 0
total_recovery = 0
total_crossing_nodes = 0
mobile_nodes_in_reg1=[]
mobile_nodes_in_reg2=[]
infected_nodes_reg1=[]
infected_nodes_reg2=[]

def update_state_based_on_neighbor(G, node_id):
    #check neighbors
    node_neighbors = list(G.neighbors(node_id))
    neighbor_status= [G.nodes[node]["status"] for node in node_neighbors]
    count_infected_neighbors = neighbor_status.count("infected")
    if len(node_neighbors) >1:
        proportion_of_infected_neighbors = count_infected_neighbors/len(node_neighbors)
        if proportion_of_infected_neighbors >= 0.5 :
            if random.random() < 0.8:
                G.nodes[node_id]["status"]="infected"
                G.nodes[node_id]["color"]="red"
                G.nodes[node_id]["counter"]+=1
              #days since infection

    return G.nodes(data=True)[node_id]

def movement(G,ref_node,crossing_nodes):

    #print(len(list(G.nodes))) 
    #if random.random() < cross_region_prob:
    updated_G = G.copy()
    #node_list = list(updated_G.nodes)
    nodes_attributes = dict(G.nodes(data=True))

    if crossing_nodes:
        #ref node is always a moving node: check the current region of the moving node
        if ref_node in crossing_nodes and G.nodes[ref_node]["cut_off"]==False:
            current_cross_node_reg = G.nodes[ref_node]["region"]
            current_nodes_attr = nodes_attributes[ref_node]
            
            print(f'current node region:{G.nodes[ref_node]["region"]}')
            #remove node from current position and region
            updated_G.remove_node(ref_node)

            #node crosses to the other region, randomly finds 20% of other nodes and reattach
            nodes_from_other_region = [node for node in G.nodes if G.nodes[node]["region"]!=current_cross_node_reg]

            #randomly select 20% node for attachment
            nodes_for_attachment= np.random.choice(nodes_from_other_region,round(0.2*len(nodes_from_other_region)), replace=False)

            #update position as the mid of its neighbors: create function for this later
            reattachment_nodes_pos= np.array([updated_G.nodes[each_node]["pos"] for each_node in nodes_for_attachment])
            #print(all_neighbors_pos)
            avg_x = reattachment_nodes_pos[:,0].mean()
            avg_y = reattachment_nodes_pos[:,1].mean()

            #reattach the node
            #first update the node pos with the new one
            current_nodes_attr["pos"]= np.array([avg_x, avg_y]) 
            #update the node region to its new region
            current_nodes_attr["region"]= G.nodes[nodes_for_attachment[0]]["region"]                           
            updated_G.add_node(ref_node, **(current_nodes_attr))
            print(f'new node region after switching:{updated_G.nodes[ref_node]["region"]}')

    else:
        #if it is not a cross node, the node moves randomly and doesnt cross the boundary
        new_x = G.nodes[ref_node]["pos"][0]+ np.random.random()
        new_y = G.nodes[ref_node]["pos"][1]+ np.random.random()
        updated_G.nodes[ref_node]["pos"]=([new_x, new_y])

    return updated_G.nodes(data=True)[ref_node]


def reconnect_node(node_id,G, mobile_nodes):
    
    if node_id not in mobile_nodes:
        node_region = G.nodes[node_id]["region"]
        #other non infected nodes from the same region
        other_non_Infected_nodes_from_same_region = [node for node in G.nodes if G.nodes[node]["region"]== node_region and\
                                                         G.nodes[node]["status"]!="infected"]
           #randomly choose connections
        potential_connections = list(np.random.choice(other_non_Infected_nodes_from_same_region, round(0.1*(len(other_non_Infected_nodes_from_same_region)))))
                                         
             #reconnect to those nodes
        for targets in potential_connections:
            G.add_edge(node_id,targets)
                                         
    else:
         for nodes in G.nodes:
            other_non_Infected_nodes = [node for node in G.nodes if G.nodes[node]["status"]!="infected"]
            choose_one_node = list(np.random.choice(other_non_Infected_nodes,1))[0]
                                    
             #randomly choose new connections in the region of chosen node
            potential_connections = [node for node in other_non_Infected_nodes if G.nodes[node]["region"]==G.nodes[choose_one_node]["region"]]
             
            #randoly select 20% of the nodes IN THE region as the chosen node
            n_selections = round(0.2*(len(potential_connections)))
            #print(f"num selections:{n_selections}")
            new_connections = list(np.random.choice(potential_connections, size=n_selections))
                                    
            for targets in new_connections:
                G.add_edge(node_id,targets)
                              
    return list(G.edges(node_id))
                               

# all_region1_nodes_list= [node for node in G.nodes() if G.nodes[node]["region"]=="1"]
# print(f"Nodes in region1: {len(all_region1_nodes_list)}")

# all_region2_nodes_list= [node for node in G.nodes() if G.nodes[node]["region"]=="2"]
# print(f"Nodes in region2: {len(all_region2_nodes_list)}")

# mobile_nodes_in_reg1 = list(np.random.choice(all_region1_nodes_list, round(0.2*len(all_region1_nodes_list)), replace=False))
# print(f"Mobile Nodes in region1: {len(mobile_nodes_in_reg1)}")

# mobile_nodes_in_reg2 = list(np.random.choice(all_region2_nodes_list, round(0.7*len(all_region2_nodes_list)), replace=False))

# print(f"Mobile Nodes in region2: {len(mobile_nodes_in_reg2)}")
# mobile_nodes = mobile_nodes_in_reg1 + mobile_nodes_in_reg2
# print(f"Mobile nodes in reg1 and 2: {mobile_nodes}")

# if mobile_nodes:
#     for node in mobile_nodes:
#         G.nodes[node]["moving"]=True                           
#         print("")

# for node in G.nodes:
#     if G.nodes[node]['region']=="1":
#         G.nodes[node]["lockdown_status"]="strict"
#     else:
#         G.nodes[node]["lockdown_status"]="none"


# for node in G.nodes:
#     if G.nodes[node]['region']=="1":
#         if random.random() < init_p_infection:
#             G.nodes[node]['status']= "infected"
#             G.nodes[node]['color']= "red"
#         else:
#             G.nodes[node]['status']= "susceptible"
#             G.nodes[node]['color']="blue"
#     else:
#         if random.random() < init_p_infection:
#             G.nodes[node]['status']= "infected"
#             G.nodes[node]['color']="red"
#         else:
#             G.nodes[node]['status']= "susceptible"
# #             G.nodes[node]['color']= "blue"

# #loop through all the nodes in the network
# for node_id, attributes in dict(G.nodes(data=True)).items():
# #if the node belongs to region 1 and is also infected
#     if attributes["region"]==1 and attributes["status"]=="infected":
#         infected_nodes_reg1.append(node_id)
#     elif attributes["region"]!=1 and attributes["status"]=="infected":
#         infected_nodes_reg2.append(node_id)

# infections_count_reg1= len(infected_nodes_reg1)
# infections_count_reg2= len(infected_nodes_reg2)
# total_infections = infections_count_reg1+infections_count_reg2

# print(f"total_infections:{total_infections}")

def initialize():
    global time, G, pos,nextNetwork,extend_,mobile_nodes_in_reg1,mobile_nodes_in_reg2,mobile_nodes,infections_count_reg1,infections_count_reg2,total_infections
    
    region1 = nx.erdos_renyi_graph(populationSize//2, 0.1, seed=0)
    for node in region1.nodes():
        region1.nodes[node]["color"]="purple"
        region1.nodes[node]["region"]="1"

    region2 = nx.barabasi_albert_graph(populationSize//2, 2, seed=0)
    for node in region2.nodes():
        region2.nodes[node]["color"]="cyan"
        region2.nodes[node]["region"]="2"

    # Combine the two regions into one graph
    G = nx.disjoint_union(region1, region2)
    
    for node in G.nodes:
        G.nodes[node]["moving"]=False
        G.nodes[node]["cut_off"]=False
        G.nodes[node]["counter"]=0
        G.nodes[node]["days_since_recovery"]=0
        
    # Assign positions to the nodes
    pos = nx.spring_layout(G)
    nx.set_node_attributes(G, pos, 'pos')
    
    all_region1_nodes_list= [node for node in G.nodes() if G.nodes[node]["region"]=="1"]
    print(f"Nodes i region1: {len(all_region1_nodes_list)}")
    all_region2_nodes_list= [node for node in G.nodes() if G.nodes[node]["region"]=="2"]
    print(f"Nodes in region2: {len(all_region2_nodes_list)}")

    mobile_nodes_in_reg1 = list(np.random.choice(all_region1_nodes_list, round(0.2*len(all_region1_nodes_list)), replace=False))
    print(f"Mobile Nodes in region1: {len(mobile_nodes_in_reg1)}")

    mobile_nodes_in_reg2 = list(np.random.choice(all_region2_nodes_list, round(0.7*len(all_region2_nodes_list)), replace=False))
    print(f"Mobile Nodes in region2: {len(mobile_nodes_in_reg2)}")
    mobile_nodes = mobile_nodes_in_reg1+mobile_nodes_in_reg2
    #print(mobile_nodes)

    if mobile_nodes:
        for node in mobile_nodes:
            G.nodes[node]["moving"]=True                           
            print("moving updated")

    for node in G.nodes:
        if G.nodes[node]['region']=="1":
            G.nodes[node]["lockdown_status"]="strict"
        else:
            G.nodes[node]["lockdown_status"]="none"


    for node in G.nodes:
        if G.nodes[node]['region']=="1":
            if random.random() < init_p_infection:
                G.nodes[node]['status']= "infected"
                G.nodes[node]['color']= "red"
                G.nodes[node]['counter']=1
            else:
                G.nodes[node]['status']= "susceptible"
                G.nodes[node]['color']="blue"
        else:
            if random.random() < init_p_infection:
                G.nodes[node]['status']= "infected"
                G.nodes[node]['color']="red"
                G.nodes[node]['counter']=1
            else:
                G.nodes[node]['status']= "susceptible"
                G.nodes[node]['color']= "blue"

    #loop through all the nodes in the network
    for node_id, attributes in dict(G.nodes(data=True)).items():
    #if the node belongs to region 1 and is also infected
        if attributes["region"]==1 and attributes["status"]=="infected":
            infected_nodes_reg1.append(node_id)
        elif attributes["region"]!=1 and attributes["status"]=="infected":
            infected_nodes_reg2.append(node_id)
        
    infections_count_reg1= len(infected_nodes_reg1)
    infections_count_reg2= len(infected_nodes_reg2)
    total_infections = infections_count_reg1+infections_count_reg2

    print(f"total_infections:{total_infections}")
    
    
    
    time = 0       
    #infections= init_p_infection/populationSize*Re #initialInfectedRatio*populationSize
    #simulation_length = 100#np.arange(0, days, 1)
    extend_=[True]
    
#     for node in G.nodes:
#         G.nodes[node]['counter']= 0  #counter counts days sinvce infection
#         #G.nodes[node]['time_since_cutoff']=0
#         G.nodes[node]['cutoff']=False
        
    nextNetwork = G.copy()
            
def observe():
    cla()
    nx.draw(G,
            pos = pos,
            node_color = [G.nodes[i]['color'] for i in G.nodes],
            cmap = cm.Wistia,
            vmin = 0,
            vmax = 1)
    axis('image')
    title('t = ' + str(time)+' '+ 'Infection_total: '+ str(total_infections) + ' '+ ' Incubation of infection: '+ str( simulation_length) +' '+ ' Effective reproduction rate: '+ str(Re) +' '+ ' Region1 infected: '+ str(infected_counts_reg1) +' '+ ' Region2 infected: '+ str(infected_counts_reg2)+' '+ 'Recovered per t: '+ str(total_recovery)+' '+ 'Total Crossing Nodes t: '+ str(total_crossing_nodes))

    

def update():
    global time, G, pos, nextNetwork,total_infections,extend_,total_recovery,total_crossing_nodes,infections_count_reg1,infections_count_reg2
#   cross_region_prob = 0.3
#     mobile_nodes_in_reg1 = list(np.random.choice(all_region1_nodes_list, round(0.2*len(all_region1_nodes_list)), replace=False))
#     print(f"Mobile Nodes in region1: {len(mobile_nodes_in_reg1)}")

#     mobile_nodes_in_reg2 = list(np.random.choice(all_region2_nodes_list, round(0.7*len(all_region2_nodes_list)), replace=False))
#     print(f"Mobile Nodes in region2: {len(mobile_nodes_in_reg2)}")
#     mobile_nodes = mobile_nodes_in_reg1.extend(mobile_nodes_in_reg2) 
    total_infections=0
    total_recovery=0
    crossing_nodes1= list(np.random.choice(mobile_nodes_in_reg1,round(0.2*len(mobile_nodes_in_reg1)), replace=False))
    crossing_nodes2= list(np.random.choice(mobile_nodes_in_reg2,round(0.2*len(mobile_nodes_in_reg2)), replace=False))
    crossing_nodes= crossing_nodes1 + crossing_nodes2 
    
    print(f"mobile nodes in 1:{mobile_nodes_in_reg1}")
    print(f"mobile nodes in 2:{mobile_nodes_in_reg2}")
    print(f"crossing nodes in 1:{crossing_nodes1}")
    print(f"crossing nodes in 2:{crossing_nodes2}")
    print(f"total crossing nodes :{crossing_nodes}")
    
    if crossing_nodes:
        total_crossing_nodes = len(crossing_nodes)     
    else:
        total_crossing_nodes=0
        
    print(f"total_crossing:{total_crossing_nodes}")
    time += 1
#     infected_nodes_reg1=[]
#     infected_nodes_reg2=[]
    #total_infections = []

    #infections = sum(infected_nodes)
    if total_infections== 0 and extend_[0]:
        simulation_length = time
        extend_[0]=False
    #'color': 'blue', 'region': '1', 'pos': array([0.17374388, 0.70365017]), 'lockdown_status': 'strict', 'status': 'susceptible
    
    for node_id in G.nodes: 
            
        if G.nodes[node_id]["status"]=="infected":
            
            if G.nodes[node_id]["counter"]>= 14:  #when infected and isolated
                nextNetwork.nodes[node_id]["status"]="recovered"
                nextNetwork.nodes[node_id]["color"]="green"
                nextNetwork.nodes[node_id]["counter"]=0
                nextNetwork.nodes[node_id]["days_since_recovery"]+=1
                total_infections-=1
                total_recovery+=1
            
            else:
                if random.random() < linkCuttingProb-0.1:
                    #print(f"link cuting prob:{linkCuttingProb-0.1}")
                    neighbors= list(G.neighbors(node_id))
                    if len(neighbors) > 0:
                        for neighbor_node in neighbors:
                            if nextNetwork.has_edge(node_id, neighbor_node):
                                nextNetwork.remove_edge(node_id, neighbor_node) #
                                #print("got here1 :cut off links")
                                nextNetwork.nodes[node_id]["cut_off"]=True
                                nextNetwork.nodes[node_id]["counter"]+=1
                else:
                    nextNetwork.nodes[node_id]["counter"]+=1
                    nextNetwork.nodes[node_id]["cut_off"]=False
                
            
#             if G.nodes[node_id]["cut_off"] and G.nodes[node_id]["status"]=="recovered":
#                 #rejoin node here
#                 #check if mobile nodes
#                 nextNetwork = reconnect_node(node_id,G, mobile_nodes_in_reg1+mobile_nodes_in_reg2).copy()
                
              
#             elif G.nodes[node_id]["counter"]>= 10:   #recover after 14 days
#                 nextNetwork.nodes[node_id]["status"]="recovered"
#                 nextNetwork.nodes[node_id]["color"]="green"
#                 nextNetwork.nodes[node_id]["counter"]=0
#                 total_infections-=1
#                 total_recovery+=1
         
#             else:
#                 if random.random() < linkCuttingProb:
#                     neighbors= list(G.neighbors(node_id))
#                     for neighbor_node in neighbors:
#                         if nextNetwork.has_edge(node_id, neighbor_node):
#                             nextNetwork.remove_edge(node_id, neighbor_node) #
#                             print("got here :cut off links")
#                             G.nodes[node_id]["cut_off"]=True               
                                                  

#                 nextNetwork.nodes[node_id]["counter"]+=1
#                 print("got here :0")
                                             
        elif G.nodes[node_id]["status"]=="susceptible":
#             print(G.nodes(data=True)[node_id])
#             print(G.nodes[node_id]["cut_off"])
#             print("got here :1")
             nextNetwork.nodes[node_id].update(update_state_based_on_neighbor(G, node_id))
           
        else:
             if G.nodes[node_id]["cut_off"]:
                 #print(f"print neighors after cutt off : {[node for node in G.neighbors(node_id)]}")
                 mobile_nodes = mobile_nodes_in_reg1 + mobile_nodes_in_reg2
                 nextNetwork.add_edges_from(reconnect_node(node_id,G,mobile_nodes))
                 nextNetwork.nodes[node_id]["days_since_recovery"]+=1
                 nextNetwork.nodes[node_id]["cut_off"]=False
                 #print(f"print neighors after reconnection : {[node for node in nextNetwork.neighbors(node_id)]}")
             
             else:  
                neighbors= list(G.neighbors(node_id))
                infected_neighbors = [node for node in neighbors if G.nodes[node]["status"]=="infected"]
                
                if infected_neighbors:
                    if random.random() < linkCuttingProb-0.1:
                        for neighbor_node in infected_neighbors:
                            if nextNetwork.has_edge(node_id, neighbor_node):
                                nextNetwork.remove_edge(node_id, neighbor_node) #remove link from infected neighbors
                                print("got here2 :cut off links")
                                nextNetwork.nodes[node_id]["cut_off"]=True


                    elif len(infected_neighbors)/len(neighbors) >= 0.1 :  #if proportion of infected is greater than 60%
                        if random.random() < 0.2:
                            nextNetwork.nodes[node_id]["status"]="infected"
                            nextNetwork.nodes[node_id]["color"]="red"
                            nextNetwork.nodes[node_id]["counter"]+=1
                            total_infections+=1
                            print("got here 3")
                                
                        else:
                            nextNetwork.nodes[node_id]["status"]="susceptible"
                            nextNetwork.nodes[node_id]["counter"]=0
                            nextNetwork.nodes[node_id]["color"]="blue"
                            print("got here 4:")
                                
#                     else:
#                         if G.nodes[node_id]["days_since_recovery"]>10:
#                             nextNetwork.nodes[node_id]["status"]="susceptible"
#                             print("recovered got here")
                        
#         print(G.nodes(data=True)[node_id])
        
    #UPDATE POS and region if the current node is a. 
    #node : i.e. part of the 20% from the strcit region and 70%

    for node in nextNetwork.nodes:
        if nextNetwork.nodes[node]["status"]=="infected":
            total_infections+=1  #from the none region
        elif nextNetwork.nodes[node]["status"]=="recovered":
            total_recovery+=1
        else:
            pass
    
    for node_id in G.nodes:
        #update nodes positions       
        if G.nodes[node_id]["moving"]==True:
            print(f"cross node: {crossing_nodes}")
            #print("get_here:last")
            print(G.nodes[node_id]["region"])
            print(G.nodes[node_id]["pos"])
            nextNetwork.nodes[node_id].update(movement(G,node_id,crossing_nodes))   
            print(nextNetwork.nodes[node_id]["region"])
            print(nextNetwork.nodes[node_id]["pos"])

                                                   
    pos = nx.spring_layout(nextNetwork, pos = pos, iterations = 2)            
                
    del G
                                                   
    G = nextNetwork.copy()       

            

pycxsimulator.GUI().start(func=[initialize, observe, update])
