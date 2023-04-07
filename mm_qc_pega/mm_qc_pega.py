import yaml
import numpy as np 
import galois 
import datetime
import time 
import graphs 
import peg_utils
import os

#Load settings
try:
    with open("config.yml", "r") as ymlfile:
        settings = yaml.safe_load(ymlfile)
        data = np.load(settings["input_fname"])
except:
    print("Config or input data not found, using default settings.")
    settings = {
                "n" : 73,
                "m" : 9,
                "scaling_factor" : 4,
                "girth_search_depth" : 1,
                "seed" : 'None',
                "gcd" : False,
                "input_fname" : None,
                "output_fname" : "test.qc"
                }
    data = {
            "vn_degree_node" : np.array([0,1]),
            "cn_degree_node" : np.array([0,0,0,0,0,0,0,0,0,1/9,0,0,0,0,0,0,8/9])
            } 


header = f"""
===================================================================
Creating ({settings["m"]},{settings["n"]},{settings["scaling_factor"]}) code with the {settings["girth_search_depth"]}-edge QC-PEGA algorithm for 
variable node edge distribution.
Input file: {settings["input_fname"]}
Output file: {settings["output_fname"]}
===================================================================
"""
print(header)
G = graphs.QC_tanner_graph(settings["m"], settings["n"], settings["scaling_factor"])
#vn_degrees = peg_utils.to_degree_distribution(data["vn_degree_node"], G.n_vn)
#cn_degrees = peg_utils.to_degree_distribution(data["cn_degree_node"], G.n_cn)

# Get the node perspective degree distributions from density evolution optimization
lam_node = data["lam_node"]
rho_node = data["rho_node"]
vn_degrees = np.flip(peg_utils.to_degree_distribution(lam_node,G.n_vn))
cn_degrees = peg_utils.to_degree_distribution(rho_node,G.n_cn)

#vn_degrees = np.repeat(vn_degrees, settings["scaling_factor"])
#cn_degrees = np.repeat(cn_degrees, settings["scaling_factor"])

if not settings["seed"] == 'None':
    np.random.seed(int(settings["seed"]))
t0 = time.time()

for current_vn_index in range(0, G.n_vn, G.N):

    d = vn_degrees[current_vn_index]
    for k in range(1, int(d+1)):
        rk = int(min(settings["girth_search_depth"], d - k +1))
        max_girth, cn_girths = peg_utils.rk_edge_local_girth(G, current_vn_index, rk, cn_degrees)
    
        ci = peg_utils.strategy1(max_girth, cn_girths, G, current_vn_index, cn_degrees)
        G.add_cyclical_edge_set(ci, current_vn_index)

    dt = time.time() - t0
    completed_progress = (current_vn_index+G.N)/(G.n_vn)
    time_left = dt/completed_progress - dt        
    status = f"{100*completed_progress:.2f}% completed, elapsed time {int(dt//60)}:{dt% 60:.2f}s. Approximate time left: {int(time_left//60)}:{time_left % 60:.2f}"
    print(status, end='\r', flush=True)
    

print("")
print(f"Edge growth finsihed. Total elapsed time: {int(dt//60)} minutes, {dt % 60:.2f} seconds.")

print(np.bincount(G.get_check_degrees()))
print(np.bincount(cn_degrees))
print(np.bincount(G.get_var_degrees()))
print(np.bincount(vn_degrees))

G_reordered = peg_utils.make_invertable(G)
print("Matrix invertion sucessfull, saving file...")

metadata = f"""
# QC-code generated with the r-edge metric constrained QC-PEG algorithm.
# 
# Name   : {settings["output_fname"]}
# R      : {settings["girth_search_depth"]}
# Metric : Minimum distance
# Date   : {datetime.datetime.now()}
"""
os.makedirs('data/', exist_ok=True)
G_reordered.save(settings["output_fname"])
f = open(settings["output_fname"], 'a')
f.write(metadata)
f.close()

print("MM-QC-PEGA completed.")

peg_utils.graph_stats(G,settings["output_fname"])