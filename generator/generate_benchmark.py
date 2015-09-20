#! /usr/bin/python

import sys
import time
import os 
import random
import math



areas = []
densities = []
sink_densities = []
path = None
grid = []
per_relays = []
clusters = []
lambdas = []
uniform = False
n_instances = 1

SEED_LOW = 10000
SEED_HIGH = 99999
def rand_int():
	return random.randint(SEED_LOW,SEED_HIGH)




if len(sys.argv) != 2:
	print "Usage: prog <config_file>"
	sys.exit(-1)

config_file = sys.argv[1]
f = open(config_file)

lines = f.readlines()

# Default values for PRR metric
ttx = 2.533
demand = 0.005

tx_range = 7.0

for line in lines:
	s = line.split()
	if len(s) == 0:
		continue
	if s[0] == "areas":
		for i in range(1,len(s),2):
			areas.append((s[i],float(s[i+1])))
	if s[0] == "densities":
		for i in range(1,len(s)):
			densities.append(float(s[i]))
	if s[0] == "sink_densities":
		for i in range(1,len(s)):
			sink_densities.append(float(s[i]))
	if s[0] == "path":
		path = s[1]
	if s[0] == "grid":
		for i in range(1,len(s),2):
			grid.append((s[i],float(s[i+1])))
	if s[0] == "relays":
		for i in range(1,len(s),2):
			per_relays.append((s[i],float(s[i+1])))
	if s[0] == "cluster":
		cluster = list()
		for i in range(2,len(s)):
			cluster.append(float(s[i]))
		clusters.append((s[1],cluster))
	if s[0] == "uniform":
		uniform = True
	if s[0] == "seed":
		seed = int(s[1])
	if s[0] == "lambdas":
		for i in range(1,len(s),2):
			lambdas.append((s[i],float(s[i+1])))
	if s[0] == "model":
		model = s[1]
	if s[0] == "number":
		n_instances = int(s[1])
	if s[0] == "ttx":
		ttx = float(s[1])
	if s[0] == "demand":
		demand = float(s[1])
	if s[0] == "tx_range":
		tx_range = float(s[1])
	if s[0] == "X":
		area_x = int(s[1])
	if s[0] == "Y":
		area_y = int(s[1])



main_str = "[Main]\nsolve = %d\ngenerate_instance = %d\n"
prob_str = "\n[ProblemInstance]\ngrid_res = %.2f\ninstance_id = %s\nstatic = %d\nmax_relays = %d\nbase = %d\nX = %.2f\nY = %.2f\ninstance_path = %s\n"
prob_const = "tx_range = %f\nmodel = %s\ninput_mps = default\nuse_mps = 1\n"%(tx_range,model)
param_str = "\n[Parameters]\noutput_path = %s\ngeneration_path = %s\nseed = %d\ngrid_strategy = %s\n"
param_const = "use_gurobi = 0\nuse_cplex = 1\nglpk_out = 1\nsolution_strategy = default\noutput_sol = 1\noutput_sim = 1\noutput_mps_file = default\nuse_heuristic = 0\nverbose = 2\ndebug_level = 0\ngen_graphs = 1\n"
prr_metric = "\n[PRR_METRIC]\nlambda = %.3f\nsigma = 4.0\nple = 3.0\n" + "ttx = %.4f\ndemand = %.4f\n"%(ttx,demand)
cluster_str = "\n[Clusters]\nuse_clusters = %d\nX = %d\nY = %d\n"
grid_h2_str = "\n[Grid_H2]\nmin_res = %.2f\nmax_res = %.2f\n"



random.seed(seed)




# Create dirs
for a in areas:
	dir1_name = "%s"%(a[0])
	os.mkdir("%s/%s"%(path,dir1_name))
	for d in range(len(densities)):
		dir2_name = "D%d"%(d)
		os.mkdir("%s/%s/%s"%(path,dir1_name,dir2_name))
		for sd in range(len(sink_densities)):
			dir3_name = "B%d"%(sd)
			os.mkdir("%s/%s/%s/%s"%(path,dir1_name,dir2_name,dir3_name))
			if uniform == True:
				nclusters = len(clusters) + 1
			else:
				nclusters = len(clusters)
			for cl in range(nclusters):
				if cl == len(clusters) and uniform:
					dir4_name = "C0"
				else:
					cname = clusters[cl][0]
					dir4_name = "C%s"%(cname)
				os.mkdir("%s/%s/%s/%s/%s"%(path,dir1_name,dir2_name,dir3_name,dir4_name))
				for rk in range(len(per_relays)+1):
					if rk == len(per_relays):
						dir5_name = "K0"
					else:
						rname = per_relays[rk][0]
						dir5_name = "K%s"%(rname)
					instance_path = "%s/%s/%s/%s/%s/%s/sol/"%(path,dir1_name,dir2_name,dir3_name,dir4_name,dir5_name)
					ini_path = "%s/%s/%s/%s/%s/%s"%(path,dir1_name,dir2_name,dir3_name,dir4_name,dir5_name) 
					os.mkdir("%s/%s/%s/%s/%s/%s"%(path,dir1_name,dir2_name,dir3_name,dir4_name,dir5_name))
					os.mkdir(instance_path)
					statics = math.ceil(1.0*(a[1]**2)*densities[d])
					bases = math.ceil(1.0*(a[1]**2)*sink_densities[sd])
					content = ""
					
					if rk == len(per_relays):
						for n_i in range(n_instances):
							instance_seed = rand_int()
							for l in range(len(lambdas)):
								content = main_str%(1,1)
								instance_id = "%s-%s-%s-%s-%s-n%d-L%s"%(dir1_name,dir2_name,dir3_name,dir4_name,dir5_name,n_i,lambdas[l][0])
								content = content + prob_str%(0,instance_id,statics,0,bases,a[1],a[1],instance_path)
								content = content + prob_const
								content = content + param_str%(instance_path,instance_path,instance_seed,"default")
								content = content + param_const
								content = content + prr_metric%(lambdas[l][1])
								if cl == len(clusters) and uniform:
									content = content + cluster_str%(0,0,0)
								else:
									clu = clusters[cl][1]
									content = content + cluster_str%(1,len(clu)**0.5,len(clu)**0.5)
									clu_d = "density = " + str(clu[0])
									for ccd in clu[1:]:
										clu_d = clu_d + "," + str(ccd)
									clu_d = clu_d + "\n"
									content = content + clu_d

								f = open(ini_path + "/" + instance_id+".ini","w")
								f.write(content)
								f.close()

							print "Generated ", instance_id
					else:
						relays = math.ceil((statics + bases)*per_relays[rk][1])
						for res in range(len(grid)+1):
							content = main_str%(1,1)

							if res != len(grid):
								instance_id = "%s-%s-%s-%s-%s-R%s"%(dir1_name,dir2_name,dir3_name,dir4_name,dir5_name,grid[res][0])
								content = content + prob_str%(grid[res][1],instance_id,statics,relays,bases,a,a,instance_path)
								content = content + prob_const
								content = content + param_str%(instance_path,instance_path,instance_seed,"default")
								content = content + param_const
						
							else:
								instance_id = "%s-%s-%s-%s-%s-RX"%(dir1_name,dir2_name,dir3_name,dir4_name,dir5_name)
								content = content + prob_str%(0,instance_id,statics,relays,bases,a,a,instance_path)
								content = content + prob_const
								content = content + param_str%(instance_path,instance_path,instance_seed,"grid_h2")
								content = content + param_const
								content = content + grid_h2_str%(grid[0][1],grid[len(grid)-1][1])
								

								
						

							if cl == len(clusters):
								content = content + cluster_str%(0,0,0)
							else:
								clu = clusters[cl][1]
								content = content + cluster_str%(1,len(clu)**0.5,len(clu)**0.5)
								clu_d = "density = " + str(clu[0])
								for ccd in clu[1:]:
									clu_d = clu_d + "," + str(ccd)
								clu_d = clu_d + "\n"
								content = content + clu_d
				
							print "Generated ", instance_id
							f = open(ini_path + "/" + instance_id+".ini","w")
							f.write(content)
							f.close()
				
					

	







	

			
