import numpy as n
import math

# connection matrix
A = n.array([[1,  0,  0,  0,  0,  0,  0], 
            [-1,  1, -1,  0,  0,  0,  0],
            [ 0, -1,  0,  1,  1,  0,  0],
            [ 0,  0,  1, -1,  0, -1,  0],
            [ 0,  0,  0,  0, -1,  1,  1],
            [ 0,  0,  0,  0,  0,  0, -1]])



# density at standard conditions
ro_stn = 890
P_stn = 101325


n_nodes = n.shape(A)[0]
n_pipes = n.shape(A)[1]
ro = n.zeros(shape=(n_nodes,1))
T = n.zeros(shape=(n_nodes,1))
z_nodes = n.zeros(shape=(n_nodes,1))  
delQ_dP = n.zeros(shape=(n_pipes,n_nodes))
P = n.zeros(shape=(n_nodes,1), dtype=int)    
delQ = n.zeros(shape=(n_nodes,1), dtype=int)    
delP = n.zeros(shape=(n_nodes,1))
lamb = n.zeros(shape=(n_pipes,1))
Q = n.zeros(shape=(n_pipes,1))
T_avg = n.zeros(shape=(n_pipes,1))
P_avg = n.zeros(shape=(n_pipes,1))
ro_avg = n.zeros(shape=(n_pipes,1))
nu_avg = n.zeros(shape=(n_pipes,1))




known_nodes = [0, 5]

P[known_nodes[0]] = 190 * 10**5
P[known_nodes[1]] = 120 * 10**5             



T[known_nodes[0]] = 23 #temperature at known nodes OLD

unknown_nodes = [i for i in range(n_nodes)]

for i in known_nodes:
    unknown_nodes.remove(i)
    
    
# parameter for thermal calculation
T_ground = -10
K_t = 10**9
K_temp = 1.7
kof_v_rs = 2 * 10**(-6)
Cp = 2300


d = n.zeros(shape=(n_pipes,1)) #pipe diameters
delta = n.zeros(shape=(n_pipes,1)) #pipe roughness
L = n.zeros(shape=(n_pipes,1)) #pipe length


for i in range(n_nodes):
    T[i] = 23 - 0.1 * i
    z_nodes[i] = 0
    ro[i] = 890


for i in range(n_pipes):
    delta[i] = 0.1 * 10**(-3)  #pipe roughnesses
    

d[0] = 0.516 #pipe diameters
d[1] = 0.410 #pipe diameters
d[2] = 0.311 #pipe diameters
d[3] = 0.313 #pipe diameters
d[4] = 0.311 #pipe diameters
d[5] = 0.410 #pipe diameters
d[6] = 0.618 #pipe diameters


L[0] = 900 * 10^2 #pipe length
L[1] = 170 * 10^2 #pipe length
L[2] = 80 * 10^2 #pipe length
L[3] = 120 * 10^2 #pipe length
L[4] = 150 * 10^2 #pipe length
L[5] = 230 * 10^2 #pipe length
L[6] = 450 * 10^2 #pipe length



# pressure in unknown nodes
for i in unknown_nodes:
      P[i] = 18 * 10**5 - i * 10**5

P_for_iter = P.copy()
SHU = n.zeros(shape=(1,n_pipes))
At = n.transpose(A)

nu_20 = 5 * 10**(-6) #kinematic viscosity at 20 degrees
nu_0 = 8 * 10**(-6) #kinematic viscosity at 0 degrees
k_nu = math.log(nu_20/nu_0)/(0-20)
T_for_calc = n.zeros(shape=(n_nodes,1))
