import numpy as np
import math
import Find_Q_Lambd as func
from Inish_data import *
import pandas as pd

# node connection vector
b = np.zeros(shape=(n_pipes,2), dtype=int)  
for i in range(n_pipes):         
      for j in range(n_nodes):
            if A[j][i] == 1:
                  b[i][0] = j
            elif A[j][i] == -1:
                  b[i][1] = j


# variables to calculate

# variables to calculate
T_for_iter = T.copy()
P = P_for_iter.copy()

T_for_calc = np.zeros(shape=(n_nodes,1))
T = T_for_iter.copy()
T_iter = 0


while sum(abs(T_for_calc - T)) > n_nodes*0.00001:
#The calculation begins, it ends when the temperatures in the nodes are correctly calculated
      
      if T_iter > 0:
            T = T_for_calc
      
      #Calculate average temperature, pressure, density and kinematic viscosity of each pipeline
      for i in range(np.shape(b)[0]):
            T_avg[i] = T_ground + (T[b[i,0]] - T[b[i,1]]) / (math.log((T[b[i,0]] - T_ground)/(T[b[i,1]] - T_ground)))
            P_avg[i] = (P[b[i,0]] + P[b[i,1]]) / 2
            ro_avg[i] = ro_stn * (1 + kof_v_rs * (20-T_avg[i]) + ((P_avg[i] - P_stn))/K_t)
            nu_avg[i] = nu_20 * math.exp(-k_nu*(T_avg[i] - 20))

      #We count the density at each node
      for i in range(n_nodes):
            ro[i] = ro_stn * (1 + kof_v_rs * (20-T[i]) + (P[i] - P_stn)/K_t)


      delQ_dP = np.zeros(shape=(n_pipes,n_nodes))
      delQ = np.zeros(shape=(n_nodes,1), dtype=int)
      # set the starting value of the flow error in the nodes
      for i in range(n_nodes):
            delQ[i] = 1.5
      
      # calculating flow value in all pipelines is begun
      while sum(abs(delQ)) > 0.0001:     
            delP = np.dot(At, P)
      
            for i in range(n_pipes):
                  Q[i], lamb[i] = func.find_q_lamb ( P[b[i][0]], P[b[i][1]], L[i], d[i], nu_avg[i], delta[i], ro[b[i][0]], ro[b[i][1]], z_nodes[b[i][0]], z_nodes[b[i][1]], ro_avg[i])
            
            # find value of the flow error in the each node
            delQ = np.dot(A, Q)
            
            #We begin determining pressure corrections in each node based on the flow error
            z_ro_g = (z_nodes.transpose()*ro.transpose()) * 9.81 
            k = np.ones(shape=(n_pipes,2))
            for i in range(n_pipes):
                  for f in range(2):
                        k[i,f] = ( math.pi * d[i]**2.5) / (math.sqrt(2 * 1.02 * lamb[i] * L[i]) * ro[b[i,f]])  / (2 * math.sqrt(abs( (P[b[i,0]]+z_ro_g[0, b[i,0]])/ro[b[i,0]] - (P[b[i,1]]+z_ro_g[0, b[i,1]])/ro[b[i,1]] )))

            matrix = np.zeros(shape=(n_nodes,n_nodes))

            # find the derivative
            for i in range(n_nodes):
                  st = A[i]
                  for j in range(n_pipes):
                        line_sign = st[j]
                        if line_sign == 0:
                              continue
                        else:
                              column_sign = A[:,j]
                              for h in range(n_nodes):
                                    if column_sign[h] == 0:
                                          continue
                                    sign = column_sign[h] * line_sign
                                    if column_sign[h] == 1:
                                          matrix[i,h] = matrix[i,h] +k[j,0]*sign
                                    else:
                                          matrix[i,h] = matrix[i,h] +k[j,1]*sign

            #removing known nodes from the matrix
            for i in range(len(known_nodes)):
                  matrix = np.delete(matrix, (known_nodes[i]-i), 0)
                  matrix = np.delete(matrix, (known_nodes[i]-i), 1)
                  delQ = np.delete(delQ, (known_nodes[i]-i), 0)   
            
            #find the corrections
            dP = np.linalg.lstsq(matrix, delQ)[0]    
            
            #correct the pressure in the node by the correction value
            for i in range(len(unknown_nodes)):
                  P[unknown_nodes[i]] = P[unknown_nodes[i]] - dP[i]

            #Pressure value sign check
            for i in range(n_nodes):
                  if  P[i] < 0:
                        print('!!!')
      

      #We calculate the Shukhov coefficient, and then we calculate the new temperature at each node
      for i in range(n_pipes):
            SHU[0, i] = ((d[i] * (-math.pi) * K_temp) / (Q[i] * ro_avg[i] * Cp)) * (L[i])
      
      T_for_calc[0] = T[0]
      for i in range(1, n_nodes):
            for j in range(n_pipes):
                  if T[b[j,1]] == T[i]:
                        if Q[j] > 0:
                              T_for_calc[b[j,1]] = T_ground + (T_for_calc[b[j,0]] - T_ground) * math.exp(SHU[0,j])
                        else: 
                              T_for_calc[b[j,1]] = T_ground + (T_for_calc[b[j,0]] - T_ground) / math.exp(SHU[0,j])
                        break
                  elif T[b[j,0]] == T[i]:
                        if Q[j] > 0:
                              T_for_calc[b[j,0]] = T_ground + (T_for_calc[b[j,1]] - T_ground) / math.exp(SHU[0,j]) 
                        else:
                              T_for_calc[b[j,0]] = T_ground + (T_for_calc[b[j,1]] - T_ground) * math.exp(SHU[0,j])
                        break

      T_iter += 1



#save results
Save_P = pd.DataFrame(data=P, index=[f'P{n}' for n in range(1, n_nodes+1)], columns=['Pressure in nodes, Pa'])
Save_T = pd.DataFrame(data=T, index=[f'T{n}' for n in range(1, n_nodes+1)], columns=['Temperature in nodes, °C'])
Save_Q = pd.DataFrame(data=abs(Q), index=[f'Q {b[n,0]+1}-{b[n,1]+1}' for n in range(n_pipes)], columns=['Flow in pipes, m^3/s'])

writer = pd.ExcelWriter('Results.xlsx')
Save_P.to_excel(writer, sheet_name='Results', startcol=1, startrow=1)
Save_Q.to_excel(writer, sheet_name='Results', startcol=4, startrow=1)
Save_T.to_excel(writer, sheet_name='Results', startcol=8, startrow=1)

writer._save()