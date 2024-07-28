import math
import numpy as np
def find_q_lamb (P_1, P_2,L,d,Nu,delta,ro_1,ro_2,z_1,z_2, ro_sr):
    
    lamb1 = 0.03
    lamb2 = 0.02
    j_lamb = 0
    while abs((lamb2-lamb1)/lamb1) > 0.0000001:

        lamb1 = lamb2
        # V = (abs(((P_1 - P_2) + ro_1 * 9.81 * z_1 - ro_2 * 9.81 * z_2) * d * 2)/(lamb1 * L * 1.02 * ro_sr))**(1/2)
        V = ((abs(P_1/(ro_1*9.81) - P_2/(ro_2*9.81) + z_1 - z_2) * d * 2 * 9.81)/(lamb1*L*1.02))**(1/2)
        
        Re = V*d/Nu
        
        if Re == 0:
            lamb2 = 0.02 
        elif Re < 2320:
            lamb2 = 64/Re
        elif Re < 10*d/delta*1000:
            lamb2 = 0.3164/Re**0.25
        elif Re < 500*d/delta*1000:
            lamb2 = 0.11*(delta/(d*1000) + 68/Re)**0.25
        else:
            lamb2 = 0.11*(68/Re)**0.25

        j_lamb += 1

        if j_lamb > 30:
            lamb2 = lamb1

    
    H_1 = P_1/(9.81*ro_1) + z_1
    H_2 = P_2/(ro_2*9.81) + z_2
    
    if H_1 > H_2:
        Q = V*math.pi*d**(2)/4
    else:
        Q = -V*math.pi*d**(2)/4
        
    if H_1 == H_2:
        Q = 0
    
    Q = np.round(Q,6)
    
    return Q, lamb1
