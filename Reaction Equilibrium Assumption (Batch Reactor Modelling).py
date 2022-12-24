#------------------------------------------------------------------------------------------------------
#  _____                 _   _               ______            _ _ _ _          _                 
# |  __ \               | | (_)             |  ____|          (_) (_) |        (_)                
# | |__) |___  __ _  ___| |_ _  ___  _ __   | |__   __ _ _   _ _| |_| |__  _ __ _ _   _ _ __ ___  
# |  _  // _ \/ _` |/ __| __| |/ _ \| '_ \  |  __| / _` | | | | | | | '_ \| '__| | | | | '_ ` _ \ 
# | | \ \  __/ (_| | (__| |_| | (_) | | | | | |___| (_| | |_| | | | | |_) | |  | | |_| | | | | | |
# |_|  \_\___|\__,_|\___|\__|_|\___/|_| |_| |______\__, |\__,_|_|_|_|_.__/|_|  |_|\__,_|_| |_| |_|
#                                                     | |                                         
#                                                     |_|   _                                     
#                   /\                                 | | (_)                                    
#                  /  \   ___ ___ _   _ _ __ ___  _ __ | |_ _  ___  _ __                          
#                 / /\ \ / __/ __| | | | '_ ` _ \| '_ \| __| |/ _ \| '_ \                         
#                / ____ \\__ \__ \ |_| | | | | | | |_) | |_| | (_) | | | |                        
#               /_/    \_\___/___/\__,_|_| |_| |_| .__/ \__|_|\___/|_| |_|                        
#                                                | |                                              
#                                                |_|                                              
#------------------------------------------------------------------------------------------------------
#
#  Code Developped & Tested by Monica-Iulia Lacatus (MSc Chemical Engineering TU Delft)
#
#------------------------------------------------------------------------------------------------------
# 
#  Example: Series Reaction A <-> B <-> C with reaction rate constants (k1,k-1,k2,k-2) in batch reactor
#
#
#  The following code aims to demonstrate the use of the reaction equilibirum assumption (REA):
#    
#  This is an assumption often used to simpliy reaction models. In the case of the series reaction in 
#  this example if k2,k-2 >> k1,k-1  meaning that the second reaction is fast enough to equilibrate 
#  immediately after any displacement from its equilibirum condition, than the full reaction model can 
#  be simplified on the basis that rate of B <-> C (r2) is equal to 0.
#
#------------------------------------------------------------------------------------------------------
#                                 IMPORTING LIBRARIES
#------------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings("ignore")

plt.rcParams["axes.labelsize"] = 14
plt.rcParams['font.family'] = 'Times New Roman'

#----------------------------------------------------------------------------------------------------
#                       Part 1 : Solving Full Reaction Model ODE Equations
#----------------------------------------------------------------------------------------------------

def series_reactions_batch(t,c_ini,k):
    
    # Batch reactor component species balances used to model series reaction.
    
    CA = c_ini[0]
    CB = c_ini[1]
    CC = c_ini[2] 
    
    k1  = k[0] 
    k_1 = k[1] 
    k2  = k[2] 
    k_2 = k[3]
    
    r1 = k1*CA - k_1*CB
    r2 = k2*CB - k_2*CC
    
    dCAdt = -r1
    dCBdt = r1 -r2
    dCCdt = r2
    
    return [dCAdt,dCBdt,dCCdt]

 
# For reaction A <-> B :
k1  = 1     #[1/s]
k_1 = 0.5   #[1/s]

# For reaction B <-> C :
k2  = 1     #[1/s]
k_2 = 1     #[1/s]

k = [k1, k_1, k2, k_2]

# Time length for solution:
t = np.linspace(0,5,110)

# Initial concentrations of A,B and C:
CA0 = 2      #[mol/L]
CB0 = 0.8    #[mol/L]
CC0 = 0      #[mol/L]

c_ini = [CA0, CB0, CC0]

# Solving the ODE reaction model describing the reaction series:
conc_sol = solve_ivp(series_reactions_batch, [0, t[-1]], c_ini, t_eval=t, method='Radau', args=[k])
CA=conc_sol.y[0]
CB=conc_sol.y[1]
CC=conc_sol.y[2]

# Plotting Results:
fig = plt.figure(figsize=(9,6))
plt.plot(t,CA,'k',label='A')
plt.plot(t,CB,'r',label='B')
plt.plot(t,CC,'g',label='C')
plt.xlabel('Time [s]')
plt.ylabel('Concentration [mol/L]')
plt.title('Concentration Evolution of Species in Batch Reactor',fontsize = 15)
plt.grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
plt.legend(loc='upper right',shadow='True',fontsize=14)
plt.show()


#----------------------------------------------------------------------------------------------------
#           Part 2 : Solving Full Reaction Model ODE Equations With Increasing k2,k_2 values
#----------------------------------------------------------------------------------------------------

# Case 1:
k2_10  = 10     #[1/s]
k_2_10 = 10     #[1/s]
k_10 = [k1, k_1, k2_10, k_2_10]

conc_sol_10 = solve_ivp(series_reactions_batch, [0, t[-1]], c_ini, t_eval=t, method='Radau', args=[k_10])
CA_10=conc_sol_10.y[0]
CB_10=conc_sol_10.y[1]
CC_10=conc_sol_10.y[2]

# Case 2:
k2_100  = 100     #[1/s]
k_2_100 = 100     #[1/s]
k_100 = [k1, k_1, k2_100, k_2_100]

conc_sol_100 = solve_ivp(series_reactions_batch, [0, t[-1]], c_ini, t_eval=t, method='Radau', args=[k_100])
CA_100=conc_sol_100.y[0]
CB_100=conc_sol_100.y[1]
CC_100=conc_sol_100.y[2]

# Case 3:
k2_1000  = 1000     #[1/s]
k_2_1000 = 1000     #[1/s]
k_1000 = [k1, k_1, k2_1000, k_2_1000]

conc_sol_1000 = solve_ivp(series_reactions_batch, [0, t[-1]], c_ini, t_eval=t, method='Radau', args=[k_1000])
CA_1000=conc_sol_1000.y[0]
CB_1000=conc_sol_1000.y[1]
CC_1000=conc_sol_1000.y[2]

# Case 4:
k2_inf  = np.inf     #[1/s]
k_2_inf = np.inf    #[1/s]
k_inf = [k1, k_1, k2_1000, k_2_1000]

conc_sol_inf = solve_ivp(series_reactions_batch, [0, t[-1]], c_ini, t_eval=t, method='Radau', args=[k_inf])
CA_inf=conc_sol_inf.y[0]
CB_inf=conc_sol_inf.y[1]
CC_inf=conc_sol_inf.y[2]


# Plotting Results:
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(20, 15), dpi=150) 
fig.suptitle('Concentration Evolution of Species in Batch Reactor With Increasing k$_2$, k$_{-2}$',fontsize = 35,y=0.95)
axs[0,0].plot(t,CA,'k-',label='A| k$_2$ = k$_{-2}$ = 1')
axs[0,0].plot(t,CB,'r-',label='B| k$_2$ = k$_{-2}$ = 1')
axs[0,0].plot(t,CC,'g-',label='C| k$_2$ = k$_{-2}$ = 1')
axs[0,0].plot(t,CA_10,'k--',label='A| k$_2$ = k$_{-2}$ = 10')
axs[0,0].plot(t,CB_10,'r--',label='B| k$_2$ = k$_{-2}$ = 10')
axs[0,0].plot(t,CC_10,'g--',label='C| k$_2$ = k$_{-2}$ = 10')
axs[0,0].set_xlabel('Time [s]',fontsize = 20)
axs[0,0].set_ylabel('Concentration [mol/L]',fontsize = 20)
axs[0,0].set_title('k$_2$ = k$_{-2}$ = 10 [1/s]',fontsize = 20)
axs[0,0].grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
axs[0,0].legend(loc='upper right',shadow='True',fontsize=14)  

axs[0,1].plot(t,CA,'k-',label='A| k$_2$ = k$_{-2}$ = 1')
axs[0,1].plot(t,CB,'r-',label='B| k$_2$ = k$_{-2}$ = 1')
axs[0,1].plot(t,CC,'g-',label='C| k$_2$ = k$_{-2}$ = 1')
axs[0,1].plot(t,CA_100,'k--',label='A| k$_2$ = k$_{-2}$ = 100')
axs[0,1].plot(t,CB_100,'r--',label='B| k$_2$ = k$_{-2}$ = 100')
axs[0,1].plot(t,CC_100,'g--',label='C| k$_2$ = k$_{-2}$ = 100')
axs[0,1].set_xlabel('Time [s]',fontsize = 20)
axs[0,1].set_ylabel('Concentration [mol/L]',fontsize = 20)
axs[0,1].set_title('k$_2$ = k$_{-2}$ = 100 [1/s]',fontsize = 20)
axs[0,1].grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
axs[0,1].legend(loc='upper right',shadow='True',fontsize=14)  

axs[1,0].plot(t,CA,'k-',label='A| k$_2$ = k$_{-2}$ = 1')
axs[1,0].plot(t,CB,'r-',label='B| k$_2$ = k$_{-2}$ = 1')
axs[1,0].plot(t,CC,'g-',label='C| k$_2$ = k$_{-2}$ = 1')
axs[1,0].plot(t,CA_1000,'k--',label='A| k$_2$ = k$_{-2}$ = 1000')
axs[1,0].plot(t,CB_1000,'r--',label='B| k$_2$ = k$_{-2}$ = 1000')
axs[1,0].plot(t,CC_1000,'g--',label='C| k$_2$ = k$_{-2}$ = 1000')
axs[1,0].set_xlabel('Time [s]',fontsize = 20)
axs[1,0].set_ylabel('Concentration [mol/L]',fontsize = 20)
axs[1,0].set_title('k$_2$ = k$_{-2}$ = 1000 [1/s]',fontsize = 20)
axs[1,0].grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
axs[1,0].legend(loc='upper right',shadow='True',fontsize=14)   

axs[1,1].plot(t,CA,'k-',label='A| k$_2$ = k$_{-2}$ = 1')
axs[1,1].plot(t,CB,'r-',label='B| k$_2$ = k$_{-2}$ = 1')
axs[1,1].plot(t,CC,'g-',label='C| k$_2$ = k$_{-2}$ = 1')
axs[1,1].plot(t,CA_inf,'k--',label='A| k$_2$ = k$_{-2}$ = $\infty$')
axs[1,1].plot(t,CB_inf,'r--',label='B| k$_2$ = k$_{-2}$ = $\infty$')
axs[1,1].plot(t,CC_inf,'g--',label='C| k$_2$ = k$_{-2}$ = $\infty$')
axs[1,1].set_xlabel('Time [s]',fontsize = 20)
axs[1,1].set_ylabel('Concentration [mol/L]',fontsize = 20)
axs[1,1].set_title('k$_2$ = k$_{-2}$ = $\infty$',fontsize = 20)
axs[1,1].grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
axs[1,1].legend(loc='upper right',shadow='True',fontsize=14)   
plt.show()

#----------------------------------------------------------------------------------------------------
#          Part 3 : Applying REA Approximation to Obtain Simplified ODE Reaction Model
#----------------------------------------------------------------------------------------------------

def series_reactions_batch_REA(t,c_ini,k):
    
    # Batch reactor component species balances obtained using REA approximation
    
    CA = c_ini[0]
    CB = c_ini[1]
    CC = c_ini[2] 
    
    k1  = k[0] 
    k_1 = k[1] 
    K2  = k[2] 
    
    r1 = k1*CA - k_1*CB
    
    dCAdt = -r1
    dCBdt = (1/(1+K2))*r1
    dCCdt = (K2/(1+K2))*r1
    
    return [dCAdt,dCBdt,dCCdt]

# Equilibrium constant of reaction B <-> C:
K2  = 1  #k2/k-2

k = [k1, k_1, K2]

# Defining new intial concentrations:

CBs = (1/(1+K2))*(CB0+CC0)
CCs = (K2/(1+K2))*(CB0+CC0)
c_ini = [CA0, CBs, CCs]

conc_sol_REA = solve_ivp(series_reactions_batch_REA, [0, t[-1]], c_ini, t_eval=t, method='Radau', args=[k])
CA_REA=conc_sol_REA.y[0]
CB_REA=conc_sol_REA.y[1]
CC_REA=conc_sol_REA.y[2]

# Plotting Results:
fig = plt.figure(figsize=(9,6))
plt.plot(t,CA_REA,'k',label='A')
plt.plot(t,CB_REA,'r',label='B')
plt.plot(t,CC_REA,'g',label='C')
plt.xlabel('Time [s]')
plt.ylabel('Concentration [mol/L]')
plt.title('Concentration Evolution of Species in Batch Reactor (REA Approximation)',fontsize = 15)
plt.grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
plt.legend(loc='upper right',shadow='True',fontsize=14)
plt.show()


#----------------------------------------------------------------------------------------------------
#               Part 4 : Comparison of REA Approximation to Full Reaction Model
#----------------------------------------------------------------------------------------------------

# Plotting Comparison:
fig = plt.figure(figsize=(9,6))
plt.plot(t,CA_REA,'k',label='A| REA approx.')
plt.plot(t,CB_REA,'r',label='B| REA approx.')
plt.plot(t,CC_REA,'g',label='C| REA approx.')
plt.plot(t,CA_10,'k--',label='A| k$_2$ = k$_{-2}$ = 10')
plt.plot(t,CB_10,'r--',label='B| k$_2$ = k$_{-2}$ = 10')
plt.plot(t,CC_10,'g--',label='C| k$_2$ = k$_{-2}$ = 10')
plt.xlabel('Time [s]')
plt.ylabel('Concentration [mol/L]')
plt.title('Comparison of REA Approximation to Full Reaction Model',fontsize = 15)
plt.grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
plt.legend(loc='upper right',shadow='True',fontsize=12)
plt.show()