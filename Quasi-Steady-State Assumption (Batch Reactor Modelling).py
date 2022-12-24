#-----------------------------------------------------------------------------------------------------------------
#   ____                  _             _____ _                 _                   _____ _        _       
#  / __ \                (_)           / ____| |               | |                 / ____| |      | |      
# | |  | |_   _  __ _ ___ _   ______  | (___ | |_ ___  __ _  __| |_   _   ______  | (___ | |_ __ _| |_ ___ 
# | |  | | | | |/ _` / __| | |______|  \___ \| __/ _ \/ _` |/ _` | | | | |______|  \___ \| __/ _` | __/ _ \
# | |__| | |_| | (_| \__ \ |           ____) | ||  __/ (_| | (_| | |_| |           ____) | || (_| | ||  __/
#  \___\_\\__,_|\__,_|___/_|          |_____/ \__\___|\__,_|\__,_|\__, |          |_____/ \__\__,_|\__\___|
#                                                                  __/ |                                   
#                                                              _  |___/                                    
#                          /\                                 | | (_)                                      
#                         /  \   ___ ___ _   _ _ __ ___  _ __ | |_ _  ___  _ __                            
#                        / /\ \ / __/ __| | | | '_ ` _ \| '_ \| __| |/ _ \| '_ \                           
#                       / ____ \\__ \__ \ |_| | | | | | | |_) | |_| | (_) | | | |                          
#                      /_/    \_\___/___/\__,_|_| |_| |_| .__/ \__|_|\___/|_| |_|                          
#                                                       | |                                                
#                                                       |_|                                                
#------------------------------------------------------------------------------------------------------------------
#
#  Code Developped & Tested by Monica-Iulia Lacatus (MSc Chemical Engineering TU Delft)
#
#------------------------------------------------------------------------------------------------------------------
#  Example: Series Reaction A -> B -> C with reaction rate constants (k1,k2) in batch reactor
#
#  The following code aims to demonstrate the use of the quasi-steady-state assumption (QSSA):
#    
#  This is an assumption often used to simpliy reaction models. In the case of the series reaction in 
#  this example, if B is a highly reactive intermediate which is able to equilibrate quickly to its 
#  quasi-steady state value, its production rate can be set equal to zero. This simplifies the chemical
#  reaction model significantly as proven by the code below.
#
#------------------------------------------------------------------------------------------------------------------
#                                 Importing Libraries
#------------------------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings("ignore")

plt.rcParams["axes.labelsize"] = 14
plt.rcParams['font.family'] = 'Times New Roman'

#------------------------------------------------------------------------------------------------------------------
#                        Part 1 : Solving Full Reaction Model ODE Equations
#------------------------------------------------------------------------------------------------------------------

def series_reactions_batch(t,c_ini,k):
    
    # Batch reactor component species balances used to model series reaction.
    
    CA = c_ini[0]
    CB = c_ini[1]
    CC = c_ini[2] 
    
    k1  = k[0] 
    k2  = k[1] 
    
    r1 = k1*CA 
    r2 = k2*CB 
    
    dCAdt = -r1
    dCBdt = r1 -r2
    dCCdt = r2
    
    return [dCAdt,dCBdt,dCCdt]

#For reaction A -> B :
k1  = 2     #[1/s]

# For reaction B -> C :
k2  = 1     #[1/s]

k = [k1, k2]

# Time length for solution:
t = np.linspace(0,5,110)

# Initial concentrations of A,B and C:
CA0 = 1      #[mol/L]
CB0 = 0      #[mol/L]
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
plt.legend(loc='center right',shadow='True',fontsize=14)
plt.show()

#----------------------------------------------------------------------------------------------------
#           Part 2 : Solving Full Reaction Model ODE Equations With Increasing k2 values
#----------------------------------------------------------------------------------------------------

# Case 1:
k2_10  = 10     #[1/s]
k_10 = [k1, k2_10]

conc_sol_10 = solve_ivp(series_reactions_batch, [0, t[-1]], c_ini, t_eval=t, method='Radau', args=[k_10])
CA_10=conc_sol_10.y[0]
CB_10=conc_sol_10.y[1]
CC_10=conc_sol_10.y[2]

# Case 2:
k2_100  = 100     #[1/s]
k_100 = [k1,k2_100]

conc_sol_100 = solve_ivp(series_reactions_batch, [0, t[-1]], c_ini, t_eval=t, method='Radau', args=[k_100])
CA_100=conc_sol_100.y[0]
CB_100=conc_sol_100.y[1]
CC_100=conc_sol_100.y[2]

# Case 3:
k2_1000  = 1000     #[1/s]
k_1000 = [k1, k2_1000]

conc_sol_1000 = solve_ivp(series_reactions_batch, [0, t[-1]], c_ini, t_eval=t, method='Radau', args=[k_1000])
CA_1000=conc_sol_1000.y[0]
CB_1000=conc_sol_1000.y[1]
CC_1000=conc_sol_1000.y[2]

# Case 4:
k2_10000  = 10000     #[1/s]
k_10000 = [k1, k2_10000]

conc_sol_10000 = solve_ivp(series_reactions_batch, [0, t[-1]], c_ini, t_eval=t, method='Radau', args=[k_10000])
CA_10000=conc_sol_10000.y[0]
CB_10000=conc_sol_10000.y[1]
CC_10000=conc_sol_10000.y[2]


# Plotting Results:
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(20, 15), dpi=150) 
fig.suptitle('Concentration Evolution of Species in Batch Reactor With Increasing k$_2$ ',fontsize = 37,y=0.95)
axs[0,0].plot(t,CA,'k-',label='A| k$_2$ = 1')
axs[0,0].plot(t,CB,'r-',label='B| k$_2$ = 1')
axs[0,0].plot(t,CC,'g-',label='C| k$_2$ = 1')
axs[0,0].plot(t,CA_10,'k--',label='A| k$_2$ = 10')
axs[0,0].plot(t,CB_10,'r--',label='B| k$_2$ = 10')
axs[0,0].plot(t,CC_10,'g--',label='C| k$_2$ = 10')
axs[0,0].set_xlabel('Time [s]',fontsize = 20)
axs[0,0].set_ylabel('Concentration [mol/L]',fontsize = 20)
axs[0,0].set_title('k$_2$ = 10 [1/s]',fontsize = 20)
axs[0,0].grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
axs[0,0].legend(loc='upper right',shadow='True',fontsize=14)  

axs[0,1].plot(t,CA,'k-',label='A| k$_2$ = 1')
axs[0,1].plot(t,CB,'r-',label='B| k$_2$ = 1')
axs[0,1].plot(t,CC,'g-',label='C| k$_2$ = 1')
axs[0,1].plot(t,CA_100,'k--',label='A| k$_2$ = 100')
axs[0,1].plot(t,CB_100,'r--',label='B| k$_2$ = 100')
axs[0,1].plot(t,CC_100,'g--',label='C| k$_2$ = 100')
axs[0,1].set_xlabel('Time [s]',fontsize = 20)
axs[0,1].set_ylabel('Concentration [mol/L]',fontsize = 20)
axs[0,1].set_title('k$_2$ = 100 [1/s]',fontsize = 20)
axs[0,1].grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
axs[0,1].legend(loc='upper right',shadow='True',fontsize=14)  

axs[1,0].plot(t,CA,'k-',label='A| k$_2$ = 1')
axs[1,0].plot(t,CB,'r-',label='B| k$_2$ = 1')
axs[1,0].plot(t,CC,'g-',label='C| k$_2$ = 1')
axs[1,0].plot(t,CA_1000,'k--',label='A| k$_2$ = 1000')
axs[1,0].plot(t,CB_1000,'r--',label='B| k$_2$ = 1000')
axs[1,0].plot(t,CC_1000,'g--',label='C| k$_2$ = 1000')
axs[1,0].set_xlabel('Time [s]',fontsize = 20)
axs[1,0].set_ylabel('Concentration [mol/L]',fontsize = 20)
axs[1,0].set_title('k$_2$ = 1000 [1/s]',fontsize = 20)
axs[1,0].grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
axs[1,0].legend(loc='upper right',shadow='True',fontsize=14)   

axs[1,1].plot(t,CA,'k-',label='A| k$_2$ = 1')
axs[1,1].plot(t,CB,'r-',label='B| k$_2$ = 1')
axs[1,1].plot(t,CC,'g-',label='C| k$_2$ = 1')
axs[1,1].plot(t,CA_10000,'k--',label='A| k$_2$ = 10000')
axs[1,1].plot(t,CB_10000,'r--',label='B| k$_2$ = 10000')
axs[1,1].plot(t,CC_10000,'g--',label='C| k$_2$ = 10000')
axs[1,1].set_xlabel('Time [s]',fontsize = 20)
axs[1,1].set_ylabel('Concentration [mol/L]',fontsize = 20)
axs[1,1].set_title('k$_2$ = 10000 [1/s]',fontsize = 20)
axs[1,1].grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
axs[1,1].legend(loc='upper right',shadow='True',fontsize=14)   
plt.show()




