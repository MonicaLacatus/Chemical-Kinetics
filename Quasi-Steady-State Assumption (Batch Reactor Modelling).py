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
    k_1 = k[1] 
    k2  = k[2] 
    k_2 = k[3]
    
    r1 = k1*CA - k_1*CB
    r2 = k2*CB - k_2*CC
    
    dCAdt = -r1
    dCBdt = r1 -r2
    dCCdt = r2
    
    return [dCAdt,dCBdt,dCCdt]




