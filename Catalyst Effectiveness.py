#---------------------------------------------------------------------------------------------------------------------
#  ______  __  __          _   _                              ____   __    _____      _        _           _       
# |  ____|/ _|/ _|        | | (_)                            / __ \ / _|  / ____|    | |      | |         | |      
# | |__  | |_| |_ ___  ___| |_ ___   ___ __   ___  ___ ___  | |  | | |_  | |     __ _| |_ __ _| |_   _ ___| |_ ___ 
# |  __| |  _|  _/ _ \/ __| __| \ \ / / '_ \ / _ \/ __/ __| | |  | |  _| | |    / _` | __/ _` | | | | / __| __/ __|
# | |____| | | ||  __/ (__| |_| |\ V /| | | |  __/\__ \__ \ | |__| | |   | |___| (_| | || (_| | | |_| \__ \ |_\__ \
# |______|_| |_| \___|\___|\__|_| \_/ |_| |_|\___||___/___/  \____/|_|    \_____\__,_|\__\__,_|_|\__, |___/\__|___/
#                                                                                                 __/ |            
#                                                                                                 |___/ 
#---------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------
#                                           Importing Libraries
#---------------------------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["axes.labelsize"] = 14

#---------------------------------------------------------------------------------------------------------------------
#               Infleunce of Thiele Modulus on Concentration Profile Through Slab
#---------------------------------------------------------------------------------------------------------------------

thiele = np.array([0.1, 1, 2, 10])
x = np.linspace(1,0,101)
C = np.zeros((len(x),4))

for i in range(len(thiele)):
    C[:,i] = np.cosh(thiele[i]*x)/np.cosh(thiele[i])
    
    
fig = plt.figure(figsize= (9,6))
plt.gca().invert_xaxis()
plt.plot(x,C[:,0],label='Φ=0.1')
plt.plot(x,C[:,1],label='Φ=1')
plt.plot(x,C[:,2],label='Φ=2')
plt.plot(x,C[:,3],label='Φ=10')
plt.grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
plt.xlabel('x$^*$ [-]')
plt.ylabel('C$^*$ [-]')
plt.title('Diffusion and Reaction Through Slab', fontsize = 16)
plt.legend(loc='lower left',shadow='True',fontsize=12)
plt.show()

#---------------------------------------------------------------------------------------------------------------------
#         Influence of Thiele Modulus on Effectiveness Factor of Different Catalyst Shapes
#---------------------------------------------------------------------------------------------------------------------

thiele = np.linspace(0,100,10001)
ni_slab = np.tanh(thiele)/thiele
ni_sphere = (1/thiele)*(1/np.tanh(3*thiele)-1/(3*thiele))
ni_cylinder = (1/thiele)*(sp.iv(1,2*thiele)/sp.iv(0,2*thiele))

fig = plt.figure(figsize= (9,6))
plt.loglog(thiele,ni_slab,label = 'slab')
plt.loglog(thiele,ni_cylinder, label = 'cylinder')
plt.loglog(thiele,ni_sphere, label = 'sphere')
plt.xlabel('Φ [-]')
plt.ylabel('η$_i$ [-]')
plt.title('Effectivness Factor for Internal Mass Transfer - Different Shape Catalysts', fontsize = 14)
plt.grid(True, color = "grey", linewidth = "1.0", linestyle = "-")
plt.legend(loc='lower left',shadow='True',fontsize=12)
plt.show()



