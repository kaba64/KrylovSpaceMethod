from matplotlib import rc
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt

#fldnum = 
datadir    = ''
dataname_c = 'C_residual'
dataname_g = 'G_residual'

# Lodad the data

fname_c   = datadir+dataname_c+'.txt'
fname_g   = datadir+dataname_g+'.txt'

data_c    = np.loadtxt(fname_c)
data_g    = np.loadtxt(fname_g)

# Find the size of each data
#plt.figure(figsize=(8.5, 6))

plt.plot(data_c[:,0],data_c[:,1],'-',color='red')
plt.plot(data_g[:,0],data_g[:,1],'-',color='green')

#plt.legend((r'CGM : $\epsilon$',r'GSM $\epsilon$'))

#plt.xlabel(r"$iteration$",fontsize=18)
#plt.ylabel(r"${\frac{r}{r{0}}$",fontsize=18)

plt.xlim((0,200))
plt.ylim((0,0.5))
#plt.xticks(fontsize=18)
#plt.yticks(fontsize=18)

plt.show()
#plt.savefig('nost_l-oh-w.png')
#plt.savefig('nost_l-oh-e1.png')
#plt.savefig('nost_l-oh-e2.png')
