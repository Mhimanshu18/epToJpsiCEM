#Plotting:
import matplotlib.pyplot as plt
import numpy as np
plt.title('Asymmetry vs pt')
plt.xlim([0, 0.8])
#plt.yscale("log")
plt.ylim(0, 0.25)
plt.xlabel(r"$p_{T}~[GeV]$")
plt.ylabel(r"$Br[J/\psi \rightarrow \mu^{+}\mu^{-}]\times d\sigma/dp_{T\psi}~[nb/GeV^{2}]$")
num=np.loadtxt(f"data/numerator1.dat")
den=np.loadtxt(f"data/denominator1.dat")
plt.plot(num[:,0], num[:,1]/(2*den[:,1]),  label=f'eptoJ/psi')

plt.legend()
plt.savefig("epASYM_p1.pdf", format="pdf", dpi=300)
plt.show()