import numpy as np
import matplotlib.pyplot as plt
import lhapdf
import csv
import vegas

#PDFs
#p = lhapdf.mkPDF("MRST2001lo", 0)
p = lhapdf.mkPDF("cteq6l1", 0)


def ggToccbar(sm,mc,alphaS):
    w=(4.0*(mc**2))/sm
    if 0<w<1:
        a0 = (np.pi*alphaS**2)/(3.0*sm)
        a1 = (1.0 + w + (w**2.0)/16.0)
        a2 = np.log( (1.0+np.sqrt(1.0-w)) / (1.0-np.sqrt(1.0-w)) )
        a3 = ((7.0/4.0) + (31.0*w)/16.0)*np.sqrt(1.0-w)
        return a0*(a1*a2 - a3)
    else:
        return 0.0

def qqbarToccbar(sm,mc,alphaS):
     w=(4.0*mc**2)/sm
     return (8.0*np.pi*alphaS**2)/(27.0*sm)*(1 + w/2.0)*np.sqrt(1.0-w)

def WWZ(y,sqs,sm):
    x2=np.sqrt(sm)*np.exp(-y)*(1/sqs)
    em=0.511*10**(-3)
    E=0.5*sqs
    if 0<x2<1:
        b0=((1+(1-x2)**2)/x2)*(np.log(E/em)-0.5)
        b1=(x2/2)*(np.log((2/x2)-2)+1)
        b2=(((2-x2)**2)/2*x2)*np.log((2-2*x2)/(2-x2))
        return (alpha/np.pi)*(b0+b1+b2)
    else:
        return 0.0


def AKPERP(qT1,kt):
    al=0.25
    bt=0.25
    TMDs= np.exp(-kt**2/al)*np.exp()
    return (np.exp((-qT1**2)/(al+bt)))/(np.pi*(al+bt))

def ppToC(x, sqs, qT1):

     #x1=x[0]
     y=x[0]
     sm=x[1]
     #q1T=x[3]
     phiq1T=x[2]
     #phikT=x[3]
     x1=np.sqrt(sm)*np.exp(y)*(1/sqs)  
     s=sqs*sqs

    #PDFs
     glx1 = p.xfxQ2(21,x1,sm)/x1   
     alphaS=p.alphasQ2(sm)
     return  (sm*qT1*glx1*AKPERP(qT1)*WWZ(y,sqs,sm)*ggToccbar(sm,mc,alphaS))/s
 
mc=1.47
mD=1.87
mC=3.096 
alpha=1/137
yLim=[0.0,1.0]
smLim=[4*mc**2, 4*mD**2]
#q1TLim=[0.0,5.0]
phiq1TLim=[0.0,2.0*np.pi]
# Integration limits for each dimension
limits = [(yLim),(phiq1TLim)]
Bmumu=0.0601
conv=3.899E5*Bmumu
data=[]
#Veribles
sqs=31.6
qT1=0.025

def f(x):
    return ppToC(x,sqs,qT1)

print("Adapting integrand")
print("Actual vegas run started:")
while qT1<1:
     integ = vegas.Integrator([yLim, smLim, phiq1TLim], nproc=4)
     result = integ(f, nitn=10, neval=50000)
     print(qT1, result.mean)
     data.append([qT1, conv*result.mean])
     with open (f"data/denominator1.dat",'w') as file:
        writer = csv.writer(file,delimiter = "\t")
        writer.writerows(data)
     qT1+=0.025

#Plotting:
plt.title('Denominator')
plt.xlim([0, 5])
plt.yscale("log")
#plt.ylim([1E-3, 2E1])
plt.xlabel(r"$p_{T}~[GeV]$")
plt.ylabel(r"$Br[J/\psi \rightarrow \mu^{+}\mu^{-}]\times d\sigma/dp_{T\psi}~[nb/GeV^{2}]$")

#Karpishkov = np.loadtxt('data/Karpishkov2020brv.dat')
#plt.step(Karpishkov[:,0], Karpishkov[:,1], where='mid', label='Karpishkov:2020brv')
dataplot=np.loadtxt(f"data/denominator1.dat")
plt.plot(dataplot[:,0], dataplot[:,1],  label=f'eptoJ/psi')

plt.legend()
plt.savefig("epCEM_d1.pdf", format="pdf", dpi=300)
plt.show()