#!/home/himanshu/anaconda3/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import lhapdf
import csv
import vegas
import pandas as pd
p = lhapdf.mkPDF("MRST2001lo", 0)

def ggToccbar(M2,mc,alphaS):
    w=(4.0*(mc**2))/M2
    alpha = 1.0/137.0
    if 0<w<1:
        a0 = 0.5*(2/3)**2*(4*np.pi*alpha*alphaS)
        a1 = (1.0 + w - 0.5*(w**2.0))
        a2 = np.log( (1.0+np.sqrt(1.0-w))/(1.0-np.sqrt(1.0-w)) )
        a3 = ((1+w)*np.sqrt(1.0-w))
        return a0*(a1*a2 - a3)/M2
    else:
        return 0.0

def WWZ(xx, sqs):
    em=0.511*10**(-3)
    E=0.5*sqs
    alpha = 1.0/137.0
    if 0<xx<1:
        b0=((1+(1-xx)**2)/xx)*(np.log(E/em)-0.5)
        b1=(xx/2)*(np.log((2/xx)-2)+1)
        b2=(((2-xx)**2)/2*xx)*np.log((2-2*xx)/(2-xx))
        return (alpha/np.pi)*(b0+b1+b2)
    else:
        return 0.0
    
def fgammae(xx, ee):
    alpha= 1.0/ 137.0
    me = 0.511e-3
    WWZ = 0.0
    if (0 < xx <1):
        WWZ = ((alpha / np.pi) * (
                ((1.0 + ((1.0 - xx) * (1.0 - xx))) / xx) * (np.log(ee / me) - 0.5)\
                + (xx * 0.50) * (np.log((2.0 / xx) - 2.0) + 1.0)\
                + (((2.0 - xx) * (2.0 - xx)) / (2.0 * xx)) * (np.log((2.0 - (2.0 * xx)) / (2.0 - xx)))
        ))
    return WWZ

def TMDpart(qT):
    al=1.23
    bt=1.23
    return (np.exp((-qT**2)/(al+bt)))/(np.pi*(al+bt))
 

#Integrand with parameters
def epToJpsi_CEM(x, sqs, mc, qT2):
        z=x[0]
        phi_qT= x[1] 
        M2 = x[2]
        MJpsi=3.096
        MT=np.sqrt(mc**2+qT2)
        x_g=(MT/sqs)*z
        x_p=(MT/sqs)  
        qT=np.sqrt(qT2)
        FJpsi=0.6
        if (0 < x_g < 1):

            s=sqs*sqs
            glx_g = p.xfxQ2(21,x_g,M2)/x_g   
            alphaS = p.alphasQ2(M2)
            energy = sqs/2

            return FJpsi*(1.0/(2.0*z*s))*glx_g*fgammae(x_p, energy)*TMDpart(qT)*ggToccbar(M2,mc,alphaS)

#main
sqs = 319
mc = 1.27
qT2 = 1.0
mD = 1.87

def denominator(x):
    return epToJpsi_CEM(x, sqs, mc, qT2)

z_limits = [0.3, 0.9]
phi_qT_limits = [0.0, 2.0*np.pi]
M2_limits = [4.0*mc**2, 4.0*mD**2]
Conv=3.89E5
data=[]
while (qT2 < 100):

    integrator = vegas.Integrator([z_limits, phi_qT_limits, M2_limits], nproc=8)
    result_den = integrator(denominator, nitn=10, neval=20000)
    print(qT2, result_den.mean*Conv)#,result_num.mean)

    data.append([qT2, result_den.mean])
    with open (f"data/gmmpr_2024.dat",'w') as file:
        writer = csv.writer(file,delimiter = "\t")
        writer.writerows(data)

    qT2+=10

#Plotting #

#HERA DATA
h1_data = pd.read_csv("data/Table1.csv", comment="#")
cem = np.loadtxt('data/gmmpr_2024.dat', delimiter='\t')
pT_low = h1_data['PT(P=3) [GEV**2] LOW']  # Lower bound of pT range
pT_high = h1_data['PT(P=3) [GEV**2] HIGH']  # Upper bound of pT range
dSigma_dPT2 = h1_data['D(SIG)/DPT [NB/GEV**2]']  # Differential cross-section
stat_plus = h1_data['stat +']  # Statistical uncertainty (positive)
stat_minus = h1_data['stat -']  # Statistical uncertainty (negative)
cem_x = cem[:, 0]
cem_y = cem[:, 1]
# Calculate the x-error (pT range uncertainty)
xerr_lower = (pT_high - (pT_low + pT_high) / 2)  # Error to the right of the central value
xerr_upper = ((pT_low + pT_high) / 2 - pT_low)  # Error to the left of the central value

# Create the plot with error bars
plt.errorbar((pT_low + pT_high) / 2, dSigma_dPT2, 
             xerr=[xerr_lower, xerr_upper],  # X-errors
             yerr=[abs(stat_minus), abs(stat_plus)],  # Y-errors (statistical uncertainties)
             fmt='o', label='H1', capsize=5)

plt.plot(cem_x, cem_y, label='cem_158.1', color='blue', marker='o', linestyle='-')


plt.xlabel(r'$p_T^2$ [GeV$^2$]')
plt.ylabel(r'$\frac{d\sigma}{dp_T^2}$ [nb/GeV$^2$]')

plt.yscale('log')
plt.xscale('log')
plt.title('HERA')
plt.legend()
plt.show()




