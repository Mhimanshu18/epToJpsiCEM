#!/home/himanshu/anaconda3/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import lhapdf
import csv
import vegas

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

def WWZ(xx, sqs, M2):
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

def Sivers_xpart(x1):
    Anu=0.4
    And=-0.97
    Alu=0.35
    Ald=0.44
    Abu=2.6
    Abd=0.90
    NuX= (Anu*x1**Alu)*((1-x1)**Abu)*(((Alu+Abu)**(Alu+Abu))/(Alu**Alu*Abu**Abu))
    NdX= (And*x1**Ald)*((1-x1)**Abd)*(((Ald+Abd)**(Ald+Abd))/(Ald**Ald*Abd**Abd))
    Ng1x=(NuX+NdX)/2
    Ng2x=NdX
    return Ng1x

def TMDpart(set, qT):
    al=0.25
    bt=0.25
    
    Am1=np.sqrt(0.19)
    aks=(al*Am1*Am1)/((Am1*Am1)+al)
    
    final = 0
    if(set=="D"):
       final = (np.exp((-qT**2)/(al+bt)))/(np.pi*(al+bt))
    elif (set=="N"):
        final = ((qT)/Am1)*((aks**2))*np.sqrt(2.*np.exp(1))\
        *np.exp((-qT**2)/(aks+bt))/(al*np.pi*(aks+bt)**2)
     
    return final


#Integrand with parameters
def epToJpsi_CEM(x, flag, sqs, mc, y):
        qT=x[0]
        phi_qT= x[1] 
        M2 = x[2]

        x_g=np.sqrt(M2)*np.exp(y)*(1/sqs)  
        x_p=np.sqrt(M2)*np.exp(-y)*(1/sqs)  

        if (0 < x_g < 1):

            s=sqs*sqs
            glx_g = p.xfxQ2(21,x_g,M2)/x_g   
            alphaS = p.alphasQ2(M2)
            energy = sqs/2

            if  (flag=="D"):
                return (2.0/s)*qT*glx_g*fgammae(x_p, energy)*TMDpart("D", qT)*ggToccbar(M2,mc,alphaS)
            elif  (flag=="N"):
                return (1.0/s)*qT*glx_g*2*Sivers_xpart(x_g)*fgammae(x_p, energy)*TMDpart("N", qT)*ggToccbar(M2,mc,alphaS)*np.sin((np.pi/2.0)-phi_qT)*np.sin(phi_qT-(np.pi/2.0))
            else:
                print('Invalid flag!')
    
#main
sqs = 7.2
mc = 1.27
y = -0.8
mD = 1.87

def denominator(x):
    return epToJpsi_CEM(x, "D", sqs, mc, y)

def numerator(x):
    return epToJpsi_CEM(x, "N", sqs, mc, y)


qT_limits = [0.0, 1.0]
phi_qT_limits = [0.0, 2.0*np.pi]
M2_limits = [4.0*mc**2, 4.0*mD**2]

data=[]
while (y < 0.6):
    integrator = vegas.Integrator([qT_limits, phi_qT_limits, M2_limits], nproc=8)
    result_den = integrator(denominator, nitn=10, neval=20000)
    result_num = integrator(numerator, nitn=10, neval=20000)
    print(y, result_den.mean,result_num.mean)
    data.append([y, result_den.mean,result_num.mean])
    with open (f"data/gmmr_2012.dat",'w') as file:
        writer = csv.writer(file,delimiter = "\t")
        writer.writerows(data)

    y+=0.1

#Plotting #
x = [x[0] for x in data]  # Extract the y values
y_d = [x[1] for x in data]  # Extract the denominator values
y_n = [x[2] for x in data]  # Extract the numerator values

# Calculate the ratio of numerator to denominator
y_ratio = [n / d for n, d in zip(y_n, y_d)]

# Plot the ratio
plt.plot(x, y_ratio)
plt.xlabel('y')
plt.ylabel('Numerator / Denominator')
plt.title('Numerator / Denominator vs y')
plt.grid(True)
plt.show()




