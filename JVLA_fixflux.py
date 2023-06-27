#!/use/bin/python3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def get_I(caldatfile):
    data = np.loadtxt(caldatfile)
    def S(f,S,alpha,beta):
        return S*(f/3.0)**(alpha+beta*np.log10(f/3.0))

    # Fit 1 - 5 GHz data points
    popt, pcov = curve_fit(S, data[0:6,0], data[0:6,1])
    print('I@3GHz', popt[0], ' Jy')
    print('alpha', popt[1])
    print('beta', popt[2])
    print( 'Covariance')
    print(pcov)
    
    plt.plot(data[0:6,0], data[0:6,1], 'ro', label='data')
    plt.plot(np.arange(1,5,0.1), S(np.arange(1,5,0.1), *popt), 'r-', label='fit')
    
    plt.title('3C286')
    plt.legend()
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Flux Density (Jy)')
    plt.show()
    return(popt)
    
    
def get_PF(caldatfile):
    data = np.loadtxt(caldatfile)
    def PF(f,a,b,c,d):
            return a+b*((f-3.0)/3.0)+c*((f-3.0)/3.0)**2+d*((f-3.0)/3.0)**3
    
    # Fit 1 - 5 GHz data points
    popt, pcov = curve_fit(PF, data[0:6,0], data[0:6,2])
    print("Polfrac Polynomial: ", popt)
    print("Covariance")
    print(pcov)
    
    plt.plot(data[0:6,0], data[0:6,2], 'ro', label='data')
    plt.plot(np.arange(1,5,0.1), PF(np.arange(1,5,0.1), *popt), 'r-', label='fit')
    
    plt.title('3C48')
    plt.legend()
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Lin. Pol. Fraction')
    plt.show()
    return(popt)


def get_PA(caldatfile):
    data = np.loadtxt(caldatfile)
    def PA(f,a,b,c,d,e):
        return a+b*((f-3.0)/3.0)+c*((f-3.0)/3.0)**2+d*((f-3.0)/3.0)**3+e**((f-3.0)/3.0)**4

    # Fit 1 - 5 GHz data points
    popt, pcov = curve_fit(PA, data[0:6,0], data[0:6,3])
    print("Polangle Polynomial: ", popt)
    print("Covariance")
    print(pcov)
    
    plt.plot(data[2:8,0], data[2:8,3], 'ro', label='data')
    plt.plot(np.arange(1,9,0.1), PA(np.arange(1,9,0.1), *popt), 'r-', label='fit')
    
    plt.title('3C48')
    plt.legend()
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Lin. Pol. Angle (rad)')
    plt.show()
    return popt
