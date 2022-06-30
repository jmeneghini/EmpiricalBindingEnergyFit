import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import elementary_charge
from scipy.optimize import curve_fit
import os 

plt.style.use(['science', 'grid']) 

data_location = os.path.dirname(os.path.realpath(__file__)) + '/OGBindingEnergy.csv'

data = pd.read_csv(data_location) # read in the data

data = data[data['Binding Energy'] > 0] # remove negative values
data['Nucleon Number'] = data['Nucleon Number'].astype(int) # convert to int
data['Proton Number'] = data['Proton Number'].astype(int) # convert to int


def semiEmpiricalBindingEnergyFit(X, a_v, a_s, a_c, a_sym, a_p): 
    A, Z = X

    delta = np.array([1 if (A[i] % 2 == 0) and (Z[i] % 2 == 0) else -1 if (A[i] % 2 != 0) and (Z[i] % 2 != 0) else 0 for i in range(len(A))]) # make array of +1 or -1 or 0 for even or odd A and Z

    binding_energy = a_v*A - a_s*A**(2/3) - a_c*Z*(Z-1)*A**(-1/3) - a_sym*(A-2*Z)**2/A + delta*a_p*A**(-3/4) # calculate binding energy
    
    return binding_energy

popt, pcov = curve_fit(semiEmpiricalBindingEnergyFit, (data["Nucleon Number"], data["Binding Energy"]), data['Binding Energy']) # fit the data
print(pcov) # print the covariance matrix
popt = popt*1e3 # convert to meV
print(f"a_v = {popt[0]}\na_s = {popt[1]}\na_c = {popt[2]}\na_sym = {popt[3]}\na_p = {popt[4]}") # print the fit parameters

plt.plot(data['Nucleon Number'], data['Binding Energy']*1e3/data['Nucleon Number'], 'o', markersize = 2, label = "Data") # plot the data
plt.plot(data['Nucleon Number'], semiEmpiricalBindingEnergyFit((data["Nucleon Number"].to_numpy(), data["Binding Energy"].to_numpy()),
                                                               popt[0], popt[1], popt[2], popt[3], popt[4])/data["Nucleon Number"], 'o', markersize = 2, label = "Fit") # plot the fit data
plt.xlabel('Nucleon Number') 
plt.ylabel('Binding Energy (MeV) per Nucleon')
plt.legend()
plt.show()

print(data["Binding Energy"])