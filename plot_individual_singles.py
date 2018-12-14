## Plot the eccentricity preference for each single planet system
##

import numpy as np
import matplotlib.pyplot as plt
import random
import scipy
from scipy import stats
import os
from scipy import interpolate
from cksecc_helper import *


cks_data, koi_data, koi_errs = load_data()
koilist, propertymatrix = divide_into_singles_and_multis(cks_data, koi_data, koi_errs, singles=True, more_properties=True) 

e, inc, like, emaxval = read_in_singles() 

propertymatrix.append(like)
koilist, propertymatrix = clean_sample(koilist, propertymatrix, customremove=True)


msun = propertymatrix[0]
rsun = propertymatrix[1]
period = propertymatrix[2]
radius = propertymatrix[3]
teff = propertymatrix[4]
metal = propertymatrix[5]
age = propertymatrix[6]
dilution = propertymatrix[7]
durations = propertymatrix[8]
rpors = propertymatrix[9]
like = propertymatrix[-1]


alltotal = np.sum(np.log(like), axis=0)
allmaxllike = max(alltotal)
allllike = alltotal-allmaxllike



def get_aors(period, rho_star):
  # [universal gravitational constant]^(1/3) * (day)^(2/3) * [solar density]^(1/3)* ( 1/(3*pi) )^(1/3)
  aorsfactor = 4.206
  aors = aorsfactor * period**(2./3.) * (rho_star*(4./3.*np.pi))**(1./3.)
  if (aors < 1):
    aors = np.nan
  return aors
def compute_circular_edgeon_duration(period, rpors, aors):
  prefactor = 1./np.pi
  T0 = prefactor * period * np.arcsin( (1.+rpors)/aors ) 
  return T0


######
kois=koilist

print "Nsystems=", len(kois)

print "koi loglike dur dur_e dur_circ Rearth"

preflist = []
d0list = []
dlist = []
radlist = []

for i in range(len(like)): 

  l = like[i]
  ll = np.log(l)
  ll = ll-max(ll)
  
  sigma = np.sqrt(-ll*2)

  #print koi[i]
  plt.figure(1)
  plt.semilogy(e, sigma)#, color='k') 
  plt.figure(2)
  plt.plot(e, sigma)#, color='k') 
  #plt.plot(e, ll, color='k') 

  #if kois[i] % 1 != 0.01:
  #  print kois[i]

  #if min(ll) < -2:
  if np.sqrt(-min(ll)*2) > 3:
    rho_stari = msun[i] / (4./3.*np.pi*rsun[i]**3.)
    aorsi = get_aors(period[i], rho_stari)
    print kois[i], np.sqrt(-min(ll)*2), durations[i], compute_circular_edgeon_duration(period[i], rpors[i], aorsi), radius[i]
 
  preflist.append(ll[-1]-ll[0])
  d0list.append(durations[i])
  radlist.append(radius[i])
  rho_stari = msun[i] / (4./3.*np.pi*rsun[i]**3.)
  aorsi = get_aors(period[i], rho_stari)
  dlist.append(compute_circular_edgeon_duration(period[i], rpors[i], aorsi))

prefsort = list(reversed(np.argsort(preflist)))
np.savetxt('koipreflist.txt',  np.transpose([np.array(kois)[prefsort], np.array(preflist)[prefsort],
    np.array(d0list)[prefsort], np.array(dlist)[prefsort], np.array(radlist)[prefsort]]), fmt=['%.2f', '%15.7f', '%15.7f', '%15.7f', '%15.7f'])


plt.figure(1)
plt.ylabel('Disfavored at $\sigma$ Level', fontsize=12)
plt.xlabel('$\sigma_e$', fontsize=12)
plt.xlim((0.,0.7))
plt.ylim((0.1, 8.))

plt.savefig('figs4/individual_singles.png')

plt.close()



plt.figure(2)
plt.ylabel('Disfavored at $\sigma$ Level', fontsize=12)
plt.xlabel('$\sigma_e$', fontsize=12)
plt.xlim((0.,0.7))
plt.ylim((0., 40.))

plt.savefig('figs4/individual_singles2.png')








