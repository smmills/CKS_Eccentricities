## Plot the eccentricity preference for each planet in every multiplanet system
##

import numpy as np
import matplotlib.pyplot as plt
import random
import scipy
from scipy import stats
import os
from scipy import interpolate



from s_g_funs import * 
from cksecc_helper import *


cks_data, koi_data, koi_errs = load_data()
koilist, propertymatrix = divide_into_singles_and_multis(cks_data, koi_data, koi_errs, singles=False, more_properties=True) 

e, inc, like = read_in_multis(alle=True)

print len(koilist)
print len(propertymatrix)
print len(propertymatrix[0])


propertymatrix.append(like)
koilist, propertymatrix = clean_sample(koilist, propertymatrix, customremove=True, singles=False)

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

point5=[2,3]
ncuts=5
ncuts2=2
nincut = 10
print e
print inc

e = e.reshape(ncuts+len(point5), nincut*ncuts2) 
inc = inc.reshape(ncuts+len(point5), nincut*ncuts2) 
like = like.reshape(len(like), ncuts+len(point5), nincut*ncuts2) 

print e
print inc
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



incdex=6
print "inc:", inc[incdex]
like = like[:,incdex,:]
print "like shape"
print like.shape

alltotal = np.sum(np.log(like), axis=0)
allmaxllike = max(alltotal)
allllike = alltotal-allmaxllike



print e[incdex]
kois=koilist

print "Nsystems=", len(kois)


for i in range(len(like)): 

  l = like[i]
  ll = np.log(l)
  ll = ll-max(ll)
  
  sigma = np.sqrt(-ll*2)

  #print koi[i]
  plt.figure(1)
  plt.semilogy(e[incdex], sigma)#, color='k') 
  plt.figure(2)
  plt.plot(e[incdex], sigma)#, color='k') 
  #plt.plot(e, sigma, color='k') 
  #plt.plot(e, ll, color='k') 

  #if kois[i] % 1 != 0.01:
  #  print kois[i]

  #if min(ll) < -2:
  if np.sqrt(-min(ll)*2) > 3:
    rho_stari = msun[i] / (4./3.*np.pi*rsun[i]**3.)
    aorsi = get_aors(period[i], rho_stari)
    llmaxi = np.argmax(ll)
    emax = e[incdex][llmaxi]
    print kois[i], np.sqrt(-min(ll)*2), durations[i], compute_circular_edgeon_duration(period[i], rpors[i], aorsi), radius[i], emax 
  if np.sqrt(-min(ll)*2) > 3.8:
    plt.figure(3)
    plt.plot(e[incdex], sigma)
    plt.xlabel('$\sigma_e$', fontsize=12)
    plt.ylabel('Disfavored at $\sigma$ Level', fontsize=12)
    plt.savefig('figs4/Kepler-10c.png')

plt.ylabel('Disfavored at $\sigma$ Level', fontsize=12)
plt.xlabel('$\sigma_e$', fontsize=12)
plt.xlim((0.,0.2))
plt.ylim((0.1, 4.))

plt.savefig('figs4/individual_multis.png')

#plt.show()
#plt.close()

  #raw_input("Press Enter to continue...")


plt.figure(2)
plt.ylabel('Disfavored at $\sigma$ Level', fontsize=12)
plt.xlabel('$\sigma_e$', fontsize=12)
plt.xlim((0.,0.20))
plt.ylim((0., 4.))

plt.savefig('figs4/individual_multis2.png')





