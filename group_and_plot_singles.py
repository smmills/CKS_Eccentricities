## Plot the overall eccentricity fit for all of the single planets
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
koilist, propertymatrix = divide_into_singles_and_multis(cks_data, koi_data, koi_errs, singles=True) 

e, inc, like, emaxval = read_in_singles() 

print "e"
print e
print "inc"
print inc

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
like = propertymatrix[-1]


alltotal = np.sum(np.log(like), axis=0)
allmaxllike = max(alltotal)
allllike = alltotal-allmaxllike

print "overall loglike"
print allmaxllike 

from s_g_funs import quadratic_intersections, savitzky_golay 


def makeplot(propertylist, cutlist, propertystring):

  nbins = len(cutlist)+1
  first = 0
  last = nbins-1

  koigroups = []

  for i in range(nbins):
    if i==first:
      koigroups.append(np.where(propertylist < cutlist[i]))
    elif i==last:
      koigroups.append(np.where(propertylist > cutlist[i-1]))
    else:
      koigroups.append(np.where(np.logical_and(cutlist[i-1] < propertylist, propertylist < cutlist[i])))

  likegroups = [like[koigroup] for koigroup in koigroups]
  liketotals = [np.sum(np.log(likegroup), axis=0) for likegroup in likegroups]
  normlikes = [liketotal - max(liketotal) for liketotal in liketotals]

  strings = []
  for i in range(nbins):
    if i==first:
      strings.append('%s < %s (N=%i)' % (propertystring, str(cutlist[i]), len(likegroups[i]))) 
    elif i==last:
      strings.append('%s > %s (N=%i)' % (propertystring, str(cutlist[i-1]), len(likegroups[i]))) 
    else:
      strings.append('%s < %s < %s (N=%i)' % (str(cutlist[i-1]), propertystring, str(cutlist[i]), len(likegroups[i]))) 

  colors = ['r', 'b', 'g', 'y', 'c', 'm', 'tan', 'maroon', 'mediumspringgreen', 'cornflowerblue', 'orange', 'gold']

  maxes = []
  maxeerrs1 = []
  maxeerrs2 = []

  #plt1, ax1 = plt.plot(e, allllike, color='k', label='All') 
  plt.figure(1)
  plt.plot(e, allllike, color='k', label='All') 

  nlall = allllike
  allfit = savitzky_golay(nlall, 7, 3)
  plt.plot(e, allfit, color='k', ls='dashed')

  for i in range(nbins):
    plt.figure(1)
    nl = normlikes[i]
    plt.plot(e, nl, color=colors[i % len(colors)], label=strings[i]) 
    nl10i = np.where(nl > -10.)
    e10 = e[nl10i]
    nl10 = nl[nl10i]
    poly = np.polyfit(e10, nl10, 2)
    quadfit = np.poly1d(poly)
    xfit = np.linspace(0., 0.4, 50)
    #plt.plot(xfit, quadfit(xfit), color=colors[i % len(colors)], ls='dashed') 
    ix = quadratic_intersections(poly, [-0.5])
    #plt.scatter(*ix, marker='x', color=colors[i % len(colors)], s=40, linewidth=2)
    
    lfit = savitzky_golay(nl, 7, 3)
    plt.plot(e, lfit, color='k', ls='dashed')
    
    xfit = np.linspace(0.01, 0.39, 1000)
    f = interpolate.interp1d(e,lfit, kind='quadratic')
    fgrid = f(xfit)
    maxi = np.argmax(fgrid)
    maxl = fgrid[maxi]
    maxe = xfit[maxi]
    fgrid1s = fgrid-maxl+0.5
    dotp = fgrid1s[0:-1]*fgrid1s[1:]
    #print maxi, maxl, maxe
    #print dotp
    errsi = np.where(dotp < 0)
    #print errsi
    evalerrs = []
    for erri in errsi[0]:
      evalerr = (xfit[erri]+xfit[erri+1])/2.
      evalerrs.append(evalerr)
    if len(evalerrs) == 1:
      if evalerrs[0] > maxe:
        evalerrs = [0.] + evalerrs
      else:
        evalerrs = evalerrs + [0.]
    if len(evalerrs) == 0:
      evalerrs = [0., 0.4]
    if len(evalerrs) > 2:
      #raise
      pass

    print "this cut"
    print cutlist, i
    print evalerrs[0], maxe, evalerrs[1]

    plt.scatter([evalerrs[0], maxe, evalerrs[1]], [maxl-0.5,maxl-0,maxl-0.5], marker='x', color=colors[i % len(colors)], s=40, linewidth=2)
    plt.plot(xfit, f(xfit), color='gray', ls='dotted')
    
    maxes.append(maxe)
    maxeerrs1.append(evalerrs[0])
    maxeerrs2.append(evalerrs[1])

  plt.figure(2)
  lc = len(cutlist)
  if lc > 1:
    xbinvals = np.array(np.array(cutlist[0:-1])+np.array(cutlist[1:]))/2.
    print xbinvals
    xbinvals = [xbinvals[0]+(cutlist[0]-xbinvals[0])*2] + [val for val in xbinvals]
    xbinvals = [val for val in xbinvals] + [xbinvals[-1]-(-cutlist[-1]+xbinvals[-1])*2]
    print xbinvals
  else:
    xbinvals = [-1,1]
  xbinvals = np.array(xbinvals)
  plt.errorbar(xbinvals, maxes, yerr=[np.array(maxes)-np.array(maxeerrs1), np.array(maxeerrs2)-np.array(maxes)], fmt='o') 
  outputname1 = propertystring+'_trend_R<15.png'
  plt.xlabel(propertystring, fontsize=12)
  plt.ylabel('$\sigma_e$', fontsize=12)
  plt.ylim((0.,0.4))
  plt.savefig('figs4/'+outputname1, bbox_inches='tight')
  plt.close()
 


  plt.figure(1)
  plt.plot([min(e), max(e)], [-2, -2], color='gray', ls='dotted')
  plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand")
  plt.xlabel('$\sigma_e$', fontsize=12)
  plt.ylabel('Natural Log Likelihood Difference', fontsize=12)
  plt.grid()

  plt.xlim((0.,0.4))
  plt.ylim((-20, 5))
  
  outputname1 = propertystring+'_R<15.png'
  plt.savefig('figs4/'+outputname1, bbox_inches='tight')
  plt.close()

  print "plotted "+propertystring

  return 





################
# Make plot for all planets, with uncertainties

plt.figure(1)

nlall = allllike
allfit = savitzky_golay(nlall, 7, 3)
offset = np.max(allfit)
allfit -= offset


plt.plot(e, allllike - offset, 'o', color='k', label='All') 

plt.plot(e, allfit, color='k')
plt.plot([min(e), max(e)], [-0.5, -0.5], color='gray', ls='dotted')
plt.plot([min(e), max(e)], [-2, -2], color='gray', ls='-.')
plt.plot([min(e), max(e)], [-4.5, -4.5], color='gray', ls='dashed')

x2 = np.linspace(0., max(e), num=1000, dtype=np.float)
y2 = np.interp(x2, e, allfit) 
plt.fill_between(x2, y2*0-40, y2, where=y2 >= -4.5, facecolor='green', interpolate=True)
plt.fill_between(x2, y2*0-40, y2, where=y2 >= -0.5, facecolor='blue', interpolate=True)

print "bestfit = ", x2[np.argmax(y2)]
print "1 sigma = ", x2[np.where(y2 >= -0.5)] 
print "2 sigma = ", x2[np.where(y2 >= -2)] 
print "3 sigma = ", x2[np.where(y2 >= -4.5)] 


#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand")
plt.xlabel('$\sigma_e$', fontsize=12)
plt.ylabel('Natural Log Likelihood Difference', fontsize=12)
plt.grid()

plt.xlim((0.,0.4))
plt.xlim((0.1,0.3))
plt.ylim((-40, 8))

outputname1 = 'All_R<15.png'
plt.savefig('figs4/'+outputname1, bbox_inches='tight')
plt.close()


# Make plot for all planets, with uncertainties
# Also add multiplanet system line

plt.figure(1)

nlall = allllike
allfit = savitzky_golay(nlall, 7, 3)
offset = np.max(allfit)
allfit -= offset


plt.plot(e, allllike - offset, 'o', color='k', label='All') 

plt.plot(e, allfit, color='k')
plt.plot([min(e), max(e)], [-0.5, -0.5], color='gray', ls='dotted')
plt.plot([min(e), max(e)], [-2, -2], color='gray', ls='-.')
plt.plot([min(e), max(e)], [-4.5, -4.5], color='gray', ls='dashed')

x2 = np.linspace(0., max(e), num=1000, dtype=np.float)
y2 = np.interp(x2, e, allfit) 
plt.fill_between(x2, y2*0-40, y2, where=y2 >= -4.5, facecolor='green', interpolate=True)
plt.fill_between(x2, y2*0-40, y2, where=y2 >= -0.5, facecolor='blue', interpolate=True)

print "bestfit = ", x2[np.argmax(y2)]
print "1 sigma = ", x2[np.where(y2 >= -0.5)] 
print "2 sigma = ", x2[np.where(y2 >= -2)] 
print "3 sigma = ", x2[np.where(y2 >= -4.5)] 


#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand")
plt.xlabel('$\sigma_e$', fontsize=12)
plt.ylabel('Natural Log Likelihood Difference', fontsize=12)
plt.grid()

gridm = np.loadtxt('./fgrid.txt')
gride = np.loadtxt('./xgrid.txt')
gridi = np.loadtxt('./ygrid.txt')

imaxi, emaxi = np.unravel_index(np.argmax(gridm, axis=None), gridm.shape)

evals = gride
llvalsm = gridm[imaxi]
llvalsm -= np.max(llvalsm)

plt.plot(evals, llvalsm, color='c', alpha=0.5)
#x2 = np.linspace(0., max(e), num=1000, dtype=np.float)
#y2 = np.interp(x2, e, allfit) 
plt.fill_between(evals, llvalsm*0-35, llvalsm, where=llvalsm >= -4.5, facecolor='green', interpolate=True, alpha=0.5)
plt.fill_between(evals, llvalsm*0-35, llvalsm, where=llvalsm >= -0.5, facecolor='blue', interpolate=True, alpha=0.5)


plt.xlim((0.,0.4))
plt.xlim((0.0,0.3))
plt.ylim((-34.5, 8))


outputname1 = 'All_R<15_plusmultis.png'
plt.savefig('figs4/'+outputname1, bbox_inches='tight')
plt.close()


