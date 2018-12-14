## Plot the e-i contour plot for all multiplanet systems which pass our quality cuts
##

import numpy as np
import matplotlib.pyplot as plt
import random
import scipy
from scipy import stats
import os
from scipy import interpolate
import scipy.signal

##########
from s_g_funs import * 
from cksecc_helper import *


cks_data, koi_data, koi_errs = load_data()
koilist, propertymatrix = divide_into_singles_and_multis(cks_data, koi_data, koi_errs, singles=False) 

e, inc, like = read_in_multis(group1=True)

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
like = propertymatrix[-1]


print len(like)
print like
point5=[2,3]
point5 = [2,3,4,5]
#point5=[]
ncuts=7
ncuts2=1
nincut = 10
e = e.reshape(ncuts+len(point5), nincut*ncuts2) 
inc = inc.reshape(ncuts+len(point5), nincut*ncuts2) 
like = like.reshape(len(like), ncuts+len(point5), nincut*ncuts2) 


alltotal = np.sum(np.log(like), axis=0)
allmaxllike = np.max(alltotal)
allllike = alltotal-allmaxllike

print "e inc:"
print e
print inc

print "llike"
print allllike

print "alltotal"
print alltotal
print "max log like"
print allmaxllike


##########
# all multis fit
nl = allllike 

print nl.shape

lfit = sgolay2d(nl, 5, 3)

print "lfit"
print lfit


imax = 0.061
emax = 0.10

xfit = np.linspace(0.0, emax-0.0001, 1000)
yfit = np.linspace(0.0, imax-0.0001, 1000)
f = interpolate.interp2d(e[0,:], inc[:,0], lfit, kind='cubic')
f1 = interpolate.interp2d(e[0,:], inc[:,0], allllike, kind='cubic')
plt.clf()
#cplotf = plt.contourf(e[0,:], inc[:,0], allllike, levels=np.linspace(-500,1,num=500))
#cplotf = plt.contourf(xfit, yfit, f, levels=np.linspace(-500,1,num=500))
#cplotf = plt.contourf(xfit, yfit, f1(xfit, yfit), levels=np.linspace(-500,1,num=500))
#plt.show()
#exit()
fgrid = f1(xfit, yfit)
maxi = np.argmax(fgrid)
maxie = maxi % len(xfit) 
maxii = maxi / len(yfit)
maxl = fgrid.flatten()[maxi]
maxe = xfit[maxie]
maxi = yfit[maxii]

print "bestfit e = ", maxe
print "bestfit i = ", maxi


fgrid -= maxl

plt.clf()

print fgrid
print fgrid.shape

cplot = plt.contour(xfit, yfit, fgrid, levels=[-4.5, -2., -.5], colors='k', linestyles=['dotted','dashdot','dashed'])
#cplot = plt.contour(xfit, yfit, fgrid, levels=[-25, -10., -2], colors='k', ls=['-','.-','.'])
#cplotf = plt.contourf(xfit, yfit, fgrid, levels=np.flip(np.arange(100)*-1,0))
cplotf = plt.contourf(xfit, yfit, fgrid, levels=np.linspace(-200,1,num=500))
bf = [maxe, maxi]
plt.plot([bf[0]], [bf[1]], marker='X')
colorbar = plt.colorbar(cplotf, shrink=0.8, extend='both')
colorbar.set_label('Log(Likelihood)')
plt.xlabel('$\sigma_e$', fontsize=12)
plt.ylabel('$\sigma_i$', fontsize=12)
    
fgrid1s = fgrid+0.5
fgrid2s = fgrid+2.0

np.savetxt("ygrid.txt", yfit)
np.savetxt("xgrid.txt", xfit)
np.savetxt("fgrid.txt", fgrid)

evalerrsmaster=[1000, -1000]
for i in range(len(yfit)):
  dotp = fgrid1s[i][0:-1]*fgrid1s[i][1:]
  errsi = np.where(dotp < 0)
  evalerrs = []
  for erri in errsi[0]:
    evalerr = (xfit[erri]+xfit[erri+1])/2.
    evalerrs.append(evalerr)
  if len(evalerrs) == 1:
    if evalerrs[0] > maxe:
      evalerrs = [0.] + evalerrs
    else:
      evalerrs = evalerrs + [0.]
  if len(evalerrs) > 2:
    raise
  if len(evalerrs) == 2:
    if evalerrs[0] < evalerrsmaster[0]:
      evalerrsmaster[0] = evalerrs[0]
    if evalerrs[1] > evalerrsmaster[1]:
      evalerrsmaster[1] = evalerrs[1]

print evalerrsmaster
plt.plot([evalerrsmaster[0], evalerrsmaster[0]], [0,imax], color='r')
plt.plot([evalerrsmaster[1], evalerrsmaster[1]], [0,imax], color='r')
evalerrsmaster=[1000, -1000]
for i in range(len(yfit)):
  dotp = fgrid2s[i][0:-1]*fgrid2s[i][1:]
  errsi = np.where(dotp < 0)
  evalerrs = []
  for erri in errsi[0]:
    evalerr = (xfit[erri]+xfit[erri+1])/2.
    evalerrs.append(evalerr)
  if len(evalerrs) == 1:
    if evalerrs[0] > maxe:
      evalerrs = [0.] + evalerrs
    else:
      evalerrs = evalerrs + [0.]
  if len(evalerrs) > 2:
    raise
  if len(evalerrs) == 2:
    if evalerrs[0] < evalerrsmaster[0]:
      evalerrsmaster[0] = evalerrs[0]
    if evalerrs[1] > evalerrsmaster[1]:
      evalerrsmaster[1] = evalerrs[1]
print "2-sigma uncertainties:"
print "e ", [evalerrsmaster[0], evalerrsmaster[1]]


ivalerrsmaster=[1000, -1000]
for i in range(len(xfit)):
  dotp = fgrid1s[:,i][0:-1]*fgrid1s[:,i][1:]
  errsi = np.where(dotp < 0)
  ivalerrs = []
  for erri in errsi[0]:
    ivalerr = (yfit[erri]+yfit[erri+1])/2.
    ivalerrs.append(ivalerr)
  if len(ivalerrs) == 1:
    if ivalerrs[0] > maxe:
      ivalerrs = [0.] + ivalerrs
    else:
      ivalerrs = ivalerrs + [0.]
  if len(ivalerrs) > 2:
    raise
  if len(ivalerrs) == 2:
    if ivalerrs[0] < ivalerrsmaster[0]:
      ivalerrsmaster[0] = ivalerrs[0]
    if ivalerrs[1] > ivalerrsmaster[1]:
      ivalerrsmaster[1] = ivalerrs[1]
print "1-sigma uncertainties:"
print "i ", [ivalerrsmaster[0], ivalerrsmaster[1]]

print ivalerrsmaster
plt.plot([0,emax], [ivalerrsmaster[0], ivalerrsmaster[0]], color='b')
plt.plot([0,emax], [ivalerrsmaster[1], ivalerrsmaster[1]], color='b')

ivalerrsmaster=[1000, -1000]
for i in range(len(xfit)):
  dotp = fgrid2s[:,i][0:-1]*fgrid2s[:,i][1:]
  errsi = np.where(dotp < 0)
  ivalerrs = []
  for erri in errsi[0]:
    ivalerr = (yfit[erri]+yfit[erri+1])/2.
    ivalerrs.append(ivalerr)
  if len(ivalerrs) == 1:
    if ivalerrs[0] > maxe:
      ivalerrs = [0.] + ivalerrs
    else:
      ivalerrs = ivalerrs + [0.]
  if len(ivalerrs) > 2:
    raise
  if len(ivalerrs) == 2:
    if ivalerrs[0] < ivalerrsmaster[0]:
      ivalerrsmaster[0] = ivalerrs[0]
    if ivalerrs[1] > ivalerrsmaster[1]:
      ivalerrsmaster[1] = ivalerrs[1]
print "2-sigma uncertainties:"
print "i ", [ivalerrsmaster[0], ivalerrsmaster[1]]
print ivalerrsmaster

print "ALL PLANETS PLOT"
plt.ylim((0.,imax))
plt.xlim((0.,emax))
plt.savefig('figs4/multis_2d.png')




