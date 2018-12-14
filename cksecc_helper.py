import numpy as np
from scipy import interpolate
from collections import Counter
from s_g_funs import *
import matplotlib.pyplot as plt
import scipy
from scipy import stats
import os
from scipy import interpolate
import radvel.utils

# read in planet data from files
def load_data(nolow=False):
  cks_fname = 'resources/CKS_GAIA_Stellar_Parameters.txt'
  cks_data = np.loadtxt(cks_fname, usecols=[0,1,2,3,4,5,6,7,8,9,10])#, usecols=[0,11,12,13,14,15,16])
  cks_data = np.transpose(cks_data)
  koi_fname = "resources/KOI_list_SNR_CONFIRMED.tsv" 
  if nolow:
    koi_fname = "resources/KOI_NOLOWSNR.tsv"
  koi_data = np.loadtxt(koi_fname, usecols=[2,3,12,18,24])
  koi_errs = np.loadtxt(koi_fname, usecols=[2,4,5,13,14,19,20])
  koi_data = np.transpose(koi_data)
  koi_errs = np.transpose(koi_errs)
  
  # Sort KOI data in order of descending period by KOI number
  koistar = np.floor(koi_data[0])
  periods = koi_data[1]
  koipers = [(koistar[i], 1./periods[i]) for i in range(len(koistar))]
  koipers2 = np.array(koipers, dtype=[('koinum', 'f'), ('invper', 'f')])
  koipersi = np.argsort(koipers2, order=('koinum', 'invper'))
  
  koi_data = np.array([ki[koipersi] for ki in koi_data])
  koi_errs = np.array([ki[koipersi] for ki in koi_errs])
  
  # process data so that only cks data with valid candidates is used
  bothkois_i1 = np.nonzero(np.in1d(cks_data[0], np.floor(koi_data[0])))
  cks_data = np.array([cks_data_column[bothkois_i1] for cks_data_column in cks_data])

  return cks_data, koi_data, koi_errs



# Returns relevant planet properties for all stars in Kepler KOI catalog
def get_kepler_transit_durations(data):

  kois = data[0]
  periods = data[1]
  durations = data[2]/24.
  rpors = data[3]
  snr = data[4]
  return kois, periods, durations, rpors, snr



# Returns planet property uncertainties for all planets in the KOI catalog
def get_kepler_transit_errors(errors):
  kois = errors[0]
  period_e_hi = errors[1]
  period_e_lo = errors[2]
  durations_e_hi = errors[3]/24.
  durations_e_lo = errors[4]/24.
  rpors_e_hi = errors[5]
  rpors_e_lo = errors[6]
  return kois, period_e_hi, period_e_lo, durations_e_hi, durations_e_lo, rpors_e_hi, rpors_e_lo


def divide_into_singles_and_multis(cks_data, koi_data, koi_errs, singles=True, more_properties=False): 
  kep_kois, kep_periods, kep_durations, kep_rpors, kep_snrs = get_kepler_transit_durations(koi_data)
  kois_int = np.floor(kep_kois)
  uniquelist, indices, inverse, counts = np.unique(kois_int, return_index=True, return_inverse=True, return_counts=True)
  kep_multisi = counts > 1
  kep_singlesi = counts == 1
  kep_singles = uniquelist[kep_singlesi]
  kep_multis = uniquelist[kep_multisi]
  
  allkois = cks_data[0]
  cks_singlesi = np.in1d(allkois, kep_singles)
  cks_singles = allkois[cks_singlesi]
  cks_multisi = np.in1d(allkois, kep_multis)
  cks_multis = allkois[cks_multisi]
 
  if singles==True:
    sset = cks_singles
    sdex1 = cks_singlesi
  else:
    sset = cks_multis
    sdex1 = cks_multisi
  
  dexi = np.in1d(np.floor(kois_int), sset).nonzero()

  # Get koi information for this group of planets
  skois = kep_kois[dexi]
  speriods = kep_periods[dexi]
  srpors = kep_rpors[dexi]
  if more_properties:
    sdurations = kep_durations[dexi]
    srpors = kep_rpors[dexi]
    ssnrs = kep_snrs[dexi]

  # Sort by KOI
  sorti = np.argsort(skois)
  
  # Sort all the koi information
  sskois = skois[sorti]
  ssperiods = speriods[sorti]
  ssrpors = srpors[sorti]
  if more_properties:
    ssdurations = sdurations[sorti]
    ssrpors = srpors[sorti]
    sssnrs = ssnrs[sorti]

  # Get the stellar information for this group of planets
  # Then sort that data by KOI as well
  ckssorti = np.argsort(sset)
  rstar = cks_data[4]
  srstar = rstar[sdex1]
  ssrstar = srstar[ckssorti]
  mstar = cks_data[1]
  smstar = mstar[sdex1]
  ssmstar = smstar[ckssorti]
  teff = cks_data[7]
  steff = teff[sdex1]
  ssteff = steff[ckssorti]
  feh = cks_data[8]
  sfeh = feh[sdex1]
  ssfeh = sfeh[ckssorti]
  logage = cks_data[9]
  slogage = logage[sdex1]
  sslogage = slogage[ckssorti]
  dilute = cks_data[10]
  sdilute = dilute[sdex1]
  ssdilute = sdilute[ckssorti]

  # If multis, must copy stellar properties for each planet
  if singles==False:
    longintkoilist_i = np.in1d(np.floor(sskois), kep_multis)
    longintkoilist = np.floor(sskois[longintkoilist_i])
    uniquelistmultis, indicesmultis, inversemultis, countsmultis = np.unique(np.floor(sskois), return_index=True, return_inverse=True, return_counts=True)
    i=0
    longintkoilistmap=[]
    for j in countsmultis:
      for k in range(j):
        longintkoilistmap.append(i)
      i+=1
    longintkoilistmap = np.array(longintkoilistmap)

    ssrstar = ssrstar[longintkoilistmap]
    ssmstar = ssmstar[longintkoilistmap]
    sslogage = sslogage[longintkoilistmap]
    ssteff = ssteff[longintkoilistmap]
    ssfeh = ssfeh[longintkoilistmap]
    ssdilute = ssdilute[longintkoilistmap]

  # Rename conveniently
  kois = sskois
  msun = ssmstar
  rsun = ssrstar
  dilution = ssdilute
  age = sslogage
  teff = ssteff
  metal = ssfeh
  period = ssperiods
  if more_properties:
    durations = ssdurations
    rpors = ssrpors
    snr = sssnrs

  # Derive Physical Radii
  REARTHRSUN = 0.009168
  radius = ssrpors*ssrstar/REARTHRSUN

  # some more properties only applicable to multis
  if singles==False:
    # Multiplicity of the system
    multiplicity = np.zeros(len(radius))
    i=0
    for j in countsmultis:
      multiplicity[i:i+j] = j
      i+=j
    # Radius of largest planet in system  
    radiusbig = np.zeros(len(radius))
    i=0
    for j in countsmultis:
      rbi = max(radius[i:i+j])
      radiusbig[i:i+j] = rbi
      i+=j
    # Period ratio of nearest interior planet to a given planet
    # and period ratio f nearest outer planet to a given planet
    periodratio_in = np.zeros(len(period))
    periodratio_out = np.zeros(len(period))
    i=0
    for j in countsmultis:
      periodshere = np.sort(period[i:i+j])
      for k in range(j-1):
        periodratio_in[i+k+1] = periodshere[k]/periodshere[k+1]
        periodratio_out[i+k] = periodshere[k+1]/periodshere[k] 
      periodratio_in[i] = np.nan #innermost planet in system
      periodratio_out[i+j-1] = np.nan #outermost planet in system
      i+=j


  propertymat = [msun, rsun, period, radius, teff, metal, age, dilution]
  if more_properties:
    propertymat.append(durations)
    propertymat.append(rpors)
    propertymat.append(snr)
  if singles==False:
    propertymat.append(multiplicity)
    propertymat.append(radiusbig)
    propertymat.append(periodratio_in)
    propertymat.append(periodratio_out)

  return kois, propertymat 
  


# Remove Bryson FPs
def remove_fps(koilist, propertymatrix, singles=True):

  len1 = len(koilist)
  
  brysondat = np.genfromtxt('resources/Morton_2016_Table_6_FPPprob.txt', dtype=None, names=('koi', 'flag', 'fpp', 'fppe'))
  b_kois = brysondat['koi']
  b_flag = brysondat['flag']
  b_fpp = brysondat['fpp']
  
  fpflagi = np.where(b_flag != 'FP')
  fpflagj = np.where(b_flag == 'FP')
  
  b_koisflagCP = b_kois[fpflagi]
  b_koisflagFP = b_kois[fpflagj]
  b_koisflagi = np.in1d(koilist, b_koisflagCP)
  b_koisflagj = np.in1d(koilist, b_koisflagFP)
  
  potentialFPcuts = np.where( np.logical_and( b_flag != 'PL', b_fpp > 0.5 ) )
  b_kois_potentialFP = b_kois[potentialFPcuts]
  
  #b_koisflag_potFP = np.in1d(koilist, b_kois_potentialFP)
  b_koisflag_SAFE = np.in1d(koilist, b_kois_potentialFP, invert=True)

  len2 = len(koilist[b_koisflag_SAFE])

  print str(len1-len2)+" False Positives removed"
  print len(koilist[b_koisflag_SAFE])
  #print [len(propertylist) for propertylist in propertymatrix]
  #print [propertylist[b_koisflag_SAFE] for propertylist in propertymatrix]

  if singles==False:
    # Make sure no multis became singles due to FPs 
    newkoilist = koilist[b_koisflag_SAFE]
    singlekois = np.array([item for item, count in Counter(np.floor(newkoilist)).iteritems() if count == 1])
    singlemask = np.in1d(np.floor(koilist), singlekois, invert=True)
    b_koisflag_SAFE = np.logical_and(b_koisflag_SAFE, singlemask)
    len3 = len(koilist[b_koisflag_SAFE])
    print "here"
    print np.shape(b_koisflag_SAFE)
    print np.shape(singlemask)
  
    print str(len1-len3)+" False Positives + neighbors removed"
    print len(koilist[b_koisflag_SAFE])

  return koilist[b_koisflag_SAFE], [propertylist[b_koisflag_SAFE] for propertylist in propertymatrix]




def make_parameter_cuts2(koilist, propertymatrix, index, lessthan=None, greaterthan=None, singles=True):

  print len(koilist)
  len1 = len(koilist)
  if lessthan==None and isinstance(greaterthan, (int, long, float, complex)):
    goodmask = np.where(propertymatrix[index] > greaterthan)
  elif isinstance(lessthan, (int, long, float, complex)) and greaterthan==None:
    goodmask = np.where(propertymatrix[index] < lessthan)
  else:
    raise

  goodkoilist = koilist[goodmask]
  goodsubset = np.in1d(koilist, goodkoilist)
  badsubset = np.in1d(koilist, goodkoilist, invert=True)
  print "removed:"
  print koilist[badsubset]

  len2 = len(koilist[goodsubset])
  print str(len1-len2)+" kois removed"
  if singles==False:
    # Make sure no multis became singles 
    newkoilist = koilist[goodsubset]
    singlekois = np.array([item for item, count in Counter(np.floor(newkoilist)).iteritems() if count == 1])
    singlemask = np.in1d(np.floor(koilist), singlekois, invert=True)
    goodsubset = np.logical_and(goodsubset, singlemask)
    len3 = len(koilist[goodsubset])
    print str(len1-len3)+" kois + neighbors removed"

  return koilist[goodsubset], [propertylist[goodsubset] for propertylist in propertymatrix]






# Clean the sample
def clean_sample(koilist, propertymatrix, customremove=True, singles=True):
  
  # Remove FPs
  new_koilist, new_propertymatrix = remove_fps(koilist, propertymatrix, singles=singles)

  # Remove some KOIs by number
  if customremove:
    if singles==True:
      notmulti = np.where(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.floor(new_koilist) != 4368, np.floor(new_koilist) != 3283), np.floor(new_koilist) != 3184), np.floor(new_koilist) != 2590), np.floor(new_koilist) != 2768), np.floor(new_koilist) != 379)) 
      new_koilist = new_koilist[notmulti]
      new_propertymatrix = [propertylist[notmulti] for propertylist in new_propertymatrix]
    if True:
      # remove KOIs where the uncertainty in rpors > 2% of the stellar radius
      # Singles:  3913.01 1800.01 4368.02 3811.01 2398.01 4546.01 5284.01 3891.01
      # Multis: 1707.02  620.03 1102.01 1203.03
      notuncert = np.where(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.floor(new_koilist) != 3913, np.floor(new_koilist) != 1800), np.floor(new_koilist) != 4368), np.floor(new_koilist) != 3811), np.floor(new_koilist) != 2398), np.floor(new_koilist) != 4546)) 
      new_koilist = new_koilist[notuncert]
      new_propertymatrix = [propertylist[notuncert] for propertylist in new_propertymatrix]
      # cont
      notuncert = np.where(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.floor(new_koilist) != 5284, np.floor(new_koilist) != 3891), new_koilist != 1707.02), new_koilist != 620.03), new_koilist != 1102.01), new_koilist != 1203.03)) 
      new_koilist = new_koilist[notuncert]
      new_propertymatrix = [propertylist[notuncert] for propertylist in new_propertymatrix]


  # remove Rp > 15
  new_koilist, new_propertymatrix = make_parameter_cuts2(new_koilist, new_propertymatrix, 3, lessthan=15., singles=singles)  
  # remove Rsun > 2
  new_koilist, new_propertymatrix = make_parameter_cuts2(new_koilist, new_propertymatrix, 1, lessthan=2., singles=singles)  
  # remove Rsun < 0.5 
  new_koilist, new_propertymatrix = make_parameter_cuts2(new_koilist, new_propertymatrix, 1, greaterthan=0.5, singles=singles)  
  # remove dilution > 5% 
  new_koilist, new_propertymatrix = make_parameter_cuts2(new_koilist, new_propertymatrix, 7, lessthan=1.05, singles=singles)  


  #return kois, [msun, rsun, period, radius, teff, metal, age, dilution]
  return new_koilist, new_propertymatrix



# Read in single-planet system data
def read_in_singles(): 
  emaxval=0.69
  nper = 100 # ndata sets per file name
  nincut = 10
  ncuts = 7
  point5 = [0,1,2]
  point25 = [1]
  point75 = [1]
  filenamepart = ['data3/singles_m12_10000_'+str(i+1)+'' for i in range(ncuts)]
  [filenamepart.append('data3/singles_m12_10000_'+str(i+1)+'.5') for i in range(ncuts)]
  [filenamepart.append('data3/singles_m12_10000_'+str(i+1)+'.25') for i in range(ncuts)]
  [filenamepart.append('data3/singles_m12_10000_'+str(i+1)+'.75') for i in range(ncuts)]
  like, koi, e, inc = [], [], [], []
  for i in range(ncuts):
    filenameparti = filenamepart[i]
    for j in range(nincut):
      likeij=[]
      incij=0.
      eij=j/100. + i/10.
      suffixj = '_%f_%f_10000.txt' % (eij, incij)
      for k in range(nper):
        kpart = '_%s.txt' % str(k)
        fname = filenameparti + kpart + suffixj
        dataijk = np.loadtxt(fname)
        likeij.append(dataijk[:,1])
      like.append(np.median(likeij, axis=0))
      e.append(eij)
      inc.append(incij)
    if i in point25:
      filenameparti = filenamepart[i+ncuts*2]
      for j in range(nincut):
        likeij=[]
        incij=0.
        eij=j/100. + i/10. + 0.0025
        suffixj = '_%f_%f_10000.txt' % (eij, incij)
        for k in range(nper):
          kpart = '_%s.txt' % str(k)
          fname = filenameparti + kpart + suffixj
          dataijk = np.loadtxt(fname)
          likeij.append(dataijk[:,1])
        like.append(np.median(likeij, axis=0))
        e.append(eij)
        inc.append(incij)
    if i in point5:
      filenameparti = filenamepart[i+ncuts*1]
      for j in range(nincut):
        likeij=[]
        incij=0.
        eij=j/100. + i/10. + 0.005
        suffixj = '_%f_%f_10000.txt' % (eij, incij)
        for k in range(nper):
          kpart = '_%s.txt' % str(k)
          fname = filenameparti + kpart + suffixj
          dataijk = np.loadtxt(fname)
          likeij.append(dataijk[:,1])
        like.append(np.median(likeij, axis=0))
        e.append(eij)
        inc.append(incij)
    if i in point75:
      filenameparti = filenamepart[i+ncuts*3]
      for j in range(nincut):
        likeij=[]
        incij=0.
        eij=j/100. + i/10. + 0.0075
        suffixj = '_%f_%f_10000.txt' % (eij, incij)
        for k in range(nper):
          kpart = '_%s.txt' % str(k)
          fname = filenameparti + kpart + suffixj
          dataijk = np.loadtxt(fname)
          likeij.append(dataijk[:,1])
        like.append(np.median(likeij, axis=0))
        e.append(eij)
        inc.append(incij)
  like = np.transpose(like)
  koi = dataijk[:,0]
  
  e = np.array(e)
  inc = np.array(inc)
  
  esorti = np.argsort(e)
  e = e[esorti]
  likesort = np.array([likei[esorti] for likei in like])
  like = likesort
 
  return e, inc, like, emaxval
 



def read_in_multis(zoom=False, alle=False, group1=False):
  #### Combining data
  nper = 19  # ndata sets per file name
  nincut = 10
  #ncuts = 7
  ncuts = 5
  if alle:
    ncuts = 5
  #point5i = [2,3,4,5]
  point5i = [2,3]
  #ncuts2 = 1
  ncuts2 = 2
  if group1:
    point5i = [2,3,4,5]
    ncuts2 = 1
    ncuts = 7
  if alle:
    ncuts2 = 2
  filenamepart = ['data3_BAD/multis_m12_10000_'+'1_'+str(i+1)+'' for i in range(ncuts)]
  filenamepart2 = ['data3_BAD/multis_m12_10000_'+'2_'+str(i+1)+'' for i in range(ncuts)]
  [filenamepart.append('data3_BAD/multis_m12_10000_'+'1_'+str(i+1)+'.5') for i in range(ncuts)]
  [filenamepart2.append('data3_BAD/multis_m12_10000_'+'2_'+str(i+1)+'.5') for i in range(ncuts)]
  filenamepart = ['data3/multis_m12_10000_'+'1_'+str(i+1)+'' for i in range(ncuts)]
  filenamepart2 = ['data3/multis_m12_10000_'+'2_'+str(i+1)+'' for i in range(ncuts)]
  [filenamepart.append('data3/multis_m12_10000_'+'1_'+str(i+1)+'.5') for i in range(ncuts)]
  [filenamepart2.append('data3/multis_m12_10000_'+'2_'+str(i+1)+'.5') for i in range(ncuts)]
  #for i in range(ncuts):
  #  filenamepart.append('data3/multis_m12_10000_'+'1_'+str(i+1)+'')
  #  filenamepart.append('data3/multis_m12_10000_'+'2_'+str(i+1)+'')
  filenameparts = [filenamepart, filenamepart2]

  if zoom:
    nper = 11
    nincut = 10
    ncuts = 6
    point5i = []
    ncuts2 = 1
    filenamepart = ['data3/multis_m12_10000_'+'1_3.'+str(i+4)+'' for i in range(ncuts-2)]
    filenamepart.append('data3/multis_m12_10000_'+'1_4')
    filenamepart.append('data3/multis_m12_10000_'+'1_4.5')
    filenameparts = [filenamepart]

  like, koi, e, inc = [], [], [], []
  for i in range(ncuts):
    for qq in range(ncuts2): 
      filenameparti = filenameparts[qq][i]
      for j in range(nincut):
        likeij=[]
        incij= i/100.
        if zoom:
          if i < 4:
            incij = 0.02+(4+i)/1000.
          elif i==4:
            incij = 0.03
          elif i==5:
            incij = 0.035
        eij= j/100. + 0.1*qq
        suffixj = '_%f_%f_10000.txt' % (eij, incij)
        for k in range(nper):
          kpart = '_%s.txt' % str(k)
          fname = filenameparti + kpart + suffixj
          dataijk = np.loadtxt(fname)
          likeij.append(dataijk[:,1])
        like.append(np.median(likeij, axis=0))
        e.append(eij)
        inc.append(incij)
    if i in point5i:
      for qq in range(ncuts2):
        filenameparti = filenameparts[qq][i+ncuts]
        for j in range(nincut):
          likeij=[]
          incij= i/100. + 0.005
          eij= j/100. + 0.1*qq
          suffixj = '_%f_%f_10000.txt' % (eij, incij)
          for k in range(nper):
            kpart = '_%s.txt' % str(k)
            fname = filenameparti + kpart + suffixj
            dataijk = np.loadtxt(fname)
            likeij.append(dataijk[:,1])
          like.append(np.median(likeij, axis=0))
          e.append(eij)
          inc.append(incij)
#      if i in point5i:
#        filenameparti = filenameparts[qq][i+ncuts]
#        for j in range(nincut):
#          likeij=[]
#          incij= i/100. + 0.005
#          eij= j/100. + 0.1*qq
#          suffixj = '_%f_%f_10000.txt' % (eij, incij)
#          for k in range(nper):
#            kpart = '_%s.txt' % str(k)
#            fname = filenameparti + kpart + suffixj
#            dataijk = np.loadtxt(fname)
#            likeij.append(dataijk[:,1])
#          like.append(np.median(likeij, axis=0))
#          e.append(eij)
#          inc.append(incij)

  like = np.transpose(like)
  koi = dataijk[:,0]
  
  e = np.array(e)
  inc = np.array(inc)

  print e
  print inc
  #exit()

  return e, inc, like





### Analysis functions:

# Find best fit e for each koi
def sg_fit(e, normlikes, emaxval, eminval=0.001):
  xfit = np.linspace(eminval, emaxval, 10000)
  maxes = []
  for i in range(len(normlikes)):
    nl = normlikes[i]
    nl10i = np.where(nl>-10.)
    e10 = e[nl10i]
    nl10 = nl[nl10i]
  
    lfit = savitzky_golay(nl, 7, 3)
    f = interpolate.interp1d(e, lfit, kind='quadratic')
    fgrid = f(xfit)
    maxi = np.argmax(fgrid)
    maxl = fgrid[maxi]
    maxe = xfit[maxi]
    maxes.append(maxe)

  return maxes


## Make property distribution plots for high and low e samples
def twosampleplot(plist, propertystring, eorderi, splitoptionsindex, isplit, lowstrengthsort, highstrengthsort, log=False, nbin=20, pctile=1, mask=False, highonly=False, singles=True, xlabel=None, legendon=False):
 
  if singles==True:
    strstart = 'singles'
  else:
    strstart = 'multis'

  if log:
    plist = np.log10(plist)
  lowedat = plist[eorderi][:splitoptionsindex[isplit]][lowstrengthsort]
  if not highonly:
    lowedat = lowedat[:int(len(lowedat)*pctile)]
  highedat = plist[eorderi][splitoptionsindex[isplit]:][highstrengthsort]
  highedat = highedat[int(len(highedat)-len(highedat)*pctile):]
  bins=np.linspace(np.nanmin(plist), np.nanmax(plist), nbin)

  xlims = [0,0]
  if np.nanmin(plist) < 0:
    xlims[0] = np.nanmin(plist)*1.02
  else:
    xlims[0] = np.nanmin(plist)*0.98
  if np.nanmax(plist) < 0:
    xlims[1] = np.nanmax(plist)*0.98
  else:
    xlims[1] = np.nanmax(plist)*1.02

  plt.figure(2)
  lowedatmask=np.isfinite(lowedat)
  lowedat=lowedat[lowedatmask]
  highedatmask=np.isfinite(highedat)
  highedat=highedat[highedatmask]

  plt.hist(lowedat, color='b', alpha=0.5, edgecolor='b', normed=True, bins=bins, label='Low e') 
  plt.hist(highedat, color='r', alpha=0.5, edgecolor='r', normed=True, bins=bins, label='High e') 
  propertystringf = propertystring
  if log:
    propertystringf = 'Log10['+propertystring+']'
    locs, labels = plt.xticks()
    plt.xticks(locs, [radvel.utils.round_sig(v) for v in [10.**loc for loc in locs]])
  if pctile < 1:
    propertystringf += '_pct_'+str(pctile)
  if not highonly:
    propertystringf += '_nho'
  outputname1 = strstart+'_e_'+propertystringf+'_2pop.png'
  if xlabel is None:
    plt.xlabel(propertystring, fontsize=15)
  else:
    plt.xlabel(xlabel, fontsize=15)
  plt.ylabel('PDF', fontsize=15)
  plt.yticks([])
  plt.xlim((xlims[0], xlims[1]))
  plt.tick_params(axis='both', which='major', labelsize=12)
  #if legendon:
  #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand", fontsize=12)
  plt.legend(loc=2)
  plt.savefig('figs4/'+outputname1, bbox_inches='tight')
  plt.close()

  lowedat = np.sort(lowedat)
  highedat = np.sort(highedat)

  d, pval = stats.ks_2samp(lowedat, highedat)
  lowcum = np.array(range(len(lowedat)))/float(len(lowedat)-1)
  highcum = np.array(range(len(highedat)))/float(len(highedat)-1)

  plt.figure(3)
  plt.plot(lowedat, lowcum, color='b', alpha=0.8, label='Low e') 
  plt.plot(highedat, highcum, color='r', alpha=0.8, label='High e') 
  texty = 0.9
  textx = (max(plist)-min(plist))/20. + min(plist)
  plt.text(textx, texty, 'KS p-value = '+str(radvel.utils.round_sig(pval)), fontsize=14, horizontalalignment='left')
  outputname1 = strstart+'_e_'+propertystringf+'_cdf_2pop.png'
  if xlabel is None:
    plt.xlabel(propertystring, fontsize=15)
  else:
    plt.xlabel(xlabel, fontsize=15)
  if log:
    #propertystring = 'Log10['+propertystring+']'
    locs, labels = plt.xticks()
    plt.xticks(locs, [radvel.utils.round_sig(v) for v in [10.**loc for loc in locs]])
  plt.ylabel('CDF', fontsize=15)
  plt.xlim((xlims[0], xlims[1]))
  plt.ylim((-0.05, 1.05))
  plt.yticks([])
  plt.tick_params(axis='both', which='major', labelsize=16)
  if legendon:
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand", fontsize=12)
  plt.savefig('figs4/'+outputname1, bbox_inches='tight')
  plt.close()

  print "plotted "+propertystring

  return 



