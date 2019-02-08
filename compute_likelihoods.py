## Script to process the Kepler duration data to determine the likelihood of different eccentricity distributions
##

import numpy as np
import random
import scipy
from scipy import stats
from scipy.stats import norm
import os


loopvec=range(100)
#loopvec=[i+54 for i in range(47)]

for loopi in loopvec:


  print "Drawing random numbers. This may take a bit"
  size=9
  #base=10
  base=6
  mult=2
  randomuscounter = 0L
  randomus = np.random.uniform(low=0, high=1, size=(base**size)*mult)
  modrandomu = base**size-1
  print "us done"
  normmult = 10
  randomnormscounter = 0L
  randomnorms = np.random.randn(normmult * base**size *mult)
  modrandomn = normmult * base**size-1
  print "norms done"
  randomraycounter = 0L
  randomray = np.random.rayleigh(scale=1, size=(base**size)*mult)
  modrandomr = base**size-1
  print "All drawn"
  ## speed things up
  
  
  
  # This assumes we know rstar, mstar in the synthetic draw
  #   then adds those uncertainties into the duration uncertainty
  #   This is a basically just a trick to make the computation of the integral via monte carlo
  #   more efficient. 
  NoStellarSigma=True
  
  
  # Load data in
  cks_fname = 'resources/CKS_GAIA_Stellar_Parameters.txt'
  cks_data = np.loadtxt(cks_fname, usecols=[0,1,2,3,4,5,6])#, usecols=[0,11,12,13,14,15,16])
  cks_data = np.transpose(cks_data)
  koi_fname = "resources/KOI_list_SNR_CONFIRMED.tsv" 
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
  
  
  
  # Returns relavent properties of all stars in CKS sample
  def get_CKS_table():
    data = cks_data
  
    kois = data[0]
    mstar = data[1]
    mstar_e_lo = data[2]
    mstar_e_hi = data[3]
    rstar = data[4]
    rstar_e_lo = data[5]
    rstar_e_hi = data[6]
  
    return kois, mstar, mstar_e_lo, mstar_e_hi, rstar, rstar_e_lo, rstar_e_hi
  
  
  
  # Returns relevant planet properties for all stars in Kepler KOI catalog
  def get_kepler_transit_durations():
    data = koi_data
  
    kois = data[0]
    periods = data[1]
    durations = data[2]/24.
    rpors = data[3]
    snr = data[4]
    return kois, periods, durations, rpors, snr
  
  
  
  # Returns planet property uncertainties for all planets in the KOI catalog
  def get_kepler_transit_errors():
    errors = koi_errs
    kois = errors[0]
    period_e_hi = errors[1]
    period_e_lo = errors[2]
    durations_e_hi = errors[3]/24.
    durations_e_lo = errors[4]/24.
    rpors_e_hi = errors[5]
    rpors_e_lo = errors[6]
  
    return kois, period_e_hi, period_e_lo, durations_e_hi, durations_e_lo, rpors_e_hi, rpors_e_lo
  
  
  
  # Returns a sample of the CKS stellar properties for a given KOI
  def sample_stellar_properties(koi):
    kois, mstar, mstar_e_lo, mstar_e_hi, rstar, rstar_e_lo, rstar_e_hi = get_CKS_table()
  
    koii = np.where(koi == kois)
    nmatch = len(koii[0])
    if (nmatch != 1):
      print "sample star"
      print koii
      print kois[koii]
      print nmatch
      raise Exception(" ")
  
    thismstar = mstar[koii]
    thismstar_e_lo = mstar_e_lo[koii]
    thismstar_e_hi = mstar_e_hi[koii]
    thisrstar = rstar[koii]
    thisrstar_e_lo = rstar_e_lo[koii]
    thisrstar_e_hi = rstar_e_hi[koii]
  
    if NoStellarSigma:
      mstardraw = thismstar[0]
      rstardraw = thisrstar[0]
    else:
      mstardraw = assymetric_draw(thismstar[0], thismstar_e_lo[0], thismstar_e_hi[0])
      rstardraw = assymetric_draw(thisrstar[0], thisrstar_e_lo[0], thisrstar_e_hi[0])
      # make sure the radii and masses are > 0 (only a problem if very uncertain)
      mstardraw = max(0., mstardraw)
      rstardraw = max(0., rstardraw)
  
    return mstardraw, rstardraw
  
  
  
  ### Draws a sample from a 2-sided Gaussian distribution 
  ##def assymetric_draw(median, sigmalow, sigmahigh):
  ##  sign = 1 if (randomus[randomuscounter] < 0.5) else -1
  ##  randomuscounter+=1
  ##  if sign==1:
  ##    gaussdraw = np.abs(np.random.randn()*sigmahigh)
  ##  else:
  ##    gaussdraw = np.abs(np.random.randn()*sigmalow)
  ##  final = median + sign*gaussdraw 
  ##
  ##  return final
  
  
  
  # Draws an e value from a Rayleigh distribution with sigma=width
  #   It uses aors to discard physically impossible values where the planet
  #   goes inside the star
  def draw_e(width, aors):
    ooaors = 1./aors
    global randomraycounter
    # why doesn't python have a do-while loop? 
    count = 0
    while True:
      e = randomray[randomraycounter]*width
      randomraycounter+=1
      #randomraycounter = randomraycounter % modrandomr 
      # make sure that the planet doesn't go inside the star
      if ((1.-e) > ooaors) and (e < 1.0):
        break
      count+=1
      # we will redraw stellar properties if a valid e can't be found
      #   this is mostly to speed up the code in edge cases
      if count>100:
        return e, 0
    return e, 1
  
  
  
  # Applies Kepler's law to get a/Rstar
  def get_aors(period, rho_star):
    # [universal gravitational constant]^(1/3) * (day)^(2/3) * [solar density]^(1/3)* ( 1/(3*pi) )^(1/3)
    aorsfactor = 4.206
    aors = aorsfactor * period**(2./3.) * (rho_star*(4./3.*np.pi))**(1./3.)
    if (aors < 1):
      aors = np.nan
    return aors
  def get_aors_sigma(period, rho_star, e_period, e_rho_star):
    # [universal gravitational constant]^(1/3) * (day)^(2/3) * [solar density]^(1/3)* ( 1/(3*pi) )^(1/3)
    aorsfactor = 4.206
    s1 = 2./3.*aorsfactor*period**(-1./3.)*(rho_star*4./3.*np.pi)**(1./3.)*e_period
    s2 = 1./3.*aorsfactor*period**(2./3.)* (4./3.*np.pi)**(1./3.) *(rho_star)**(-2./3.)*e_rho_star 
    e_aors = np.sqrt(s1**2 + s2**2)
    return e_aors
  
  
  
  # Creates a synthetic transit duration for a given koi by drawing i and e 
  #   and rpors, etc. from the Kepler posteriors
  def draw_transit_duration(period, rpors, koi, esigma, iprev, isigma, snr, tobs, sigobs, aors, cosinom, cosinom_e):
    global randomuscounter
    global randomnormscounter
    e, status = draw_e(esigma, aors)
    if not status:
      print "Can not find a valid e"
      print koi, period
      raise
  
    omega = randomus[randomuscounter]*2*np.pi
    randomuscounter+=1
    #randomuscounter = randomuscounter % modrandomu 
    rho_c = (1.-e**2)/(1.+e*np.sin(omega))
    #maxcosi = min(1.0, (1.+rpors)/(aors*rho_c))
    maxcosi = (1.+rpors)/(aors*rho_c)
    sigfrac = sigobs/tobs
  
    # we will loop over cosi and e until we find suitable values only for the first planet
    #   in the system. Otherwise we must go with whatever is given to us. Otherwise, we seek out certain
    #   e and omega for a given inclination for a given system. 
    firstpl = False
    if iprev < 0:
      firstpl = True
    if firstpl:
      counthere=0
      while firstpl:
        counthere+=1
        if counthere % 1000 == 0:
          print "counthere=%i" % counthere
          #print koi, tdur, tobs, draw_snr, period, aors, rpors, rho_c 
          print koi, tobs, period, aors, rpors, rho_c, e, omega, cosi, cosinom_e, cosinom 
        ## Importance Sampling
        ## iprev = randomnorms[randomnormscounter]*cosinom_e + cosinom
        ## randomnormscounter+=1
        ## randomnormscounter = randomnormscounter % modrandomn 
        ### if the first planet doesn't transit, then redraw inclination plane and inclination
        ##if cosi > maxcosi:
        ##  continue
        iprev = randomus[randomuscounter]*maxcosi
        randomuscounter+=1
        #randomuscounter = randomuscounter % modrandomu 
        cosi = iprev
        tdur = transit_duration(period, rho_c, e, aors, cosi, rpors) 
        # equation S11
        tdur *= (1. + randomnorms[randomnormscounter]*sigfrac)
        randomnormscounter+=1
        #randomnormscounter = randomnormscounter % modrandomn 
        draw_snr =  snr * np.sqrt(tdur/tobs)
        # if it passes the SNR test move on 
        if draw_snr > 7.1:
          break
        # otherwise, need to redraw e, omega, and restart, redrawing i. 
        e, status = draw_e(esigma, aors)
        if not status:
          print "Can not find a valid e 2"
          print koi, period
          raise
        omega = randomus[randomuscounter]*2*np.pi
        randomuscounter+=1
        #randomuscounter = randomuscounter % modrandomu 
        rho_c = (1.-e**2)/(1.+e*np.sin(omega))
        #maxcosi = min(1.0, (1.+rpors)/(aors*rho_c))
        maxcosi = (1.+rpors)/(aors*rho_c)
    else:
      cosi = iprev+randomnorms[randomnormscounter]*isigma
      randomnormscounter+=1
      #randomnormscounter = randomnormscounter % modrandomn 
      tdur = transit_duration(period, rho_c, e, aors, cosi, rpors)
      # equation S11
      tdur *= (1. + randomnorms[randomnormscounter]*sigfrac)
      randomnormscounter+=1
      #randomnormscounter = randomnormscounter % modrandomn 
      draw_snr =  snr * np.sqrt(tdur/tobs)
    ## Importance Sampling
    ###This is wrong for multiplanet system
    ##normalpdf = norm.pdf(iprev, cosinom, cosinom_e)
    ##adjustedpdf = normalpdf / (norm.cdf(maxcosi, cosinom, cosinom_e) - norm.cdf(0, cosinom, cosinom_e))
    ##weight = 1./maxcosi / adjustedpdf
    weight=1.
  
    return tdur, iprev, draw_snr, weight
    
  
  
  # Computes transit duration for a set of stellar and planet parameters
  def transit_duration(P, rho_c, e, aors, cosi, rpors):
    #sini = np.sin(np.arccos(cosi))
    sini = np.sqrt(1.-cosi**2) # the same but probably faster 
  
    prefac = P/np.pi * rho_c**2 / np.sqrt(1.-e**2)
    argument = (1.+rpors)**2 - aors**2 * rho_c**2 * cosi**2
    # Check if duration is zero (this happens in the multi-planet cases sometimes)
    if (argument <= 0.):
      return 0.
    numerator = np.sqrt(argument)
    denominator = aors * rho_c * sini
    # this min is in there due to very close to star transits sometimes have value slightly above 1.0
    #T = prefac * np.arcsin( min(numerator / denominator, 1.0))
    T = prefac * np.arcsin( numerator / denominator )
  
    return T
  
  
  
  # Finds and symmetrizes the uncertainties for all relevant parameters in a KOI list
  def get_symmetric_uncertainties_list(koi_decimal_list):
    lenlist = len(koi_decimal_list)
    period_e_list=np.empty(lenlist)
    durations_e_list=np.empty(lenlist)
    rpors_e_list=np.empty(lenlist)
    rhostar_e_list=np.empty(lenlist)
    for i in range(lenlist):
      period_e, durations_e, rpors_e, rhostar_e = get_symmetric_uncertainties(koi_decimal_list[i])
      period_e_list[i] = period_e
      durations_e_list[i] = durations_e
      rpors_e_list[i] = rpors_e
      rhostar_e_list[i] = rhostar_e
    return period_e_list, durations_e_list, rpors_e_list, rhostar_e_list
  
  
  # Finds and symmetrizes the uncertainties for all relevant parameters for a KOI
  #   Both from the CKS data and the Kepler candidate catalog
  def get_symmetric_uncertainties(koi_decimal):
    kois, period_e_hi, period_e_lo, durations_e_hi, durations_e_lo, rpors_e_hi, rpors_e_lo = get_kepler_transit_errors()
  
    koii = np.where(koi_decimal == kois)
    nmatch = len(koii[0])
    if nmatch != 1:
      print "ERROR"
      print nmatch
      print koi_decimal
      print kois[0:30]
      raise
    period_e = np.mean([abs(period_e_hi[koii]), abs(period_e_lo[koii])])
    durations_e = np.mean([abs(durations_e_lo[koii]), abs(durations_e_hi[koii])]) 
    rpors_e = np.mean([abs(rpors_e_hi[koii]), abs(rpors_e_lo[koii])])
  
    kois, mstar, mstar_e_lo, mstar_e_hi, rstar, rstar_e_lo, rstar_e_hi = get_CKS_table()
    koii = np.where(np.floor(koi_decimal) == kois)
    nmatch = len(koii[0])
    if nmatch != 1:
      raise
    mstar_e = np.mean([mstar_e_lo[koii], mstar_e_hi[koii]])
    rstar_e = np.mean([rstar_e_lo[koii], rstar_e_hi[koii]])
    mstari = mstar[koii][0]
    rstari = rstar[koii][0]
    s1 = mstar_e / (4./3.*np.pi*(rstari**3))
    s2 = rstar_e*(-3.)*mstari / (4./3.*np.pi*(rstari**4))
    rhostar_e = np.sqrt(s1**2 + s2**2)
  
    return period_e, durations_e, rpors_e, rhostar_e
  
  
  
  
  # Computes and returns a synthetic duration draw and the nominal (b=0, e=0) transit duration
  #   for the list of inputted KOIs
  # NOTE: the koilist MUST be sorted.
  def generate_synthetic_distribution(koilist, esigma, isigma, kep_kois, kep_periods, kep_rpors, kep_snrs, kep_durations, \
      duration_e, rhostar_e, rpors_e, aors, rhostar, dur0s, cosinom, cosinom_e):
    global randomuscounter
    global randomnormscounter
  
    nplanets = len(kep_kois)
  ## can edit this
    nkois = len(koilist)
    synthetic_durations = np.empty(nplanets)
    nominal_durations = np.empty(nplanets)
    snrlist = np.empty(nplanets)
    weightlist = np.empty(nplanets)
    # iprev is a given system's inclination plane around the same star. If its the first planet we've looked at
    #   in the system we set iprev=-1 and draw it uniformly from cos(i)
    iprev = -1
    plinsys = 0
    koi_i = 0
    #snrlist=[]
    #synlist=[]
    #nomlist=[]
    # loop through list of planets in the list
    i=0
    counttrys=1
    while i <= nplanets:
      if counttrys % 10000 == 0 and plinsys==1:
        print "retry system %s %i time so far, iplane=%f" % (koilist[koi_i], counttrys, iprev*180/np.pi)
      if (i==nplanets) or (koilist[koi_i] != np.floor(kep_kois[i])):
        # Check that the snr of all of the planets in the previous system are detectable
        if all(snrlist[i-plinsys:i]): #and (plinsys > 0):
          if counttrys > 100:
            print "had to retry system %s %i times, iplane=%f" % (koilist[koi_i], counttrys, iprev*180/np.pi)
          # append the good values to the list
          #[synthetic_durations.append(di) for di in synlist]
          #[nominal_durations.append(d0) for d0 in nomlist]
          if (i==nplanets):
            break
          # advance the stellar properties list index
          koi_i += 1
          # and set the system's overall inclination plane to uninitialized 
          iprev = -1
          # and reset the number of planets in this system 
          plinsys = 0
          # reset the list of snrs
          #snrlist = []
          #synlist = []
          #nomlist = []
          counttrys=1
        else:
          print "This condition should never be reached now that I shortcircuit earlier"
          print  kep_kois[i], kep_periods[i] 
          # retry this whole system
          # go back to first planet
          i -= plinsys
          # reset to first planet
          plinsys = 0
          # draw new inclination plane
          iprev = -1
          # reset the list of snrs
          #snrlist = []
          #synlist = []
          #nomlist = []
          counttrys+=1
      # if it's just another planet around the same star, draw its duration as usual 
      else:
        pass
     
      koii = koilist[koi_i]
      duration_ei, rpors_ei, rhostar_ei = duration_e[i], rpors_e[i], rhostar_e[i]
      # This needs to be inside while loop if stellarsigma=true
      if NoStellarSigma == False:
        print "a/Rstar calculation must be moved inside the loop b/c it depends on the stellar properties"
        print "  which should be resampled for each draw (but only once per system!)"
        raise 
      aorsi = aors[i]
      # Draw a transit duration with a given e, i, distribution, and other planet properties
      duri, iprev, snri, weight = draw_transit_duration(kep_periods[i], kep_rpors[i], koii,  esigma, iprev, isigma, kep_snrs[i], kep_durations[i], duration_ei, aorsi, cosinom[i], cosinom_e[i])
      # If any planet fails, must restart whole system
      if (not np.isfinite(snri)) or (snri < 7.1):
        # go back to first planet
        i -= plinsys
        plinsys = 0 
        # draw new inclination plane
        iprev = -1
        # reset lists of draws
        #snrlist = []
        #synlist = []
        #nomlist = []
        counttrys += 1
        continue
  
      snrlist[i] = snri >= 7.1
      # Eq S12
      ### can do this outside of loop for speedup
      dur0i = dur0s[i] 
      sigratio = np.sqrt( (rhostar_ei/(3*rhostar[koi_i]))**2 + (rpors_ei/(1.+kep_rpors[i]))**2 ) 
      dur0i *= (1. + randomnorms[randomnormscounter]*sigratio)
      randomnormscounter+=1
      #randomnormscounter = randomnormscounter % modrandomn 
  
      #print duri, dur0i, kep_periods[i], kep_rpors[i], koii 
  
      synthetic_durations[i] = duri
      nominal_durations[i] = dur0i
      weightlist[i] = weight
      #synlist.append(duri)
      #nomlist.append(dur0i)
  
      plinsys += 1
      i+=1
  
    #synthetic_durations=np.array(synthetic_durations)
    #nominal_durations=np.array(nominal_durations)
  ####
  #  for i in range(len(nominal_durations)):
  #    #if not np.isfinite(synthetic_durations[i]) or not np.isfinite(nominal_durations[i]):
  #    print i, synthetic_durations[i], nominal_durations[i], kep_kois[i]
  
    return synthetic_durations, nominal_durations, weightlist
  
  
  
  # get aors for a list of KOIs 
  def get_aors_list(koi_decimal_list, period_list):
    aors_list=[]
    rho_star_list=[]
    for i in range(len(koi_decimal_list)):
      mstari, rstari = sample_stellar_properties(np.floor(koi_decimal_list[i]))
      rho_stari = mstari / (4./3.*np.pi*rstari**3.)
      aorsi = get_aors(period_list[i], rho_stari)
      aors_list.append(aorsi)
      rho_star_list.append(rho_stari)
    return np.array(aors_list), np.array(rho_star_list)
  
  
  ###
  def get_aors_err_list(koi_decimal_list, period_list):
    aors_list=[]
    period_e, duration_e, rpors_e, rhostar_e = get_symmetric_uncertainties_list(koi_decimal_list)
    for i in range(len(koi_decimal_list)):
      mstari, rstari = sample_stellar_properties(np.floor(koi_decimal_list[i]))
      rho_stari = mstari / (4./3.*np.pi*rstari**3.)
      aorsi = get_aors_sigma(period_list[i], rho_stari, period_e[i], rhostar_e[i])
      aors_list.append(aorsi)
    return np.array(aors_list)
  
  
  
  # Computes nominal duration (b=0,e=0) for list of planet properties
  def compute_circular_edgeon_duration_list(kep_periods, kep_rpors, aors):
    npl = len(kep_periods)
    dur0_list=np.empty(npl)
    for i in range(npl):
      dur0i = compute_circular_edgeon_duration(kep_periods[i], kep_rpors[i], aors[i])
      dur0_list[i] = dur0i
    return np.array(dur0_list)
  
  
  
  # Computes the nominal duration (b=0,e=0) for a planet with the given properties
  def compute_circular_edgeon_duration(period, rpors, aors):
    prefactor = 1./np.pi
    T0 = prefactor * period * np.arcsin( (1.+rpors)/aors ) 
    return T0
  
  
  
  # Function to compute the likelihood for each koi in koilist
  #   A single draw from the i, e, rstar, and mstar posteriors is made
  def compute_duration_likelihood(koilist, esigma, isigma, kep_kois, kep_periods, kep_rpors, kep_snrs, kep_durations, \
        duration_e, rhostar_e, rpors_e, aors, rhostar, dur0s, kep_durations_use_e, cosinom, cosinom_e):
    synth_durations, nom_durations, weights = \
        generate_synthetic_distribution(koilist, esigma, isigma, kep_kois, kep_periods, kep_rpors, kep_snrs, kep_durations, \
        duration_e, rhostar_e, rpors_e, aors, rhostar, dur0s, cosinom, cosinom_e)
    #np.savetxt('durdraws/draws3.txt', np.transpose([synth_durations, nom_durations]))
    #exit()
    numerator = (kep_durations/nom_durations - synth_durations/nom_durations)**2
    denominator = 2 * kep_durations_use_e**2
    chi_sq = numerator/denominator
    likelihood = np.exp(-numerator/denominator)
    #####
    #for i in range(len(likelihood)):
    #  print i, kep_kois[i], likelihood[i] 
    return likelihood, weights
  
  
  
  # Function to take N draws from the i, e, rstar, and mstar priors and find the 
  #   integral in Section 4 of the Supplement to Xie+2017
  #   The mean, axis=0 line weights the draws appropriately (i.e. relative to their frequency)
  #   and the total log-likelihood for all KOIs in koilist is combined and returned
  def compute_duration_likelihood_integral(koilist, esigma, isigma, N, txtstr):
    kep_kois, kep_periods, kep_durations, kep_rpors, kep_snrs = get_kepler_transit_durations()
    thesekois = np.nonzero(np.in1d(np.floor(kep_kois), koilist))  
    kep_kois, kep_periods, kep_durations, kep_rpors, kep_snrs = \
        kep_kois[thesekois], kep_periods[thesekois], kep_durations[thesekois], kep_rpors[thesekois], \
        kep_snrs[thesekois]
    kois, mstar, mstar_e_lo, mstar_e_hi, rstar, rstar_e_lo, rstar_e_hi = get_CKS_table()
    thesekois = np.nonzero(np.in1d(kois, koilist))  
    kois, mstar, mstar_e_lo, mstar_e_hi, rstar, rstar_e_lo, rstar_e_hi = \
        kois[thesekois], mstar[thesekois], mstar_e_lo[thesekois], \
        mstar_e_hi[thesekois], rstar[thesekois], rstar_e_lo[thesekois], rstar_e_hi[thesekois]
    
    period_e, duration_e, rpors_e, rhostar_e = get_symmetric_uncertainties_list(kep_kois)
    if NoStellarSigma==False:
      raise
    aors, rhostar = get_aors_list(kep_kois, kep_periods) 
    aors_e = get_aors_err_list(kep_kois, kep_periods) 
    dur0s = compute_circular_edgeon_duration_list(kep_periods, kep_rpors, aors)
    kep_durations_use_e = get_sigma_duration_list(kep_kois)
    # This is approximate
    #cosinom = np.sqrt(1. - (kep_durations/dur0s)**2 ) / aors
    cosinom = np.sqrt( (1./aors**2)*(1+kep_rpors)**2 - np.sin(np.pi*kep_durations/kep_periods)**2 )
    maxcosi = np.array([min(1.0, (1.+kep_rpors[i])/(aors[i])) for i in range(len(aors))])
    cosinom = np.nan_to_num(cosinom)
    cosinom = np.array([min(maxcosi[i], cosinom[i]) for i in range(len(cosinom))])
    #sigma_cosinom = kep_durations_use_e / dur0s / (aors*cosinom) 
    #sigma_cosinom = (kep_durations_use_e / (kep_durations/dur0s)) 
    sigma_cosinom = np.sqrt( (1./cosinom * -2 / aors**3 * aors_e)**2 ) 
    minsig=0.003
    sigma_cosinom = np.array([max(min(maxcosi[i]/2, sigma_cosinom[i]), minsig) for i in range(len(sigma_cosinom))])*3*(1. + esigma/0.05)
    likelihoods=[]
    weights=[]
    for i in range(N):
      if i % 10 == 0:
        print "n = %i" % i
      triali, weightsi = compute_duration_likelihood(koilist, esigma, isigma, kep_kois, kep_periods, kep_rpors, kep_snrs, kep_durations, \
          duration_e, rhostar_e, rpors_e, aors, rhostar, dur0s, kep_durations_use_e, cosinom, sigma_cosinom)
      likelihoods.append(triali)
      weights.append(weightsi)
    likelihoods = np.asarray(likelihoods)
    weights = np.asarray(weights)
    #likelihoods = np.mean(likelihoods, axis=0)
    print cosinom[2], sigma_cosinom[2] 
    #print likelihoods[:,2]
    i=0
    for wl in zip(weights[:,2], likelihoods[:,2]):
      print i, wl
      i+=1
    i=0
    for wl in zip(weights[:,0], likelihoods[:,0]):
      print i, wl
      i+=1
    likelihoods = 1./N*np.sum(likelihoods*weights, axis=0)
    endstr ='_%f_%f_%i.txt' % (esigma, isigma, N)
    np.savetxt('data3/'+txtstr+endstr, np.transpose([kep_kois, likelihoods]))
    print randomuscounter
    print randomnormscounter
    print randomraycounter
    total_loglike = np.sum(np.log(likelihoods))
    return total_loglike
  
  
  
  # Returns a list of the duration ratio uncertainties for each koi in koilist
  def get_sigma_duration_list(koi_decimal_list):
    siglist = []
    for i in range(len(koi_decimal_list)):
      siglist.append(get_sigma_duration(koi_decimal_list[i]))
    return np.transpose(np.array(siglist))
  
  
  
  # Computes the uncertainty in the measured duration ratio for a given koi
  def get_sigma_duration(koi_decimal):
    kep_kois, kep_periods, kep_durations, kep_rpors, kep_snrs = get_kepler_transit_durations()
    kois, mstar, mstar_e_lo, mstar_e_hi, rstar, rstar_e_lo, rstar_e_hi = get_CKS_table()
  
    kepkoii = np.where(koi_decimal == kep_kois)
    ckskoii = np.where(np.floor(koi_decimal) == kois)
  
    duration = kep_durations[kepkoii][0]
    rpors = kep_rpors[kepkoii][0]
    mstar = mstar[ckskoii][0]
    rstar = rstar[ckskoii][0]
    if NoStellarSigma==False:
      raise
    rhostar = mstar / ( (4./3.*np.pi)* rstar**3)
  
    period_e, duration_e, rpors_e, rhostar_e = get_symmetric_uncertainties(koi_decimal)
    s1 = duration_e / duration 
    s2 = rhostar_e / (3.*rhostar)
    s3 = rpors_e / (1.+rpors)
    sigma = np.sqrt(s1**2 + s2**2 + s3**2)
    return sigma
  
  
  
  ## Computes the distribution of p-values for N realizations of the koilist's synthetic durations
  ##   for a given sigma_e, compared to the measured durations 
  #### This is no longer used
  #def get_pval_distribution(koilist, esigma, N):
  #  pdist = np.array([compare_duration_distributions_ks(koilist, esigma) for i in range(N)])
  #  mn = np.mean(pdist)
  #  sig = np.std(pdist)
  #
  #  return mn, sig
  
  
  
  ## Grid over a bunch of e values, comparing duration ratio distributions via KS-test
  #### This is no longer used
  #def egrid(koilist, estep, nperstep):
  #  nsteps = int((1.-0.)/estep)
  #  evals = np.arange(nsteps, dtype='f')/nsteps
  #  pvals=[]
  #  for i in range(nsteps):
  #    evali = evals[i]
  #    mni, sigi = get_pval_distribution(koilist, evali, nperstep) 
  #    pvals.append(mni)
  #  return pvals
  
  
  
  # Grid over a bunch of e values (separated by estep), and compute the likelihood integrals
  #   with N draws from P(TDUR), for a given KOI list
  #   The e-grid and log-likelihoods for each e are returned
  def egrid_like(koilist, emax, estep, imax, istep, nperstep, savefile=None):
    if savefile != None:
      if os.path.exists('data3/'+savefile):
        os.rename('data3/'+savefile, 'data3/'+savefile+'.backup')
        print "WARNING: moved existing file to %s.backup" % 'data3/'+savefile
      with open('data3/'+savefile, 'w') as svf: 
        svf.write("loglike\t\t\te\t\t\t\ti\n")
    if isinstance(estep, list):
      n_esteps = len(estep)
      evals = estep
    else:
      maxe = emax
      n_esteps = int((maxe-0.)/estep)
      evals = np.arange(n_esteps, dtype='f')*maxe/n_esteps
    if isinstance(istep, list):
      n_isteps = len(istep)
      ivals = istep
    else:
      maxi = imax
      n_isteps = int((maxi-0.0)/istep)
      ivals = np.arange(n_isteps, dtype='f')*maxi/n_isteps
    likes = np.zeros((n_esteps, n_isteps), dtype='f')
    for k in range(n_isteps):
      ivalk = ivals[k]
      print "i = %f" % ivalk
      for j in range(n_esteps):
        evalj = evals[j]
        print "e = %f" % evalj
        evalj = evals[j]
        like = compute_duration_likelihood_integral(koilist, evalj, ivalk, nperstep, savefile) 
        print like
        likes[j,k] = like
        if savefile != None:
          with open('data3/'+savefile, 'ab') as svf:
            np.savetxt(svf, np.asarray([[like, evalj, ivalk]]))
    return evals, ivals, likes
  
  
  
  ## select stars with a single planet or multiple planets
  kep_kois, kep_periods, kep_durations, kep_rpors, kep_snrs = get_kepler_transit_durations()
  
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
  
  ## Run the fit 
  #evals, ivals, llike = egrid_like(cks_multis, 0.3, [i/100.+0.0+0.000 for i in range(10)],\
  #    0.12, [0.06], 10000, savefile='multis_m12_10000_1_7_' +str(loopi)+'.txt') 
  evals, ivals, llike = egrid_like(cks_singles, 0.3, [i/100.+0.1+0.000 for i in range(10)],\
      0.12, [0.0], 10000, savefile='singles_m12_10000_2_' +str(loopi)+'.txt') 
  


  del randomray
  del randomus 
  del randomnorms

