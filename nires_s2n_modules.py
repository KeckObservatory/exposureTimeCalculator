#!/usr/bin/env python
# coding: utf-8
 
# In[1]:
 
 
import numpy as np
import scipy.constants as sc
from astropy.io import ascii
import math
 
SourceTemplate = '/home/filer/elucas/exposureTimeCalculator/ascii/cdbs'
InstrThroughputPath = '/home/filer/elucas/exposureTimeCalculator/ascii/'
IRSkyPath = '/home/filer/elucas/exposureTimeCalculator/ascii/'
 
import os
os.environ['PYSYN_CDBS']
SourceTemplate
 
import pysynphot as S
import matplotlib.pyplot as plt
 
 
# In[3]:
 
 
#plt.ion()
#fig = plt.figure()
#fig.canvas.draw()
 
mag_src = 18.2
src_type = "sy1"
redshift = 2
teff = 5000
exp_time = 3600 #sec
coadds = 1
dither = "AB"
dither_repeat = 1
obs_wave = 2.2 #um
seeing = 0.8 #arcsec
num_reads = 16 

qso = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'qso_template.fits')
qso_sp = S.FileSpectrum(qso).redshift(redshift)
sy1 = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'seyfert1_template.fits')
sy1_sp = S.FileSpectrum(sy1).redshift(redshift)
sy2 = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'seyfert2_template.fits')
sy2_sp = S.FileSpectrum(sy2).redshift(redshift)

## fixed instrument parameters
slit_width = 0.55 #arcsec
slit_lenth = 18.0 #arcsec


# In[13]:


def calculate_nires_s2n_pointsource(mag_src, obs_wave, exp_time, coadds, dither, dither_repeat, seeing, num_reads):
    
    dark_current = 0.01 #e-/s
    throughput = ascii.read(InstrThroughputPath+"nires_throughput_18_edit.dat")
    tpt_wave = throughput[0][:]
    tpt_val = throughput[1][:]
    tpt_val_interp = np.interp(obs_wave, tpt_wave, tpt_val)
    gain = 3.8 #e-/ADU
    pix_size = 0.123 #arcsec/pix
    collect_area = 76*1e4 # 76 m^2 = 76e4 cm^2
    del_lmbd = obs_wave/2700/1e4 # cm
    spatial_cov = (0.55*seeing)/(math.pi*seeing**2)
    #num_pix = math.pi*(seeing/2)**2 / pix_size**2
    num_pix_slt = 18*0.55/pix_size**2
    num_pix_src = 0.55*seeing

    ## mauna kea IR sky bkg at airmass 1, water vapor col 1 
    IR_sky = ascii.read(IRSkyPath+"mk_skybg_zm_10_10_ph.dat")
    sky_wave = IR_sky[0][:]/1000 # um
    sky_bkg  = 1000*IR_sky[1][:] #photon/s/arcsec^2/cm/cm^2
    sky_bkg_interp = np.interp(obs_wave, sky_wave, sky_bkg)

    flux_src = np.power(10,-0.4*mag_src)/1e8 #input source flux in erg/s/cm^2/cm

    photon_energy = sc.h*1e7*(sc.c*1e2)/(obs_wave*1e4) #single photon energy in erg

    flux_src_phot = flux_src/photon_energy #source flux in photon/s/cm^2/cm

    sig_src = flux_src_phot*collect_area*tpt_val_interp*del_lmbd*spatial_cov # source signal in e-/s 


    if num_reads == 2:
        read_noise = 15 #CDS
    elif num_reads == 16:
        read_noise = 5 #MCDS

    if dither == "AB":
        num_dith = 2
    elif dither == "ABBA":
        num_dith = 4
        
    
    sig_src_int = sig_src*exp_time #source flux integrated over time, single frame
    noise_int = (sig_src*exp_time + num_pix_slt*dark_current*exp_time + (num_pix_slt-num_pix_src)*sky_bkg_interp*exp_time + num_pix_slt*read_noise**2*exp_time)**(1/2) #noise, single frame

    s2n = sig_src_int / noise_int * math.sqrt(num_dith) # s/n 

    return s2n
   # print("src flux", flux_src)
   # print("photon energy=", photon_energy)    
   # print("signal =", sig_src_int)     
   # print("noise =", noise_int)
   # print("S/N =", s2n)


# In[14]:


def calculate_nires_s2n_qso(src_type, redshift, obs_wave, exp_time, coadds, dither, dither_repeat, seeing, num_reads):
    ## fixed instrument parameters
    slit_width = 0.55 #arcsec
    slit_lenth = 18.0 #arcsec
    dark_current = 0.01 #e-/s
    pix_size = 0.123 #arcsec/pix
    collect_area = 76*1e4 # 76 m^2 = 76e4 cm^2
    obs_wave_angstrom = obs_wave*1e4 #convert to angstrom
    del_lmbd = obs_wave_angstrom/2700 # unit: A
    spatial_cov = 0.55*seeing
    num_pix_slt = 18*0.55/pix_size**2
    num_pix_src = 0.55*seeing

    # throughput wave unit in A 
    #tpt_tab = S.FileBandpass('/Users/syeh/Dropbox/keck/etc/nires/nires_throughput_18_long.dat')
    tpt_tab = S.FileBandpass(InstrThroughputPath+"nires_throughput_18_edit.dat")
    #convert input wavelength cm to A then interp
    tpt_val_interp = np.interp(obs_wave_angstrom, tpt_tab.wave, tpt_tab.throughput) 
    print(tpt_val_interp)

    #sky value unit is phot/s/arcsec^2/cm^2/A
    #sky_am10_pw10 = S.FileBandpass('/Users/syeh/Dropbox/keck/etc/IRsky/IRsky_am10_pw10.dat') 

    #sky_bkg_interp = np.interp(obs_wave, sky_am10_pw10.wave, sky_am10_pw10.throughput) #unit: phot/s/arcsec^2/cm^2/A 
    ## mauna kea IR sky bkg
    IR_sky = ascii.read(IRSkyPath+"IRsky_am10_pw10.dat")
    sky_wave = IR_sky[0][:]/1000 # um
    sky_bkg  = 1000*IR_sky[1][:] #photon/s/arcsec^2/cm/cm^2
    sky_bkg_interp = np.interp(obs_wave_angstrom, sky_wave, sky_bkg)
   
    src_sp = {"qso":qso_sp, "sy1":sy1_sp, "sy2":sy2_sp}
    
    ## convolve source spectra and sky and plot, estimate source s2n 

    ## source signal flux at given obs. wavelength
    sig_src_obs = np.interp(obs_wave_angstrom, src_sp[src_type].wave, src_sp[src_type].flux) #unit: phot/s/cm^2/A 

    sig_src = sig_src_obs*collect_area*tpt_val_interp*del_lmbd*spatial_cov # source signal in e-/s     

    if num_reads == 2:
        read_noise = 15 #CDS
    elif num_reads == 16:
        read_noise = 5 #MCDS
    
    sig_src_int = sig_src*exp_time #source flux integrated over time, single frame
    noise_int = (sig_src*exp_time + num_pix_slt*dark_current*exp_time +              (num_pix_slt-num_pix_src)*sky_bkg_interp*exp_time +              num_pix_slt*read_noise**2*exp_time)**(1/2) #noise, single frame

    s2n = sig_src_int / noise_int # s/n, single frame
 
    #print("signal =", sig_src_int)     
    #print("noise =", noise_int)
    #print("S/N =", s2n)    
    
    return s2n
    


# In[15]:


if src_type == "PointSource":
    s2n=calculate_nires_s2n_pointsource(mag_src, obs_wave, exp_time, coadds, dither, dither_repeat, seeing, num_reads)
    print ("S/N = %f" % s2n)
    
elif src_type == "qso":
    s2n=calculate_nires_s2n_qso(src_type, redshift, obs_wave, exp_time, coadds, dither, dither_repeat, seeing, num_reads)
    print ("S/N = %f" % s2n)
    
elif src_type == "sy1":
    s2n=calculate_nires_s2n_qso(src_type, redshift, obs_wave, exp_time, coadds, dither, dither_repeat, seeing, num_reads)
    print ("S/N = %f" % s2n)
    
elif src_type == "sy2":
    s2n=calculate_nires_s2n_qso(src_type, redshift, obs_wave, exp_time, coadds, dither, dither_repeat, seeing, num_reads)
    print ("S/N = %f" % s2n)
    
#elif src_type == "BlackBody":
#    s2n=nires_s2n_BlackBody(teff, redshift, obs_wave, exp_time, coadds, dither, dither_repeat, seeing, num_reads)
#    print ("S/N = %f" % s2n)


# In[ ]: