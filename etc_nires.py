#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scipy.constants as sc
from astropy.io import ascii
import math

import os
os.environ['PYSYN_CDBS']='/home/filer/elucas/Flask/ETC/ascii/cdbs'
import pysynphot as S
import matplotlib.pyplot as plt
from nires_s2n_modules import *

def do_calc(src_type,mag_src, obs_wave, exp_time, coadds, dither, dither_repeat, seeing, num_reads, redshift):
	qso = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'qso_template.fits')
	qso_sp = S.FileSpectrum(qso).redshift(redshift)
	sy1 = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'seyfert1_template.fits')
	sy1_sp = S.FileSpectrum(sy1).redshift(redshift)
	sy2 = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'agn', 'seyfert2_template.fits')
	sy2_sp = S.FileSpectrum(sy2).redshift(redshift)

	## fixed instrument parameters
	slit_width = 0.55 #arcsec
	slit_lenth = 18.0 #arcsec
	print(src_type)
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

	else:
		s2n = 'N/A'
	return s2n
    
#elif src_type == "BlackBody":
#    s2n=nires_s2n_BlackBody(teff, redshift, obs_wave, exp_time, coadds, dither, dither_repeat, seeing, num_reads)
#    print ("S/N = %f" % s2n)