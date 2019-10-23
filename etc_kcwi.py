import numpy as np
from astropy.io import fits
from bokeh.plotting import figure,output_file,show,save,ColumnDataSource
from bokeh.embed import file_html,components
from bokeh.resources import CDN
from bokeh.models import CrosshairTool,HoverTool
from scipy.interpolate import interp1d
import json
import sys

global datadir
datadir = './datafiles/kcwi/'

def wpA(flux):
	#get wavelength and photons/Angstrom
	w = np.arange(3000.0,6000.0,1.0)
	pA = flux/(2.0*10**-8/w) #for linear f power law
	return w,pA

def obj_cts(w,pA,grating,exptime):
	#get object counts
	A_geo = np.pi/4.0*(1000.0)**2
	eff = get_params(w,grating)
	cts = eff*A_geo*exptime*pA

	return cts

def sky_cts(w,grating,exptime,airmass=None,area=None):
	#get sky counts
	if not airmass:
		airmass = 1.2
	if not area:
		area = 1.0
	A_geo = np.pi/4.*(1000.0)**2
	print("AGEOSKY",A_geo)
	eff = get_params(w,grating)
	cts = eff*A_geo*exptime*sky_mk(w)*airmass*area
	print('EFF',eff[1500])
	print('CTS',cts[1500])
	return cts

def sky_mk(w):
	#load dat file array
	data = np.loadtxt(datadir+'mk_sky.dat',skiprows=14)
	ws = data[:,0]
	fs = data[:,1]
	#open fits file
	file = fits.open(datadir+'lris_esi_skyspec_fnu_uJy.fits')
	header = file[0].header
	f_nu = file[0].data
	#grab data
	dw = header['CDELT1']
	w0 = header['CRVAL1']
	ns = len(fs)
	ws = np.arange(ns)*dw+w0
	f_nu = f_nu[0:len(ws)]
	f_lam = f_nu*10**-29*3.0*10**18/ws**2
	p_lam = f_lam/(2.*10**-8.0/ws)
	psinterpol = interp1d(ws,p_lam,fill_value='extrapolate')
	ps_int = psinterpol(w)
	return ps_int

def get_params(w,grating):
	eff_bl = np.array([0.1825,0.38,0.40,0.46,0.47,0.44])
	eff_bm = np.array([0.1575, 0.33, 0.36, 0.42, 0.48, 0.45])
	eff_bh1 = np.array([0., 0.0, 0.0, 0.0, 0.0, 0.])
	eff_bh2 = np.array([0.,  0.18, 0.3, 0.4, 0.28, 0.])
	eff_bh3 = np.array([0., 0., 0., 0.2, 0.29, 0.31])
	wave_0 = np.array([355.,380.,405.,450.,486.,530.])*10.
	wave_bl = np.array([355., 530.])*10.
	wave_bm = np.array([355., 530.])*10.
	wave_bh1 = np.array([350., 450.])*10.
	wave_bh2 = np.array([405., 486.])*10.
	wave_bh3 = np.array([405., 530.])*10.
	trans_atmtel = np.array([0.54, 0.55, 0.56, 0.56, 0.56, 0.55])

	#get parameters from grating
	if grating == 'BL':
		eff = eff_bl*trans_atmtel
		wave_range = wave_bl
	elif grating == 'BM':
		eff = eff_bm*trans_atmtel
		wave_range = wave_bm
	elif grating == 'BH1':
		eff = eff_bh1*trans_atmtel
		wave_range = wave_bh1
	elif grating == 'BH2':
		eff = eff_bh2*trans_atmtel
		wave_range = wave_bh2
	elif grating == 'BH3':
		eff = eff_bh3*trans_atmtel
		wave_range = wave_bh3
	else:
		print('*** ERROR ***')
		print('INVALID GRATING')
		return 0.0

	#interpolate
	intrpf = interp1d(wave_0,eff,kind='linear',fill_value='extrapolate')
	eff_int = intrpf(w)
	print("EFFINT",eff_int[256])
	print(wave_range[0],wave_range[1])
	indexlt = np.transpose(np.array(np.where(w < wave_range[0])))
	indexgt = np.transpose(np.array(np.where(w > wave_range[1])))
	print(indexlt.shape,indexgt.shape)
	index = np.concatenate((indexlt,indexgt))
	print(index[0])
	if index[0] >= 0:
		eff_int[index] = 0.0
	print(eff_int[1500])
	return eff_int


def do_calc(slicer,grating,gratwave,seeing,exptime,ccdbin,
		    magAB=None,flux=None,nframes=None,ccdspeed=None,
		    spat_bin=None,spec_bin=None,surf_bright=None,
		    emline_w=None):
	
	#establish binning
	bin_factor = 1.0
	if ccdbin == '2x2':
		bin_factor = 0.25
		if slicer == 'Small':
			print('*** WARNING ***')
			print('DO NOT USE 2x2 BINNING WITH SMALL SLICER')
	read_noise = 2.7 #electrons
	if nframes == 'None' or nframes == None:
		nframes = 1.0
	chsz = 3
	seeing1,seeing2 = seeing,seeing
	pixels_per_arcsec = 1.0/0.147
	print(slicer,grating,gratwave)

	#define attributes for specific slicer size
	if slicer == 'Large':
		seeing2 = 1.38
		snr_spatial_bin = seeing1*seeing2
		pixels_spectral = 8
		arcsec_per_slice = 1.35
	elif slicer == 'Medium':
		seeing2 = np.maximum([0.69,seeing])
		snr_spatial_bin = seeing1*seeing2
		pixels_spectral = 4
		arcsec_per_slice = 0.69
	elif slicer == 'Small':
		seeing2 = seeing
		snr_spatial_bin = seeing1*seeing2
		pixels_spectral = 2
		arcsec_per_slice = 0.35
	nslices = seeing/arcsec_per_slice

	#if spatial bin keyword set
	print('spatbin',type(spat_bin[0]))
	if not isinstance(spat_bin[0],str) and not isinstance(spat_bin[1],str):
		nslices = spat_bin[1]/arcsec_per_slice
		snr_spatial_bin = spat_bin[0]*spat_bin[1]
	pixels_spat_bin = pixels_per_arcsec*nslices

	#define angstroms per pixel for specific grating
	if grating == 'BL':
		A_per_pixel = 0.625
	elif grating == 'BM':
		A_per_pixel = 0.28
	elif grating in ['BH1','BH2','BH3']:
		A_per_pixel = 0.125
	#find angstroms per spectral bin
	A_per_specbin = pixels_spectral*A_per_pixel
	print("A",pixels_spectral,A_per_pixel)
	if spec_bin != 'None' and spec_bin != None:
		snr_spectral_bin = spec_bin
	else:
		snr_spectral_bin = A_per_specbin
	pixels_per_snr_specbin = snr_spectral_bin/A_per_pixel

	#tweak flux if emission line width set
	flux1 = 0
	if flux != 'None':
		print("FLUX SET")
		flux1 = flux
		if emline_w != 'None' and emline_w != None:
			flux1 = flux/emline_w
		print("FLUX1",flux1)
	else:
		print("CANNOT SET FLUX")
		if emline_w != 'None' and emline_w != None:
			print('*** WARNING ***')
			print('DO NOT USE MAG_AB FOR EMISSION LINE')
			return
	if magAB != 'None':
		print("MAGAB SPECIFIED,SETTING FLUX")
		print("gratwave",gratwave)
		print("magAB",magAB)
		print("3.0e18/gratwave**2",3.0*10**18/gratwave**2)
		flux1 = (10**(-0.4*(magAB+48.6)))*(3.0*10**18/gratwave**2)
		if surf_bright:
			flux_input = 'mag_AB/arcsec^2'
		else:
			flux_input = 'magAB'
		print('OBJECT MAGNITUDE: ',magAB,' ',flux_input)
	#find wavelength and photons/Angstrom from flux
	print("flux",flux1)
	w,pA = wpA(flux1)

	#if flux set, establish units
	if flux != 'None':
		if surf_bright != 'None':
			if emline_w != 'None' and emline_w != None:
				flux_input = 'erg cm^-2 s^-1 arcsec^-2 in '+str(np.round(emline_w,1))+' Å'
			else: 
				flux_input = 'erg cm^-2 s^-1 arcsec^-2 Å^-1'
		else:
			if emline_w != 'None' and emline_w != None:
				flux_input = 'erg cm^-2 s^-1 in '+str(np.round(emline_w,1))+' Å'
			else:
				flux_input = 'erg cm^-2 s^-1 Å^-1'
		print('OBJECT FLUX: ',flux,' ',flux_input)
		if emline_w != 'None' and emline_w != None:
			print('EMISSION LINE OBJECT: flux is not per unit Å')
	
	#calculate Read Noise
	#need to determine pixels per spectral and spatial bin
	#1 unbinned pixel = 0.147 arcsec in both dimensions (no significant difference in grating angles
	#1 pixel spectral depends on grating: slicer BH L = 8 CCD pixels per spectral element = 4500/4500 = 1Å; M = 0.5Å; S = 0.25Å
	#1 pixel spectral depends on grating: slicer BM L = 8 CCD pixels per spectral element = 4500/2000 = 2.4Å; M = 1.2Å; S = 0.6Å
	#1 pixel spectral depends on grating: slicer BL L = 8 CCD pixels per spectral element = 4500/900 = 5Å; M = 2.5Å; S = 1.25Å
	print('SPAT',snr_spatial_bin,'SPEC',snr_spectral_bin)
	c_o = obj_cts(w,pA,grating,exptime)*snr_spatial_bin*snr_spectral_bin
	c_s = sky_cts(w,grating,exptime)*snr_spatial_bin*snr_spectral_bin

	np.set_printoptions(threshold=sys.maxsize)
	print(type(nframes),type(read_noise),type(pixels_per_snr_specbin),type(pixels_spat_bin),type(bin_factor))
	c_r = nframes*(read_noise**2)*pixels_per_snr_specbin*pixels_spat_bin*bin_factor
	snr = c_o/np.sqrt(c_s+c_o+c_r)

	#Signal to Noise Ratio
	p1 = figure(title='SNR',plot_width=350,plot_height=225,tools='hover')
	p1.line(w,snr)
	hover=p1.select(dict(type=HoverTool))
	hover.tooltips=[("Wavelength","@x"),("SNR","@y")]
	p1.add_tools(CrosshairTool())
	p1.xaxis.axis_label='Wavelength (Å)'
	p1.yaxis.axis_label='SNR / '+str(np.round(snr_spectral_bin,1))+' Å'
	script_snr,div_snr = components(p1)

	#Object Counts
	p2 = figure(title='SNR',plot_width=350,plot_height=225,tools='hover')
	p2.line(w,c_o)
	hover=p2.select(dict(type=HoverTool))
	hover.tooltips=[("Wavelength","@x"),("Object Counts","@y")]
	p2.add_tools(CrosshairTool())
	p2.xaxis.axis_label='Wavelength (Å)'
	p2.yaxis.axis_label='Object Counts / '+str(np.round(snr_spectral_bin,1))+' Å'
	script_obj,div_obj = components(p2)

	#Sky Counts
	p3 = figure(title='SNR',plot_width=350,plot_height=225,tools='hover')
	p3.line(w,c_s)
	hover=p3.select(dict(type=HoverTool))
	hover.tooltips=[("Wavelength","@x"),("Sky Counts","@y")]
	p3.add_tools(CrosshairTool())
	p3.xaxis.axis_label='Wavelength (Å)'
	p3.yaxis.axis_label='Sky Counts / '+str(np.round(snr_spectral_bin,1))+' Å'
	script_sky,div_sky = components(p3)
	
	#Read Noise Counts
	p4 = figure(title='SNR',plot_width=350,plot_height=225,tools='hover')
	p4.line(w,c_r*np.ones(len(w)))
	hover=p4.select(dict(type=HoverTool))
	hover.tooltips=[("Wavelength","@x"),("Read Noise Counts","@y")]
	p4.add_tools(CrosshairTool())
	p4.xaxis.axis_label='Wavelength (Å)'
	p4.yaxis.axis_label='Read Noise Counts / '+str(np.round(snr_spectral_bin,1))+' Å'
	script_rn,div_rn = components(p4)

	#Obj/Sky Counts
	p5 = figure(title='SNR',plot_width=350,plot_height=225,tools='hover')
	p5.line(w,c_o/c_s)
	p5.add_tools(CrosshairTool())
	hover=p5.select(dict(type=HoverTool))
	hover.tooltips=[("Wavelength","@x"),("Object/Sky Counts","@y")]
	p5.xaxis.axis_label='Wavelength (Å)'
	p5.yaxis.axis_label='Object Counts / Sky Counts'
	script_os,div_os = components(p5)
	
	#Flux
	p6 = figure(title='SNR',plot_width=350,plot_height=225,tools='hover')
	p6.line(w,pA)
	hover=p6.select(dict(type=HoverTool))
	hover.tooltips=[("Wavelength","@x"),("Flux","@y")]
	p6.add_tools(CrosshairTool())
	p6.xaxis.axis_label='Wavelength (Å)'
	p6.yaxis.axis_label='Flux [ph cm^-2 s^-1 Å^-1]'
	script_flux,div_flux = components(p6)

	return script_snr,div_snr,script_obj,div_obj,script_sky,div_sky,script_rn,div_rn,script_os,div_os,script_flux,div_flux

if __name__ == '__main__':
	do_calc('Small','BH2',4500.,0.75,3600.,'1x1',magAB=20)
	print('FINISHED')

