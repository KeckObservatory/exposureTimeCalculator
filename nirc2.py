import numpy as np
import datetime
from flask import Flask, render_template, request, redirect, url_for, session
from etc_nirc2 import do_calc as nirc2_docalc
import json

def instr_calc(allowed):
	instr = 'NIRC2'
	magnitude = request.form.get("magnitude",20)
	tpe = request.form.get("tpe",10)
	ndithers = request.form.get("ndithers",5)
	strehl = request.form.get("strehl",0.3)
	coadds = request.form.get("coadds",10)
	rpd = request.form.get("rpd",1)
	camera = request.form.get("camera",'Narrow')
	fltr = request.form.get("filter",'Kp')
	nreads = request.form.get("nreads",2)
	window = request.form.get("window",1024)
	ao = request.form.get("ao",'NGS')
	laser = request.form.get("laser",'Dither to Object')
	calculate = request.form.get("calculate")
	if calculate:
	    if camera == 'Narrow':
	        cam = 'narrow'
	    elif camera == 'Wide':
	        cam = 'wide'

	    if ao == 'NGS':
	        ao_bin = 0
	    elif ao == 'LGS':
	        ao_bin = 1

	    if laser == 'Fixed to Center':
	        laser_bin = 0
	    elif laser == 'Dither to Object':
	        laser_bin = 1

	    fltmag = float(magnitude)
	    fltstrehl = float(strehl)
	    inttpe = int(tpe)
	    intcoadds = int(coadds)
	    intndithers = int(ndithers)
	    intrpd = int(rpd)
	    intwin = int(window)
	    intnr = int(nreads)
	    result = nirc2_docalc(magnitude   = fltmag, 
	                          strehl      = fltstrehl,  
	                          exp_time    = inttpe,  
	                          coadds      = intcoadds, 
	                          num_dith    = intndithers,  
	                          num_repeats = intrpd,  
	                          x_extent    = intwin, 
	                          y_extent    = intwin, 
	                          camera      = cam, 
	                          img_filter  = fltr, 
	                          num_read    = intnr, 
	                          ao_mode     = ao_bin,  
	                          laser_dith  = laser_bin)

	    result = json.loads(result)
	    result = {k:round(v,2) for k, v in result.items()}
	    # print(result)
	    snr = result.get('s2n')
	    tot_signal = result.get('signal_tot')
	    aperature_area = result.get('ap_area')
	    tot_noise = result.get('noise_tot')
	    bkg_per_frame = result.get('background_per_frame')
	    efficiency = result.get('nirc2_eff')
	    # ao_tel_ovrhd = 
	    # nirc2_ovrhd = result.get('nirc2_eff') 
	    tot_int = result.get('tot_exp_time')
	    tot_elaptime = result.get('tot_elps_obs_time')
	    return render_template('etc_nirc2.html',instrument=instr,magnitude=magnitude,tpe=tpe,ndithers=ndithers,
	                        strehl=strehl,coadds=intcoadds,rpd=rpd,camera=camera,filter=fltr,
	                        nreads=intnr,window=intwin,ao=ao,laser=laser,snr=snr,tot_signal=tot_signal,
	                        aperature_area=aperature_area,tot_noise=tot_noise,bkg_per_frame=bkg_per_frame,
	                        efficiency=efficiency,tot_int=tot_int,tot_elaptime=tot_elaptime,calculate=calculate,
	                        allowed=allowed)
	print(instr)
	return render_template('etc_nirc2.html',instrument=instr,magnitude=magnitude,tpe=tpe,ndithers=ndithers,
	                        strehl=strehl,coadds=coadds,rpd=rpd,camera=camera,filter=fltr,
	                        nreads=nreads,window=window,ao=ao,laser=laser,calculate=calculate,allowed=allowed)