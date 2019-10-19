def instr_calc():
	instr='DEIMOS'
	grating = request.form.get("grating",'600ZD')
    wavec = request.form.get("wavec",'7000')
    slitwidth = request.form.get("slitwidth",'0.75 arcsec')
    binning = request.form.get("binning",'1x1 pixels')
    exptime = request.form.get("exptime",1200.0)
    seeing = request.form.get("seeing",0.75)
    airmass = request.form.get("airmass",1.1)
    magnitude = request.form.get("magnitude",22.0)
    ffilter = request.form.get("ffilter",'r')
    mtype = request.form.get("mtype",'AB')
    template = request.form.get("template",'flat')
    redshift = request.form.get("redshift",0.0)
    
    calculate = request.form.get("calculate")
    if calculate:
        ###### FILL WITH HIRES DOCALC ######
        pass
        ###### FILL WITH HIRES DOCALC ######


        # fltmag = float(magnitude)
        # fltstrehl = float(strehl)
        # inttpe = int(tpe)
        # intcoadds = int(coadds)
        # intndithers = int(ndithers)
        # intrpd = int(rpd)
        # intwin = int(window)
        # intnr = int(nreads)
        # result = nirc2_docalc(magnitude   = fltmag, 
        #                       strehl      = fltstrehl,  
        #                       exp_time    = inttpe,  
        #                       coadds      = intcoadds, 
        #                       num_dith    = intndithers,  
        #                       num_repeats = intrpd,  
        #                       x_extent    = intwin, 
        #                       y_extent    = intwin, 
        #                       camera      = cam, 
        #                       img_filter  = fltr, 
        #                       num_read    = intnr, 
        #                       ao_mode     = ao_bin,  
        #                       laser_dith  = laser_bin)

        # result = {k:round(v,2) for k, v in result.items()}

        # snr = result.get('s2n')
        # tot_signal = result.get('signal_tot')
        # aperature_area = result.get('ap_area')
        # tot_noise = result.get('noise_tot')
        # bkg_per_frame = result.get('background_per_frame')
        # efficiency = result.get('nirc2_eff')
        # # ao_tel_ovrhd = 
        # # nirc2_ovrhd = result.get('nirc2_eff') 
        # tot_int = result.get('tot_exp_time')
        # tot_elaptime = result.get('tot_elps_obs_time')
        # return render_template('etc.html',instrument=instr,magnitude=magnitude,tpe=tpe,ndithers=ndithers,
        #                     strehl=strehl,coadds=intcoadds,rpd=rpd,camera=camera,filter=fltr,
        #                     nreads=intnr,window=intwin,ao=ao,laser=laser,snr=snr,tot_signal=tot_signal,
        #                     aperature_area=aperature_area,tot_noise=tot_noise,bkg_per_frame=bkg_per_frame,
        #                     efficiency=efficiency,tot_int=tot_int,tot_elaptime=tot_elaptime,calculate=calculate)

    return render_template('etc_deimos.html',instrument=instr,grating=grating,wavec=wavec,
                            slitwidth=slitwidth,binning=binning,exptime=exptime,
                            seeing=seeing,airmass=airmass,magnitude=magnitude,ffilter=ffilter,mtype=mtype,
                            template=template,redshift=redshift,calculate=calculate)