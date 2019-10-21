import numpy as np
import datetime
from flask import Flask, render_template, request, redirect, url_for, session
from etc_nires import do_calc as nires_docalc
import json

def instr_calc(allowed):
    instr = 'NIRES'
    exptime = request.form.get("exptime",3600.0)
    seeing = request.form.get("seeing",0.8)
    teff = request.form.get("teff",5000)
    coadds = request.form.get("coadds",1)
    obswave = request.form.get("obswave",2.2)
    dither = request.form.get("dither",'AB')
    rpd = request.form.get("rpd",1)
    nreads = request.form.get("nreads",16)
    magnitude = request.form.get("magnitude",18.2)
    src_type = request.form.get("src_type","Point Source")
    redshift = request.form.get("redshift",2)
    print(src_type)
    
    calculate = request.form.get("calculate")
    if calculate:
        if "Point" in src_type:
          src_tcalc = "PointSource"
        elif src_type == "QSO":
          src_tcalc = "qso"

        redshift = int(redshift)
        seeing = float(seeing)
        obswave = float(obswave)
        magnitude = float(magnitude)
        exptime = float(exptime)
        nreads = int(nreads)
        rpd = int(rpd)
        coadds = int(coadds)
        teff = int(teff)
        ###### FILL WITH ESI DOCALC ######
        result = nires_docalc(src_type      = src_tcalc,
                              mag_src       = magnitude, 
                              obs_wave      = obswave, 
                              exp_time      = exptime, 
                              coadds        = coadds, 
                              dither        = dither, 
                              dither_repeat = rpd,
                              seeing        = seeing, 
                              num_reads     = nreads,
                              redshift      = redshift)

        s2n = np.round(result,2)
        print("S2N:",s2n)

        return render_template('etc_nires.html',instrument=instr,exptime=exptime,seeing=seeing,teff=teff,
                               obswave=obswave,dither=dither,rpd=rpd,nreads=nreads,magnitude=magnitude,
                               src_type=src_type,redshift=redshift,coadds=coadds,calculate=calculate,s2n=s2n,
                               allowed=allowed)

    return render_template('etc_nires.html',instrument=instr,exptime=exptime,seeing=seeing,teff=teff,
                            obswave=obswave,dither=dither,rpd=rpd,nreads=nreads,magnitude=magnitude,
                            src_type=src_type,redshift=redshift,coadds=coadds,calculate=calculate,allowed=allowed)