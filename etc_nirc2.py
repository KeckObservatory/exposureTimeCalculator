import numpy as np
import json

def do_calc(magnitude,          # source magnitude
            strehl,             # estimated strehl ratio
            exp_time,           # time per exposure, tint
            coadds,
            num_dith,           # number of dithers, 5, 9, etc.
            num_repeats,        # number of repeats in dithering pattern
            camera,             # camera size: wide or narrow
            img_filter,         # filter J H K Kp Lp Ms
            num_read,           # number of reads, 2 (CDS), 8, 16, 32, 64, 128
            x_extent,
            y_extent,
            ao_mode,            # AO: 0 NGS, 1 LGS
            laser_dith          # Laser: 0 fixed or 1 slaved
            ):


    samp_rate   = 200.      # unit in kherz, slow read out. normally 250 and can be as fast as 350 khz 
    overhead    = 6.001     # in percent

    ### lines below need to be clarified
    ### some variable definitions are unclear

    ## y_extent
    #todo: jriley: this introduces a bug below so i'm commenting it out for now
    #y_extent = y_extent + 8

    ## ticks??? ##
    ticks = 5e6 / (samp_rate * 1000.) # 5e6 seconds

    ## num_pause??? ##
    num_pause = ticks - 8.

    start_time = 240. * (1324. - y_extent) 

    ## rows_time ??? ##
    rows_time = 25. * y_extent * (num_pause * (8 + x_extent / 4) + (864 + 1.25 * x_extent)) 

    min_time_tmp = start_time + rows_time
    min_time_sampmode2 = (1 + overhead / 100) * min_time_tmp / 1e9 #6.001 percent overhead , convert from nanoseconds to seconds

    ## if samp_mode = 3, min_time_final = min_time* num_read 

    if num_read == 2:
        samp_mode = 2
        min_time = min_time_sampmode2
    elif num_read > 2:
        samp_mode = 3
        min_time = min_time_sampmode2*num_read


    ## if min itime < 2.5 miliseconds, set min itime to 0.0025 (sec)    
    if min_time < 2.5e-3:
        min_time = 2.5e-3

    read_min_time = min_time


    ### nirc2 efficiency ###
    ## function of: exp_time, coadds, num_repeats, num_read, samp_mode, x_extent, y_extent, ao_mode, num_dith, lgs_dith

    ## calculate overheads, which includes 
    ## integration time, coadds, number of images (num_dith * num_repeats)
    ## sampling mode, number of reads, window size (x_extent y_extent)
    ## ao_mode (NGS or LGS), lgs dither (fixed on target or dither slaved)

    ## overhead A associated with gathering and writing FITS header data
    overhead_write_fits_head = 3.05 #secconds
    ## overhead B from opening the file and writing data per pixel
    overhead_write_data = 4.65e-6 #seconds

    #todo: jriley: bug here due to "y_extent + 8"
    if x_extent == 1024 and y_extent == 1024 :
        read_time = 0.18 ## for CDS mode, x_extent, y_extent, samp_mode = 2, num_read = 2
    elif x_extent < 1024 and y_extent < 1024 :
        read_time = 0.05

    overhead_ao_dith = 6 # seconds
    overhead_lgs_dith = 15 # seconds

    if num_dith == 5 or num_dith == 9:
        num_ao_moves = num_dith
    elif num_dith == 1:
        num_ao_moves = 0
    else:
        num_ao_moves = num_dith + 1
        
       
    # 1 if laser slaved to dither, 0 if laser fixed on target
    if ao_mode == 0:
        lgs_dith = 0
    elif ao_mode == 1:
        lgs_dith = 1
        
    overhead_dith = overhead_ao_dith + overhead_lgs_dith * lgs_dith  

    time_per_coadd = exp_time + (num_read * read_min_time)

    overhead_write = overhead_write_fits_head + (overhead_write_data * x_extent * y_extent)

    overhead_move = overhead_dith * num_ao_moves

    ## number of exposure = coadds * number of dither * dither repeat

    num_exp = coadds*num_dith*num_repeats
    num_img = num_dith*num_repeats

    tot_exp_time = exp_time * num_exp
    tot_dith_time = num_img * (time_per_coadd * coadds + overhead_write) 
    tot_clock_time = overhead_move + tot_dith_time 

    tot_elps_obs_time = 6 * (num_dith+1) + (12*num_dith*num_repeats) + num_dith*num_repeats*coadds*(exp_time + read_time*(num_read-1))

    nirc2_eff = (tot_exp_time / tot_elps_obs_time) * 100
    overhead_nirc2 = tot_clock_time - overhead_move - tot_exp_time

    print(tot_exp_time)
    print(tot_elps_obs_time)
    print(tot_clock_time)


    ## minimum read time is a function of window (narrow, medium, wide), sampling mode, number of reads 

    ## effieicncy is a funciton of tint, coadds, number of images sampling mode number of reads, window, 
    ##  ao mode, number of dither, lgs dither


    ### nirc2 S/N calculation ###
    ## function of: magnitude, img_filter, camera, num_exp, exp_time, strehl, num_read

    noise_read = 56. ## 38 e-, drops as sqrt(reads) or 56 or 40 ???
    if num_read > 2:
        noise_read = noise_read / np.sqrt(num_read)

    num_pix = 50. ## 
    background = 634. ##background rate electrons per second, j = 634
    mag_zero = 25.1 
    gain = 4.0 ## electrons per DN



    #strehl = {'J':0.1, 'H':0.2, 'K':0.4, 'Kp':0.4, 'Lp':0.7, 'Ms':0.7}

    if camera == "wide":
        background = {'J':0.5, 'H':4.0, 'K':5.7, 'Kp':5.6, 'Lp':18535, 'Ms':18535 } 
        zero_point = {'J':26.9, 'H':26.96, 'K':26.18, 'Kp':26.30, 'Lp':25.08, 'Ms':22.87}
        num_pix = {'J':12.5, 'H':12.5, 'K':12.5, 'Kp':12.5, 'Lp':28.4, 'Ms':38.5}
    elif camera == "narrow":
        background = {'J':0.5, 'H':4.0, 'K':5.7, 'Kp':5.6, 'Lp':18535, 'Ms':18535 } 
        zero_point = {'J':26.9, 'H':26.96, 'K':26.18, 'Kp':26.30, 'Lp':25.08, 'Ms':22.87}
        num_pix = {'J':78.95, 'H':50.2, 'K':95.2, 'Kp':95.2, 'Lp':283.5, 'Ms':490.8}

    if camera == "wide":
        bg = background[img_filter] * gain * 16
    elif camera == "narrow":
        bg = background[img_filter] * gain

    mag_zero = zero_point[img_filter] + 2.5 * np.log10(strehl)

    signal = num_exp * exp_time * np.power(10,0.4*(mag_zero - magnitude))

    noise = np.sqrt(num_exp*np.square(noise_read)*num_pix[img_filter] + num_pix[img_filter]*background[img_filter]*num_exp*exp_time + signal)

    s2n = signal / noise 

    ap_area = num_pix[img_filter]

    noise_tot = noise / gain

    background_per_frame = bg * exp_time / gain
    if background_per_frame > 10000:
        print("Warning, background level in nonlinear regime.")
        
    signal_tot = signal / gain


    #debug
    print("Total noise is", noise_tot, "DN")
    print("Total signal is", signal_tot, "DN, in", ap_area, "pixels")
    print("S/N =", s2n)
    print("Aperature area is", ap_area, "pixels")
    print("Background per frame is", background_per_frame, "DN")
    print("Efficiency = ", nirc2_eff, "%")
    print("Total integration time is", tot_exp_time, "sec")
    print("Total elapsed observing time is", tot_elps_obs_time, "sec")
    print("Total clock time is", tot_clock_time, "sec" )


    #return final data dictionary result
    result = {}
    result['noise_tot'] = noise_tot
    result['signal_tot'] = signal_tot
    result['ap_area'] = ap_area
    result['s2n'] = s2n
    result['ap_area'] = ap_area
    result['background_per_frame'] = background_per_frame
    result['nirc2_eff'] = nirc2_eff
    result['tot_exp_time'] = tot_exp_time
    result['tot_elps_obs_time'] = tot_elps_obs_time
    result['tot_clock_time'] = tot_clock_time
    return json.dumps(result)


#This code below will only run if this module is run from the command line
#example: python etc_nirc2.py
if __name__ == "__main__":

    #test
    data = do_calc( magnitude   = 15, 
                    strehl      = 0.3,  
                    exp_time    = 10,  
                    coadds      = 10, 
                    num_dith    = 5,  
                    num_repeats = 1,  
                    x_extent    = 512, 
                    y_extent    = 512, 
                    camera      = "narrow", 
                    img_filter  = "J", 
                    num_read    = 16, 
                    ao_mode     = 0,  
                    laser_dith  = 0
                    )
    print ('result = ', data)
