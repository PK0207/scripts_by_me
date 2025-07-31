;; Comparing COS H2 emission lines in UV to Rotating H2 disk Models
;;
;;
;; Keri Hoadley
;; Outline - 16 December 2013
;; Assumptions + Approximations - 11 April 2013
;; Parameters and Limits of Parameters -
;;
;
; Purpose
; -------
; With models created, we have to compare via chi-squared minimization and Eliot's 
; v_rotation code to determine the velocity shift in the model code, to match the centers
; of the emission lines (model and data). This code will store chi-sqr values for a given
; parameter set of the model-to-data comparison.
; The chi-sqr values will be stored in an N-dimensional array, representative of the
; parameter space of the models. The location of a chi-sqr value in the scheme of of the N-dim
; array will give info on which parameters it is associated with.
;
;
; Outline
; -------
; 1. Read in the star directory and search through the directory for info on the rest emission wavelength
;    of the H2.
;    
;    b. Using the velocity shift code Eliot wrote: For each wavelength, open the data file for the target
;       (actually, the velocity-space of the same wavelength) and determine the velocity shift from v=0.
;       
;    c. Store the velocity shift value to later shift all model data by that value.
;    
;    d. In the same directory, there will be a SAV file with all the parameter information used in the models.
;       Using this information, make an N-dim chi-sqr array to store all calculated chi-squared values
;       into this array that can give us info on which parameters fit well with our data.
;       
; 2. Within a given rest wavelength directory, again search through the directory for a complete list of
;    all the models for that wavelength. Using the lengths of the parameter space of the models (see 1.d.),
;    write a nested FOR loop to search through all the parameter space.
;    
;    b. Convolve the model data to the COS-LSF function.
;    
;    c. Calculate the chi-squared value for each set of parameters, storing the chi-square value in the
;       appropriate spot within the chi-sqr array.
; 
; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Revisions - October 2014
; ---------
; To accomodate the new file structure of the SAVE file read in (produced in disk_profile_v2),
; we need a revised chi2 calculation program to read in each model, pair the calculated models
; with the observed H2 emission wavelengths, and compare in the chi2 calculation.
; 
; 
; 
;






PRO H2models_v_COSdata_chisqr_v2, arg=arg
;device,decompose=0
;setplotcolors


;;;;;;;;;;;;;;;;;;;;;;; FOR pc ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; List all reference paths and files relevant to restore.
; Main reference path to follow:
;ref_path = 'C:\Users\keho8439\Documents\Comps2\'
;ref_path = 'F:\H2disks_models\'
;
;; Observation file paths:
;ref_dir_obs = 'C:\Users\keho8439\Documents\Comps2\obs\'+['DFTau\']
;file_obs_vel = ['V-DF-TAU_H2velocity.sav']
;file_obs_flux = ['V-DF-TAU_H2flux.sav']
;file_obs_fluxerr = ['V-DF-TAU_H2fluxerr.sav']
;file_obs_wave = ['V-DF-TAU_H2wavelengths.sav']
;parm_list = 'parameter_list.sav'
;
;; Model file paths:
;;ref_dir_models = 'models\'+['BPTau\']
;ref_dir_models = ['DFTau\']
;ref_dir_wavelengths = ['1467.08']
;file_models = '\'+['V-AA-TAU000000.sav']

;ref_path = 'C:\Users\keho8439\Documents\Comps2\'
;targ = 'LkCa15';'AATau'
;
;; Observation file paths:
;ref_dir_obs = 'obs\'+['LkCa15\']
;file_obs_vel = ['EM-LKCA-15_H2velocity.sav'];['V-DF-TAU_H2velocity.sav']
;file_obs_flux = ['EM-LKCA-15_H2flux.sav'];['V-DF-TAU_H2flux.sav']
;file_obs_fluxerr = ['EM-LKCA-15_H2fluxerr.sav'] ;['V-DF-TAU_H2fluxerr.sav']
;file_obs_wave = ['EM-LKCA-15_H2wavelengths.sav'];['V-DF-TAU_H2wavelengths.sav']
;
;; Model file paths:
;ref_dir_models = 'F:\H2disks_models\'+['LkCa15\'];'AATau\']
;ref_dir_model_res = 'models\'+['LkCa15\']
;parm_list = 'parameter_list.sav'

;
;ref_path = 'C:\Users\keho8439\Documents\Comps2\'
;targ = 'AATau'
;
;; Observation file paths:
;ref_dir_obs = 'obs\'+['AATau\']
;file_obs_vel = ['V-AA-TAU_H2velocity.sav']
;file_obs_flux = ['V-AA-TAU_H2flux.sav']
;file_obs_fluxerr = ['V-AA-TAU_H2fluxerr.sav']
;file_obs_wave = ['V-AA-TAU_H2wavelengths.sav']
;
;; Model file paths:
;ref_dir_models = 'F:\H2disks_models\'+['AATau\']
;ref_dir_model_res = 'models\'+['AATau\']+'disk_profiles_v2\'
;parm_list = 'parameter_list.sav'

; List of progressions, by rest wavelength (in 2d array)
progression = [ [1338.56, 1398.95, 1460.17, 1521.59, 1338.57],$        ; [0,1] excited state
                [1342.26, 1393.96, 1402.65, 1454.97, 1463.83],$        ; [0,2]
                [1446.12, 1489.57, 1504.76, 0, 0],$              ; [1,4]
                [1442.87, 1467.08, 1500.45, 1524.65, 1580.67] ]        ; [1,7]
                
prog_name = [ '[0,1]','[0,2]','[1,4]','[1,7]' ]


;;;;;;;;;;;;;;;;;;;;; FOR arlsrv2 ;;;;;;;;;;;;;;;;;;;;;;;;;
;; List all reference paths and files relevant to restore.
;; Main reference path to follow:
;targ = ['AATau','CSCha','DFTau','DMTau','BPTau','GMAur']
;ref_path_home = '/arlstorage/home/student/keho8439/Research/Comps2/'
;ref_path_data = '/arlstorage/data/student/keho8439/'+targ+'/disk_profiles_v2/'
;
;;
                 
;targ = ['AATau','CSCha','DFTau','DMTau','BPTau','GMAur','LKCA15','UXTau','V4046SGR']
;
;;;; Observation file paths:
;ref_dir_obs = 'obs/'+targ+'/'
;file_obs_vel = ['V-AA-TAU_H2velocity.sav','V-CS-CHA_H2velocity.sav','V-DF-TAU_H2velocity.sav',$
;                'V-DM-TAU_H2velocity.sav','V-BP-TAU_H2velocity.sav','V-GM-AUR_H2velocity.sav',$
;                'EM-LKCA-15_H2velocity.sav','NAME-UX-TAU-A_H2velocity.sav','V4046SGR_H2velocity.sav']
;                
;file_obs_flux = ['V-AA-TAU_H2flux.sav','V-CS-CHA_H2flux.sav','V-DF-TAU_H2flux.sav',$
;                 'V-DM-TAU_H2flux.sav','V-BP-TAU_H2flux.sav','V-GM-AUR_H2flux.sav',$
;                 'EM-LKCA-15_H2flux.sav','NAME-UX-TAU-A_H2flux.sav','V4046SGR_H2flux.sav']
;                 
;file_obs_fluxerr = ['V-AA-TAU_H2fluxerr.sav','V-CS-CHA_H2fluxerr.sav','V-DF-TAU_H2fluxerr.sav',$
;                    'V-DM-TAU_H2fluxerr.sav','V-BP-TAU_H2fluxerr.sav','V-GM-AUR_H2fluxerr.sav',$
;                    'EM-LKCA-15_H2fluxerr.sav','NAME-UX-TAU-A_H2fluxerr.sav','V4046SGR_H2fluxerr.sav']
;                    
;file_obs_wave = ['V-AA-TAU_H2wavelengths.sav','V-CS-CHA_H2wavelengths.sav','V-DF-TAU_H2wavelengths.sav',$
;                 'V-DM-TAU_H2wavelengths.sav','V-BP-TAU_H2wavelengths.sav','V-GM-AUR_H2wavelengths.sav',$
;                 'EM-LKCA-15_H2wavelengths.sav','NAME-UX-TAU-A_H2wavelengths.sav','V4046SGR_H2wavelengths.sav']
;
;parm_list = 'parameter_list.sav'
;ref_path_home = '/arlstorage/home/student/keho8439/Research/Comps2/'
;ref_path_data = '/arlstorage/data/student/keho8439/'+targ+'/disk_profiles_v2/'



targ = ['AATau','HNTau','RECX11','RECX15','SUAur']    ; No TW Hya yet

;;; Observation file paths:
ref_dir_obs = 'obs/'+targ+'/'
file_obs_vel = ['V-AA-TAU-2013_H2velocity.sav','V-HN-TAU_H2velocity.sav','RECX-11_H2velocity.sav',$
  'RECX-15_H2velocity.sav','V-SU-AUR_H2velocity.sav']
  
file_obs_flux = ['V-AA-TAU-2013_H2flux.sav','V-HN-TAU_H2flux.sav','RECX-11_H2flux.sav',$
  'RECX-15_H2flux.sav','V-SU-AUR_H2flux.sav']
  
file_obs_fluxerr = ['V-AA-TAU-2013_H2fluxerr.sav','V-HN-TAU_H2fluxerr.sav','RECX-11_H2fluxerr.sav',$
  'RECX-15_H2fluxerr.sav','V-SU-AUR_H2fluxerr.sav']
  
file_obs_wave = ['V-AA-TAU-2013_H2wavelengths.sav','V-HN-TAU_H2wavelengths.sav','RECX-11_H2wavelengths.sav',$
  'RECX-15_H2wavelengths.sav','V-SU-AUR_H2wavelengths.sav']
  
parm_list = 'parameter_list.sav'
ref_path_home = '/arlstorage/home/student/keho8439/Research/Comps2/'
ref_path_data = '/arlstorage/data/student/keho8439/'+targ+'/disk_profiles_v2/'



;
;; Observation file paths:
;ref_dir_obs = 'obs/'+targ+'/'
;file_obs_vel = ['V-DF-TAU_H2velocity.sav']
;file_obs_flux = ['V-DF-TAU_H2flux.sav']
;file_obs_fluxerr = ['V-DF-TAU_H2fluxerr.sav']
;file_obs_wave = ['V-DF-TAU_H2wavelengths.sav']
parm_list = 'parameter_list.sav'
;;
;; Model file paths:
ref_dir_models = 'models/'+targ+'/'
;


; FOR loop: Open the model directory, make a list of all the wavelengths calculated for the models,
; and then open the velocity curve SAV file of the wavelengths to make a velocity offset array for
; the given emission wavelength.

If keyword_set(arg) THEN BEGIN
  start_targ = arg
  end_targ = arg
ENDIF ELSE BEGIN
  start_targ = 0
  end_targ = N_ELEMENTS(targ) - 1
ENDELSE

FOR target = start_targ, end_targ DO BEGIN
  print,'Starting chi^2 statistics for full model-data comparison for '+targ[target]
  
  ; Open an ASCII file to save the min(chi-sqr) parameters for wach wavelength
  openw,lun3,ref_path_home+ref_dir_models[target]+'min_chisqr_redwing.dat',/get_lun,/APPEND
  openw,lun2,ref_path_home+ref_dir_models[target]+'min_chisqr_per_emission.dat',/get_lun,/APPEND
  openw,lun1,ref_path_home+ref_dir_models[target]+'min_chisqr_per_progression.dat',/get_lun,/APPEND
  ; Print the opening line of this file
  opening_string_1 = 'Wave (A)'+STRING(9B)+'z (Hp)'+STRING(9B)+'gamma'+STRING(9B)+'T1AU (K)'+STRING(9B)+$
                     'q'+STRING(9B)+'r_char (AU)'+STRING(9B)+'M_H2 (Msol)'+STRING(9B)+'Min. chi sqr'
  opening_string_1a = 'Progression'+STRING(9B)+'z (Hp)'+STRING(9B)+'gamma'+STRING(9B)+'T1AU (K)'+STRING(9B)+$
                     'q'+STRING(9B)+'r_char (AU)'+STRING(9B)+'M_H2 (Msol)'+STRING(9B)+'Min. chi sqr'
  opening_string_2 = '--------'+STRING(9B)+'------'+STRING(9B)+'-----'+STRING(9B)+'--------'+STRING(9B)+$
                     '-'+STRING(9B)+'-----------'+STRING(9B)+'-----------'+STRING(9B)+'------------'
  opening_string_2a = '------------'+STRING(9B)+'------'+STRING(9B)+'-----'+STRING(9B)+'--------'+STRING(9B)+$
                     '-'+STRING(9B)+'-----------'+STRING(9B)+'-----------'+STRING(9B)+'------------'
  printf,lun2,opening_string_1
  printf,lun2,opening_string_2
  
  printf,lun3,opening_string_1
  printf,lun3,opening_string_2
  
  printf,lun1,opening_string_1a
  printf,lun1,opening_string_2a
  

;  stop
  
  ; Put all files into an array to loop through
  list_model_files = FILE_SEARCH(ref_path_data[target]+'*-*')
  
  
  
  ; Using Eliot's velocity shift code:
  ; Open up the observation file with the velocity cruves.
  ; Locate the appropriate wavelength, its index, and the resulting velocity/flux arrays
  ; Using the velocity shift function, store the offsets we need to add to each model
  ; at the appropriate emission line, such that the emission_wave index of a given rest wavelength matches the
  ; velocity offset we need to use.
  restore, ref_path_home+ref_dir_obs[target]+file_obs_vel[target]         ; Array: velocity
  restore, ref_path_home+ref_dir_obs[target]+file_obs_flux[target]        ; Array: subflux
  restore, ref_path_home+ref_dir_obs[target]+file_obs_fluxerr[target]     ; Array: suberr
  restore, ref_path_home+ref_dir_obs[target]+file_obs_wave[target]        ; Array: labwave
  
  
  ; Open the parameter list file, noting the names of each parameter.
  ; In the next loop (and subsequant loops), we will use the number of elements in each array
  ; to determine where which chi-squared values are stored, to represent the chi-squared
  ; fit to the data for a given set of parameters.
  restore, '/arlstorage/data/student/keho8439/'+targ[target]+'/'+parm_list    ; arrays (in order): T_1AU, q, gamma, 
                                              ;                    r_char, M_H2, Hp_height
  ; q values at index 4 is messed up (is "0.05" instead of "-0.05", so define that here)
;  q[4] = -0.05
  
  ; With the velocity shift stored, create a large chisqrmin array to store all chi squared values
  ; associated with a given parameter set.
  ; This will change with changing labwave observed, so that is why I'm make it now and then
  ; saving the chi-sqr arrays by wavelength in the /stats folder.
  chisqr_models_ind = DBLARR( N_ELEMENTS(labwave), N_ELEMENTS(M_H2), N_ELEMENTS(r_char), $
    N_ELEMENTS(q), N_ELEMENTS(T_1AU), N_ELEMENTS(gamma), N_ELEMENTS(Hp_height) )      ; chisqr by emission line
  chisqr_models_prog = DBLARR( N_ELEMENTS(progression[0,*]), N_ELEMENTS(M_H2), N_ELEMENTS(r_char), $
    N_ELEMENTS(q), N_ELEMENTS(T_1AU), N_ELEMENTS(gamma), N_ELEMENTS(Hp_height) )      ; chisqr by progression
  chisqr_models_red = DBLARR( N_ELEMENTS(labwave), N_ELEMENTS(M_H2), N_ELEMENTS(r_char), $
    N_ELEMENTS(q), N_ELEMENTS(T_1AU), N_ELEMENTS(gamma), N_ELEMENTS(Hp_height) )      ; chisqr by emission, red wing only

  
  ; Calculate the velocity shift from center of all the observed emission lines
  ; from fit_velocity_gauss: [height, center, sigma, fwhm]
  start = 50
  finish = 350
  parms_vshift = DBLARR( 4, N_ELEMENTS(labwave) )     ; 4 = the number of elements returned from fit_vel_gauss
  FOR i = 0, N_ELEMENTS(labwave) - 1 DO parms_vshift[*,i] = fit_velocity_gauss(velocity[start:finish,i], subflux[start:finish,i], suberr[start:finish,i])

;  stop
    
  ; Store the v_shift value before moving on to compare the observations to the models
  v_offset = parms_vshift[1,*]    ; Gives the center (centroid) of the observed emission line.ar
    
  ;;;;;;;;;;;;;;;;;;;;;; SUBTRACT OFF CONTINUUM ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;
  ; Each wavelength in the observations is offset from 0 flux by some continuum value.
  ; First, look at [-200,-100] - [100,200] km/s and fit a linear fit to said data.
;  v_linear = [velocity[0:100,wave_index],velocity[300:400,wave_index]]
;  f_linear = [subflux[0:100,wave_index],subflux[300:400,wave_index]]
  v_linear = [velocity[30:150,*],velocity[260:390,*]]
  f_linear = [subflux[30:150,*],subflux[260:390,*]]
  e_linear = [suberr[30:150,*],suberr[260:390,*]]
  c_fit = DBLARR( 2,N_ELEMENTS(labwave) )
  FOR i = 0, N_ELEMENTS(labwave) - 1 DO c_fit[*,i] = poly_fit(v_linear[*,i],f_linear[*,i],1);,MEASURE_ERRORS=e_linear)
  linefit = DBLARR( N_ELEMENTS(velocity[*,0]),N_ELEMENTS(labwave) )
  FOR i = 0, N_ELEMENTS(labwave) - 1 DO linefit[*,i] = poly(velocity[*,i],c_fit[*,i])
;  stop
    
  ; Subtract the linefit array from the data array
  subflux = subflux - linefit
;  stop
    
    ; FOR loop: Loop through all files in list_model_files and determine where which chi-squared results
    ; will be places based on matching up the chi-squared index with each parameter index.
    
    FOR i = 0, N_ELEMENTS(list_model_files) - 1 DO BEGIN
      print,'Iterating i = '+strtrim(i,2)+' of '+strtrim(N_ELEMENTS(list_model_files)-1,2)
      ; Restore the model file
      restore,list_model_files[i]   ; Arrays: starname, lambda_ref, temperature, gradient, 
                                    ;         density_gamma, disk_height, characteristic_radius, 
                                    ;         totalH2mass, v_obs_bin, intensity_binned
      
      ; Determine the indexing of the models to the observations (ONLY for first model in list:
      ; this will be the same for all other models)

      IF i EQ 0 THEN BEGIN
        emission_waves = INTARR( N_ELEMENTS(labwave) )
        ; Do chi-sqr minimization for each model/flux combo. First, add the the offset velocity
        ; shift to the model velocity array.
        v_model = DBLARR( N_ELEMENTS(v_obs_bin),N_ELEMENTS(labwave) )
        
        FOR j = 0, N_ELEMENTS(labwave) - 1 DO emission_waves[j] = WHERE(labwave[j] EQ lambda_ref)
        for j = 0, N_ELEMENTS(labwave) - 1 DO v_model[*,j] = v_obs_bin+v_offset[j]
      ENDIF
      
      ; Make an intensity array that matches the same array placements of subflux to intensity_binned
      intensity = DBLARR( N_ELEMENTS(v_obs_bin), N_ELEMENTS(labwave) )
      for j = 0, N_ELEMENTS(labwave) - 1 DO intensity[*,j] = intensity_binned[*,emission_waves[j]]
      
;      stop
      
      ; Determine each parameter index value
      ; Will have ~ 8 parameters to do this for.
      nparms = 6    ; m_h2, r_c, q, t_1au, gamma, hp
      i_t1au = WHERE(temperature EQ T_1AU)
      i_q = WHERE(gradient EQ q)
      i_gamma = WHERE(density_gamma EQ gamma)
      i_hp = WHERE(disk_height EQ Hp_height)
      i_rch = WHERE(characteristic_radius EQ r_char)
      i_mh2 = WHERE(totalH2mass EQ M_H2)
        
;      stop
;      
     
      intensity[N_ELEMENTS(intensity[*,0])-1,*] = 0.0     ; Remove last array element with random errors
      
;      flux_model = interpol(intensity,v_model,velocity)
      

      ; Convole the COS LSF with the interpolated model.
      ; Convolve the model v_obs_bin and intensity data with Eliot's
      ; convolve_data function
      flux_model_conv = DBLARR(N_ELEMENTS(v_model[*,0]),N_ELEMENTS(labwave) )
      for j = 0, N_ELEMENTS(labwave) - 1 DO flux_model_conv[*,j] = convolve_data(intensity[*,j],labwave[j])
      
 ;     stop
      ; Also, need to interpolate the flux_model values, based on the COS velocity values.
      ; Make a new model array with the same number of elements as the COS subflux array, and
      ; using the COS velocity solution, interpolate the model flux values, based on the
      ; velocity match-up we have, to find the model flux at a given COS velocity position.
   ;   flux_model = INTERPOL(flux_model_conv, v_model, velocity);, /QUADRATIC)
      flux_model = DBLARR(N_ELEMENTS(velocity[*,0]),N_ELEMENTS(labwave) )
      for j = 0, N_ELEMENTS(labwave) - 1 DO flux_model[*,j] = INTERPOL(flux_model_conv[*,j], v_model[*,j], velocity[*,j])
  ;    stop
      

      ; Make an error bar checker to avoid including error bars that are stupid low 
      ; out in the continue.
      err_mx = INTARR( N_ELEMENTS(labwave) )
      avg_err = DBLARR( N_ELEMENTS(labwave) )
      FOR j = 0, N_ELEMENTS(labwave) - 1 DO BEGIN 
        err_mx[j] = WHERE(subflux[*,j] EQ max(subflux[180:240,j]))
        avg_err[j] = avg(suberr[err_mx[j]-10:err_mx[j]+10,j])
      ENDFOR
  
      
      ; Make a chisqr_sum counter to calculate the chi squared min for each parameter.
      chisqr_sum = DBLARR( N_ELEMENTS(labwave) )
      cs_sum_red = DBLARR( N_ELEMENTS(labwave) )
      cs_prog = DBLARR( N_ELEMENTS(progression[0,*]) )
      
      
      start = 160       ; index in subflux to start at for emission feature comparison to models
      start_red = INTARR( N_ELEMENTS(labwave) )
      FOR j = 0, N_ELEMENTS(labwave) - 1 DO BEGIN
        temp =  WHERE(velocity[*,j] GE v_offset[j]);+15.0)    ; red wind of emission line only
        start_red[j] = temp[0]      ; first instance of velocity meeting criteria
      ENDFOR 
      finish = 250      ; ending index
      
      ; Loop through the model,cos data arrays and calculate the chi-squared
      ; INCLUDES ONLY v ~ [-200,200]
      FOR k = start, finish DO BEGIN
        ;chisqr_sum += (subflux[k,wave_index] - flux_model_conv[k])^2 / (suberr[k,wave_index])^2
        FOR j = 0, N_ELEMENTS(labwave) - 1 DO BEGIN
          p_ind = array_indices(progression,WHERE(labwave[j] EQ progression))     ; use p_ind[1] to store into cs_prog
;          stop
          IF (subflux[k,j]/suberr[k,j]) GE (subflux[err_mx[j],j]/avg_err[j])*2. THEN BEGIN
            cs_prog[p_ind[1]] += (subflux[k,j] - flux_model[k,j])^2 / (avg_err[j]) ;(flux_model[k,j])
            chisqr_sum[j] += (subflux[k,j] - flux_model[k,j])^2 / (avg_err[j]) ;(flux_model[k,j])
            IF k GE start_red[j] THEN cs_sum_red[j] += (subflux[k,j] - flux_model[k,j])^2 / (avg_err[j])
          ENDIF ELSE BEGIN
            cs_prog[p_ind[1]] += (subflux[k,j] - flux_model[k,j])^2 / (suberr[k,j])^2
            chisqr_sum[j] += (subflux[k,j] - flux_model[k,j])^2 / (suberr[k,j])^2
            IF k GE start_red[j] THEN cs_sum_red[j] += (subflux[k,j] - flux_model[k,j])^2 / (suberr[k,j])^2
          ENDELSE
        ENDFOR
      ENDFOR
;      stop
      ; Store the sum in the appropriate array position
      chisqr_models_ind[*,i_mh2,i_rch,i_q,i_t1au,i_gamma,i_hp] = chisqr_sum / (N_ELEMENTS(subflux[start:finish,0])-nparms-1)
      chisqr_models_prog[*,i_mh2,i_rch,i_q,i_t1au,i_gamma,i_hp] = cs_prog / (3*(N_ELEMENTS(subflux[start:finish,0])-nparms-1))
      FOR j = 0, N_ELEMENTS(labwave) - 1 DO $
        chisqr_models_red[j,i_mh2,i_rch,i_q,i_t1au,i_gamma,i_hp] = cs_sum_red[j] / (N_ELEMENTS(subflux[start_red[j]:finish,0])-nparms-1)
;      stop
      
    ENDFOR
    
    ; Save the chisqr models results in a SAVE file, to play around with later.
    FILE_MKDIR,ref_path_home+ref_dir_models[target]+'chi_square_v2/'
    file = ref_path_home+ref_dir_models[target]+'chi_square_v2/'+'chisqr_array_1.sav'
    SAVE,chisqr_models_ind,chisqr_models_prog,chisqr_models_red,filename=file
    
    ; In the open ascii file, save the parameters that give us the minimum chi-square value for a given wavelength
    ; Save the wavelength: parameters: min(chi-sqr)
;    restore,ref_path+ref_dir_model_res[target]+parm_list
    
    ; Find where chi-sqr is minimized, and determine the array indices of the chi-sqr min.
    ; separate each labwave row in text file
    FOR k = 0, N_ELEMENTS(labwave) - 1 DO BEGIN
      mn_chi = min(chisqr_models_red[k,*,*,*,*,*,*],location)
      ind = array_indices(chisqr_models_red[k,*,*,*,*,*,*],location)
      string_3 = STRTRIM(lambda_ref[emission_waves[k]],2)+STRING(9B)+STRING(9B)+$
        STRTRIM(Hp_height[ind[6]],2)+STRING(9B)+$
        STRTRIM(gamma[ind[5]],2)+STRING(9B)+$
        STRTRIM(T_1AU[ind[4]],2)+STRING(9B)+$
        STRTRIM(q[ind[3]],2)+STRING(9B)+STRING(9B)+$
        STRTRIM(r_char[ind[2]],2)+STRING(9B)+$
        STRTRIM(M_H2[ind[1]],2)+STRING(9B)+$
        STRTRIM(mn_chi,2)
      printf,lun3,string_3
      ;      stop
    ENDFOR
    ; Close the min chis squared parms file for the target
    close,lun3
    free_lun, lun3
    
    FOR k = 0, N_ELEMENTS(labwave) - 1 DO BEGIN
      mn_chi = min(chisqr_models_ind[k,*,*,*,*,*,*],location)
      ind = array_indices(chisqr_models_ind[k,*,*,*,*,*,*],location)
      string_3 = STRTRIM(lambda_ref[emission_waves[k]],2)+STRING(9B)+STRING(9B)+$
          STRTRIM(Hp_height[ind[6]],2)+STRING(9B)+$
          STRTRIM(gamma[ind[5]],2)+STRING(9B)+$
          STRTRIM(T_1AU[ind[4]],2)+STRING(9B)+$
          STRTRIM(q[ind[3]],2)+STRING(9B)+STRING(9B)+$
          STRTRIM(r_char[ind[2]],2)+STRING(9B)+$
          STRTRIM(M_H2[ind[1]],2)+STRING(9B)+$
          STRTRIM(mn_chi,2) 
      printf,lun2,string_3
;      stop
    ENDFOR
    ; Close the min chis squared parms file for the target
    close,lun2
    free_lun, lun2
    
    FOR k = 0, N_ELEMENTS(progression[0,*]) - 1 DO BEGIN
      mn_chi = min(chisqr_models_prog[k,*,*,*,*,*,*],location)
      ind = array_indices(chisqr_models_prog[k,*,*,*,*,*,*],location)
      string_3 = prog_name[k]+STRING(9B)+STRING(9B)+$
        STRTRIM(Hp_height[ind[6]],2)+STRING(9B)+$
        STRTRIM(gamma[ind[5]],2)+STRING(9B)+$
        STRTRIM(T_1AU[ind[4]],2)+STRING(9B)+STRING(9B)+$
        STRTRIM(q[ind[3]],2)+STRING(9B)+$
        STRTRIM(r_char[ind[2]],2)+STRING(9B)+$
        STRTRIM(M_H2[ind[1]],2)+STRING(9B)+$
        STRTRIM(mn_chi,2)
      printf,lun1,string_3
    ENDFOR
    
    close,lun1
    free_lun, lun1
;    stop
  ENDFOR
  
  stop




END