;;;;;;; GRID TEST ;;;;;;;
;
; I think that, to accurately model our disks, we need to take "slices" of said
; protoplanetary disk in an (x,y) grid space, where the disk center is located
; at (x,y) = (0,0) and x, y axes are measured in arcseconds away from the center.
;
; This is where the line of sight measurement becomes important and makes sense for
; the disk orientation in the view of the observer.
;

FUNCTION log_grid, r_inner, r_mid, r_outer, distance, sample_inner, sample_outer

  ; Calculate the grid (inner, outer) in arcseconds
  grid_inner = sample_inner / distance
  grid_outer = sample_outer / distance

  ; Calculate the grid size of each region
  range_inner = (r_mid - r_inner)/sample_inner
  range_outer = (r_outer - r_mid)/sample_outer

  ; Define the following grids:
  grid1 = FINDGEN((range_outer+2*grid_outer+1))*grid_outer - (r_outer/distance)   ; r_outer --> r_mid
  grid2 = FINDGEN((range_inner+2*grid_inner+1))*grid_inner - (r_mid/distance)   ; r_mid --> r_inner
  grid3 = FINDGEN((range_inner+2*grid_inner+1))*grid_inner + (r_inner/distance)   ; r_inner --> r_mid
  grid4 = FINDGEN((range_outer+2*grid_outer+1))*grid_outer + (r_mid/distance)   ; r_mid --> r_outer

  ; Smush all the grids together
  grid_array = [grid1, grid2, grid3, grid4]

  RETURN, grid_array

END



FUNCTION define_xygrid, range, step, inc
  ;
  ; Input:
  ;   range - doubles, represent the min/max values in our grid presentation arrays
  ;   step - doubles, step size of the spacing between each grid point
  ; Output:
  ;   grid - array, represents both (x,y) grids in our disk analysis
  ;
  ; *Note: the number of grids and where they start/end will be the same for (x,y).
  ;        (x,y) basically make a square of grids around the disk.
  ;

  ; Create 2 FINDGEN arrays, one from -range --> 0-step, the other from 0+step --> range
  grid = FINDGEN((2*range+2*step)/step)*step*COS(inc*!pi/180.) - range*COS(inc*!pi/180.)
;  grid2 = FINDGEN(range/step)*step
;  grid = [grid1, grid2]




  RETURN, grid

END


;
; Sightline equation to solve:
;
;   sightline = SQRT( d^2 + x^2 + y^2 - 2*d*SQRT(x^2 + y^2)*SIN(i_disk)*COS(ATAN(x/y)) )
;
FUNCTION define_sightline, xgridarray, ygridarray, distance, inc
  ;
  ; Inputs:
  ;   gridarray - array, represents both (x,y) grid space on the disk
  ;   distance - double, distance to the target star
  ;   inc - double, angle at which the observer views the disk
  ; Output:
  ;   sightline - 2d double, represents distance to points on the disk in the (x,y) grid
  ;


  ; Define the sightline 2d matrix
  sightline = DBLARR( N_ELEMENTS(xgridarray),N_ELEMENTS(ygridarray) )

  ; Two FOR loops: One over rows of sightline, and the second over columns of sighline.
  ; Solve for the sightline, as defined above.
  sz = -distance*COS(inc*!pi/180)
  FOR i = 0, N_ELEMENTS(xgridarray) - 1 DO BEGIN     ; represents x grid positions
    FOR j = 0, N_ELEMENTS(ygridarray) - 1 DO BEGIN   ; represents y grid positions
      sy = 0.0
 ;     sx = SQRT(xgridarray[i]^2 + (ygridarray[j])^2)*SIN(ATAN((ygridarray[j])/xgridarray[i]))
      sx = -SQRT(xgridarray[i]^2 + (ygridarray[j])^2)*COS(ATAN((ygridarray[j])/xgridarray[i])) - distance*SIN(inc*!pi/180)
 ;     sightline[i,j] = SQRT( distance^2 + xgridarray[i]^2 + ygridarray[j]^2 + 2*distance*SQRT(xgridarray[i]^2 + ygridarray[j]^2)*SIN(inc*!pi/180)*COS(ATAN(ygridarray[j]/xgridarray[i])) )
      sightline[i,j] = SQRT(sx*sx + sy*sy + sz*sz)
    ENDFOR
  ENDFOR

;  plot,sightline[*,0],ygridarray,xr=[139.9998,140.00002]
;  stop

  RETURN, sightline

END


FUNCTION define_sightline_cylindrical_coords, rH2, phi, distance, inc
  ;
  ; Inputs:
  ;   gridarray - array, represents both (x,y) grid space on the disk
  ;   distance - double, distance to the target star
  ;   inc - double, angle at which the observer views the disk
  ; Output:
  ;   sightline - 2d double, represents distance to points on the disk in the (x,y) grid
  ;

  ; Define the sightline 2d matrix
  sightline = DBLARR( N_ELEMENTS(rH2[*,0]),N_ELEMENTS(rH2[0,*]) )

  ; Two FOR loops: One over rows of sightline, and the second over columns of sighline.
  ; Solve for the sightline, as defined above.
  sz = -distance*COS(inc*!pi/180)
  FOR i = 0, N_ELEMENTS(rH2[*,0]) - 1 DO BEGIN     ; represents x grid positions
    FOR j = 0, N_ELEMENTS(rH2[0,*]) - 1 DO BEGIN   ; represents y grid positions
      sy = (rH2[i,j]/206265.)*COS(phi[i,j]*!pi/180.) ;+ distance*COS(inc*!pi/180)
      ;     sx = SQRT(xgridarray[i]^2 + (ygridarray[j])^2)*SIN(ATAN((ygridarray[j])/xgridarray[i]))
      sx = -(rH2[i,j]/206265.)*SIN(phi[i,j]*!pi/180.) - distance*SIN(inc*!pi/180)
      ;     sightline[i,j] = SQRT( distance^2 + xgridarray[i]^2 + ygridarray[j]^2 + 2*distance*SQRT(xgridarray[i]^2 + ygridarray[j]^2)*SIN(inc*!pi/180)*COS(ATAN(ygridarray[j]/xgridarray[i])) )
      sightline[i,j] = SQRT(sx*sx + sy*sy + sz*sz)
    ENDFOR
  ENDFOR

  ;  plot,sightline[*,0],ygridarray,xr=[139.9998,140.00002]
  ;  stop

  RETURN, sightline

END




FUNCTION define_radius, xgridarray, ygridarray, inc, distance, r_inner, r_outer
  ;
  ; Input:
  ;   gridarray - represents the (x,y) grid positions
  ; Output:
  ;   r_H2 - radial component in the disk that we are sampling in our line of sight
  ;

  ; From each grid point, we can determine the radius of H2 in the disk structure
  ; Define r_H2
  r_H2 = DBLARR( N_ELEMENTS(xgridarray),N_ELEMENTS(ygridarray) )

  ; Loop through the x-grids and store each y-grid column
  FOR i = 0, N_ELEMENTS(xgridarray) - 1 DO BEGIN     ; x grid positions
    FOR j = 0, N_ELEMENTS(ygridarray) - 1 DO BEGIN
      r_dummy = SQRT( xgridarray[i]*xgridarray[i] + ygridarray[j]*ygridarray[j] )*distance  ; AU
      IF (r_dummy LT r_inner) OR (r_dummy GT r_outer) THEN BEGIN
        r_H2[i,j] = 0.0
      ENDIF ELSE BEGIN
        r_H2[i,j] = r_dummy
      ENDELSE
    ENDFOR
  ENDFOR

;  stop
;  plot,xgridarray,r_H2[0,*]
;    stop

  RETURN, r_H2

END




FUNCTION define_phi, xgridarray, ygridarray
  ;
  ; Input:
  ;   gridarray - represents the (x,y) grid positions
  ; Output:
  ;   phi - angular component in the disk that we are sampling in our line of sight
  ;

  ; From each grid point, we can determine the radius of H2 in the disk structure
  ; Define phi
  phi = DBLARR( N_ELEMENTS(xgridarray),N_ELEMENTS(ygridarray) )

  ; Define the half-way point of the grid array (where grid array = 0)
;  xygrid_half = WHERE(xgridarray EQ 0.0)

  ; Loop through the x-grids and store each y-grid column
  FOR i = 0, N_ELEMENTS(xgridarray) - 1 DO BEGIN     ; x grid positions
    FOR j = 0, N_ELEMENTS(ygridarray) - 1 DO BEGIN
    phi[i,j] = ATAN( (ygridarray[j])/xgridarray[i] )*180/!pi      ; degrees

    ; Test to see where where we are in the disk, and how phi needs to be
    ; shifted to fit our idea of the model (phi = 0, 180 at y = 0
    IF ygridarray[j] LT 0.0 THEN BEGIN
      IF xgridarray[i] LT 0.0 THEN BEGIN
        phi[i,j] += 180.0
      ENDIF ELSE BEGIN
        phi[i,j] += 360.0
      ENDELSE
    ENDIF ELSE BEGIN
      IF xgridarray[i] LT 0.0 THEN BEGIN
        phi[i,j] += 180.0
      ENDIF
    ENDELSE

;    IF j LT xygrid_half THEN BEGIN
;      IF i LT xygrid_half THEN phi[i,j] += 90.0
;      ELSE phi[i,j] += 270.0
;    ENDIF
;    ELSE BEGIN
;      IF i LT xygrid_half THEN phi[i,j] += 90.0
;      ELSE phi[i,j] += 270.0
;    ENDELSE

    ENDFOR
  ENDFOR

;  contour,phi,xgridarray,ygridarray
;  print,phi[0,*]
;  stop

  RETURN, phi

END



;
;;;;;;;; NEW SIGHTLINE MEASUREMENT USING r, phi ;;;;;;;;
;
;
; Sightline equation to solve:
;
;   sightline = SQRT( d^2 + x^2 + y^2 - 2*d*SQRT(x^2 + y^2)*SIN(i_disk)*COS(ATAN(x/y)) )
;
FUNCTION define_sightline_r_phi_z, r_H2, phi, distance, inc
  ;
  ; Inputs:
  ;   gridarray - array, represents both (x,y) grid space on the disk
  ;   distance - double, distance to the target star
  ;   inc - double, angle at which the observer views the disk
  ; Output:
  ;   sightline - 2d double, represents distance to points on the disk in the (x,y) grid
  ;

  ; First, define out zgridarray, which will just be xgridarray*COS(inclination)
;  zgridarray = WHERE( index le 950 AND index gt 1050  xgridarray*COS(inc*!pi/180)

  ; Define the sightline 3d matrix (x,y,z)
  ;sightline = DBLARR( N_ELEMENTS(r_H2[*,0]),N_ELEMENTS(r_H2[0,*]) )
  sightline = DBLARR( N_ELEMENTS(r_H2),N_ELEMENTS(phi) )

  ; Three FOR loops: One over rows of sightline, and the second over columns of sighline last over z.
  ; Solve for the sightline, as defined above.
  sz = -distance*COS(inc*!pi/180)
  FOR i = 0, N_ELEMENTS(r_H2) - 1 DO BEGIN     ; represents x grid positions
    FOR j = 0, N_ELEMENTS(phi) - 1 DO BEGIN   ; represents y grid positions
;      FOR k = 0, N_ELEMENTS(zgridarray) - 1 DO BEGIN
        sx = (r_H2[i]/distance)*COS(phi[j]*!pi/180)
;        sx = (r_H2[i,j]/distance)*COS(phi[i,j]*!pi/180)
        sy = -(r_H2[i]/distance)*SIN(phi[j]*!pi/180) - distance*SIN(inc*!pi/180)
 ;       sightline[i,j] = SQRT( distance^2 + gridarray[i]^2 + gridarray[j]^2 + 2*distance*SQRT(gridarray[i]^2 + gridarray[j]^2)*SIN(inc*!pi/180)*COS(ATAN(gridarray[j]/gridarray[i])) )
        sightline[i,j] = SQRT(sx*sx + sy*sy + sz*sz)
;      ENDFOR
    ENDFOR
  ENDFOR

;  plot,sightline[0,*],ygridarray
;  stop

  RETURN, sightline

END


;
;
;;;;;;;;;;;;;;;;;;; NEED NEW GRID MEASUREMENTS ;;;;;;;;;;;;;;;;;;;;
;
FUNCTION grid_radial_gas, rH2
  ; We need to calculate the grid size of the radial parcel of H2 gas
  ; we are observing, which will
  ;

  ; Make an empty array for the radial grid points
  radial_grid = DBLARR( N_ELEMENTS(rH2[*,0]),N_ELEMENTS(rH2[0,*]) )

  ; Loop through rad and phi positions on the disk, and determine the grid size
  ; of our grids, defining it by the angular distance between r_d in the radius
  FOR i = 0, N_ELEMENTS(rH2[*,0]) - 1 DO BEGIN
    FOR J = 0, N_ELEMENTS(rH2[0,*]) - 1 DO BEGIN
      ; If the first element, you have to determine a little differently. Otherwise
      ; do [i] - [i-1]
      IF i EQ 0 THEN BEGIN
        radial_grid[i,j] = ABS(rH2[i,j] - rH2[i+1,j])
      ENDIF ELSE BEGIN
        radial_grid[i,j] = ABS(rH2[i,j] - rH2[i-1,j])
      ENDELSE
    ENDFOR
  ENDFOR

  RETURN, radial_grid

END






FUNCTION define_radial_velocity, rH2, Mstar
  ;
  ; Input:
  ;   rH2 - 2d array, radial positions in the disk based on (x,y) location on the grid
  ;   Mstar - double, stellar mass of target
  ; Output:
  ;   v_k - 2d array, represents the radial velocity of the H2 gas parcel at a given radius
  ;         in the disk.
  ;

  ; Define constants/conversions
  G = 6.67259E-8           ; Gravitational constant [cm+3 g-1 s-2]
  AUtoCM = 1.496E+13       ; [AU --> cm] - conversion for r_obs from cm to AU
  CMtoKM = 1.0E-5          ; [cm --> km]
  MSUNtoG = 1.99E+33       ; [M_solar --> g] - conversion for M_star

  ; Define v_k array
  ;v_k = DBLARR( N_ELEMENTS(rH2[*,0]),N_ELEMENTS(rH2[0,*]) )
  v_k = DBLARR( N_ELEMENTS(rH2) )

  ; FOR loop through the rows (x-pixels in the grid) and add the elements of v_k
  ; along the y pixels
;  FOR i = 0, N_ELEMENTS(rH2[*,0]) - 1 DO BEGIN
;    FOR j = 0, N_ELEMENTS(rH2[0,*]) - 1 DO BEGIN
;      v_k[i,j] = SQRT(G*Mstar*MSUNtoG / (rH2[i,j]*AUtoCM)) * CMtoKM
;      ;IF v_k[i,j] GT 200.0 THEN v_k[i,j] = 0.0
;    ENDFOR
;  ENDFOR

  FOR i = 0, N_ELEMENTS(rH2) - 1 DO BEGIN
    v_k[i] = SQRT(G*Mstar*MSUNtoG / (rH2[i]*AUtoCM)) * CMtoKM
    ;IF v_k[i,j] GT 200.0 THEN v_k[i,j] = 0.0
  ENDFOR


;  stop

  RETURN, v_k

END



FUNCTION vrad_to_vobs, r_H2, M_star, i_disk, phi
  ;
  ; Input:
  ;   v_rad - 2d array, radial velocity in disk over (x,y) grid
  ;   inc - double, inclination angle (viewing angle) of the disk in our sight line
  ;   phi - 2d array, angular velocity in disk over (x,y) grid
  ; Output:
  ;   vobs - 2d array, observed velocity of the gas in the disk at (x,y)
  ;

  ; Define v_k
  v_k = define_radial_velocity(r_H2,M_star)
;  stop

  ; Define vobs array
  vobs = DBLARR( N_ELEMENTS(r_H2[*,0]),N_ELEMENTS(r_H2[0,*]) )
  vobs = DBLARR( N_ELEMENTS(r_H2),N_ELEMENTS(phi) )


  ; FOR loop over x-pixels and calculate over y-pixels the resulting shift in
  ; velocity from our observing frame
;  FOR i = 0, N_ELEMENTS(r_H2[*,0]) - 1 DO BEGIN
;    FOR j = 0, N_ELEMENTS(r_H2[0,*]) - 1 DO BEGIN
;      vobs[i,j] = v_k[i,j]*SIN(i_disk*!pi/180.)*COS(phi[i,j]*!pi/180.)
;    ENDFOR
;  ENDFOR

  FOR i = 0, N_ELEMENTS(r_H2) - 1 DO BEGIN
    FOR j = 0, N_ELEMENTS(phi) - 1 DO BEGIN
      vobs[i,j] = v_k[i]*SIN(i_disk*!pi/180.)*SIN(phi[j]*!pi/180.)
    ENDFOR
  ENDFOR

;  FOR i = 0, N_ELEMENTS(vobs[*,0]) - 1 DO BEGIN
;    FOR j = 0, N_ELEMENTS(vobs[0,*]) - 1 DO BEGIN
;      IF vobs[i,j] LT 0.0 THEN vobs_neg += 1
;      IF vobs[i,j] GE 0.0 THEN vobs_pos += 1
;    ENDFOR
;  ENDFOR

;  print,'Number of positive observed velocity points: ', vobs_pos
;  print,'Number of negative observed velocity points: ', vobs_neg
;  plot,phi[*,0],vobs[*,0]
;  stop
;  stop


  RETURN, vobs

END




FUNCTION define_radial_temperature, r_H2, T_1AU, q
  ;
  ; Inputs:
  ;   r_H2 - 2d array, describes radius looking into disk by observed in (x,y) grid
  ;   T_1AU - double, temperature defined at 1AU
  ;   q - double, temperature gradient of our temperature power law profile
  ; Output:
  ;   T_H2 - 2d array, describes the temperature profile of the disk in the (x,y) grid
  ;
  ; *Note: The temperature profile is described as:
  ;      T_H2 = T_1AU * ( rH2[i,*] / r (1AU) )^q
  ;

  ; Define the H2 temperature profile array
  T_H2 = DBLARR( N_ELEMENTS(r_H2) )

  ; FOR loop for the x-pixels and solve for the temperature profile along the y-pixels
  FOR i = 0, N_ELEMENTS(r_H2[*,0]) - 1 DO BEGIN
;    FOR j = 0, N_ELEMENTS(r_H2[0,*]) - 1 DO BEGIN
      T_H2[i] = T_1AU*(r_H2[i])^q
      IF (T_H2[i] GT 5000.0) OR (T_H2[i] LT 1000.0) THEN T_H2[i] = 0.0
;    ENDFOR
  ENDFOR

;  plot,r_H2[800,*],T_H2[800,*]
;  stop
  RETURN, T_H2

END



FUNCTION define_doppler_broadening, T_H2, v_turb
  ;
  ; Input:
  ;   T_H2 - 2d array, describes the radial temperature profile in the line of sight on the
  ;         (x,y) grid
  ;   v_turb - double, turblent velocity of the gas in the disk (assumed constant at all radii)
  ; Output:
  ;   dop_v - 2d array, the doppler broadening term that describes the width of the line
  ;           line profile of the emission line
  ;

  ; Define any constants/conversions
  k_b = 1.38E-16        ; Boltzmann constant [erg K-1]
  m_H2 = 2.*1.67E-24    ; mean molecular weight of H2, = 2* m_H [g]
  KMtoCM = 1.0E+5       ; [km --> cm] - conversion for del_v

  ; Define the dop_v array
  dop_v = DBLARR( N_ELEMENTS(T_H2[*,0]),N_ELEMENTS(T_H2[0,*]) )

  ; FOR loop over T_H2 in x-pixel and solve for dop_v in y-pixel space
  FOR i = 0, N_ELEMENTS(T_H2[*,0]) - 1 DO BEGIN
    FOR j = 0, N_ELEMENTS(T_H2[0,*]) - 1 DO BEGIN
      dop_v[i,j] = SQRT( ( (2.0*k_b*T_H2[i,j]) / m_H2) + (v_turb*KMtoCM)^2 ) / KMtoCM   ; km s-1
    ENDFOR
  ENDFOR

  RETURN, dop_v

END




FUNCTION define_lineprofile, lambda, TrH2, v_turb, v_doppler, vobs
  ;
  ; Inputs:
  ;   lambda - double, rest wavelength of the H2 emission line
  ;   delV - 2d array, Doppler Broadening term found in function define_doppler_broadening
  ;   v_doppler - double, offset in central position of the emission line from the
  ;               doppler speed of the star/disk
  ;   vobs - 2d array, observed velocity of the gas parcels measured in the (x,y) plane
  ; Output:
  ;   line_prof - 2d array, the resulting emission feature of the H2 gas given the
  ;                  (x,y) position in the disk.
  ;                  Will probably look at this over x values at a given y pixel.
  ;
  ; *Note: Line profile is defined as:
  ;                         lambda          ( [v_doppler - vobs)^2 )
  ;      Line profile = -------------- * EXP( -------------------- )
  ;                     SQRT(!pi)*delV      (        delV^2        )

  ; Define constant/conversions to be used
  KMtoCM = 1.0E+5          ; [km --> cm] - conversion for delV (in denominator)
  ANGtoCM = 1.0E-8         ; [Ang --> cm] - conversion for lambda_0

  ; Define delV
  delV = define_Doppler_broadening(TrH2,v_turb)

  ; Make a v_obs tolerance level around v_obs values that have steep transitions
  ; So, say if dv = 5 km s-1, then skip that calculation?
  dv = 20.0     ; km s-1

  ; Define the line_profile array
  line_prof = DBLARR( N_ELEMENTS(vobs[*,0]),N_ELEMENTS(vobs[0,*]) )

  ; FOR loop through the x-pixels and solve for all y-pixels at once
  FOR i = 0, N_ELEMENTS(line_prof[*,0]) - 1 DO BEGIN
    FOR j = 0, N_ELEMENTS(line_prof[0,*]) - 1 DO BEGIN
      ; First, make sure delV =/= 0
      IF TrH2[i] EQ 0.0 THEN BEGIN
        line_prof[i,j] = 0.0
      ENDIF ELSE BEGIN
      ; short circuit evaluate out (i-1) if i == 0
;      IF i EQ 0 OR i GT 0 AND ABS(vobs[i,j] - vobs[i-1,j]) LT dv THEN BEGIN
        line_prof[i,j] = ( lambda*ANGtoCM / ( SQRT(!pi)*delV[i]*KMtoCM ) )*EXP( -(v_doppler - vobs[i,j])^2 / (delV[i]*delV[i]) )
;      ENDIF ELSE BEGIN
;        line_prof[i,j] = 0.0
      ENDELSE
    ENDFOR
  ENDFOR

;  plot,vobs[*,1000],line_prof[*,1000],xr=[-40,40]
;  stop

  RETURN, line_prof

END





;;;;;;;;;;;;;;;; COLLAPSE LINE PROFILE INTO 1D TOTAL LINE PROFILE ;;;;;;;;;;;;;;;;;;
;
;
FUNCTION define_collapsed_line_profile, binned_vobs, bin_factor, lineprof, vobs
  ;
  ; Inputs:
  ;   lineprof - 2d array (representing the x,y grid points on our disk), presents the emission profile of the H2
  ;               gas located in the section of the disk sampled over the grid points
  ;   vobs - 2d array, velocity of the gas parcel in the disk sampled at x,y
  ; Output:
  ;   total_line_profile - 1d array, sampled only over steps of 1 km s-1 in v_obs, to give us the total line
  ;                        profile from all the emitting regions of the disk
  ;
  stop
  ; Define the binned v_obs array, empty total_line_profile
  total_line_profile = DBLARR( N_ELEMENTS(binned_vobs),N_ELEMENTS(lineprof[0,0,*]) )

  ; Collapse v_obs, line_prof down into 1d arrays for ease
;  vobs_1d = REFORM(vobs,1,N_ELEMENTS(vobs[0,*])*N_ELEMENTS(vobs[*,0]))
;  lineprof_1d = REFORM(lineprof,1,N_ELEMENTS(lineprof[0,*])*N_ELEMENTS(lineprof[*,0]))

  nbins = N_ELEMENTS(total_line_profile) - 1
  binmin = min(binned_vobs)
  binmax = max(binned_vobs)
  binrange = (binmax - binmin)/nbins
  FOR i = 0, N_ELEMENTS(vobs[*,0]) - 1 DO BEGIN ;
    FOR p = 0, N_ELEMENTS(lineprof[0,0,*]) - 1 DO BEGIN
    FOR j = 0, N_ELEMENTS(vobs[0,*]) - 1 DO BEGIN
      IF vobs[i,j] GE binmin OR vobs[i,j] LE binmax THEN BEGIN
        ; Store v_index
;        v_index = WHERE( FIX(vobs[i,j]) EQ binned_vobs )
        v_index = WHERE( FIX(vobs[i,j]) EQ binned_vobs )
        ;IF v_index EQ -1 THEN v_index = 300
        IF v_index EQ -1 THEN total_line_profile[300] += 0.0
        IF FINITE(lineprof[i,j,p]) EQ 1 THEN BEGIN   ; Returns only real, finite numbers
        ;  IF vobs_1d[i] LE binmax AND vobs_1d[i] GE binmin THEN BEGIN
          total_line_profile[v_index,p] += lineprof[i,j,p]
        ENDIF ELSE total_line_profile[v_index] += 0.0
      ENDIF
    ENDFOR
    ENDFOR
  ENDFOR
  stop


  ; FOR loop: loop through the i elements of vobs to look for a match in binned_vobs (k)
  ; If there is a match, then use that v_obs index to add the line_profile element at that index
  ; to the total_line_profile array
;  FOR i = 0, N_ELEMENTS(vobs_1d) - 1 DO BEGIN           ; for data accumulated thus far in our grid
;    FOR k = 0, N_ELEMENTS(binned_vobs) - 1 DO BEGIN     ; for binned data to come
;      ; Test to see if i is in k
;      IF UINT(vobs_1d[i]) EQ binned_vobs[k] THEN BEGIN
;      ;IF vobs_1d[i] GE binned_vobs[k]-(bin_factor/2.) AND vobs_1d[i] LE binned_vobs[k]+(bin_factor/2.) THEN BEGIN
;        total_line_profile[k] += lineprof[i]
;        BREAK
;      ENDIF
;    ENDFOR
;  ENDFOR

;  plot,binned_vobs,total_line_profile,xr=[-30,30]
;  stop


  RETURN, total_line_profile

END





;pro grid_disks, r_grid, radius_H2, phi_grid, phi_disk, distance, restwave, LyAwave, M_star, i_disk, T_1AU, grad, v_turbulent, v_Doppler, sight_line, v_obs, T_rH2, H2_line_profile
pro grid_disks_v2, r_inner, r_mid, r_outer, AUsize_inner, AUsize_outer, phi_density, distance, restwave, LyAwave, M_star, i_disk, T_1AU, grad, v_turbulent, v_Doppler, r_grid, phi_grid, sight_line, v_obs, T_rH2, H2_line_profile
;pro grid_disks


;
; INPUTS:
;   distance - double (parsecs)
;   restwave - double (Angstroms)
;   M_star - double (Solar masses)
;   i_disk - double (degrees)
;   grad - negative double (no units)
;   v_turbulent - double (km s-1)
;   v_Doppler - double (km s-1)
;
; OUTPUTS:
;   radius_H2 - 2d double matrix (AU)
;   phi_disk - 2d double matrix (degrees)
;   sight_line - 2d double matrix (parsecs)
;   v_obs - 2d double matrix (km s-1)
;   T_rH2 - 2d double matrix (Kelvin)
;   H2_line_profile - 2d double matrix (no units)
;
;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Call xgrid, ygrid, r_H2, and phi_disk from the grid_from__radius.pro program,
; which determines the x,y grid and r,phi locations on the disk of the H2
grid_from_radius, i_disk, distance, r_inner, r_mid, r_outer, AUsize_inner, AUsize_outer, phi_density, xgrid, ygrid, r_grid, phi_grid, radius_H2, phi_disk
;
r_grid = r_grid*distance
phi_disk = phi_disk*180/!pi
phi_grid = phi_grid*180/!pi

sight_line = define_sightline_r_phi_z(r_grid, phi_grid, distance, i_disk)
;stop

; Call conversion from v_k to v_obs
v_obs = vrad_to_vobs(r_grid, M_star, i_disk, phi_grid)
;stop

; Call the radial temperature profile of the disk.
T_rH2 = define_radial_temperature(r_grid, T_1AU, grad)


RETURN
END
