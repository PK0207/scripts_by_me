;; Model Line Profiles for H2 disks around T Tauri Stars
;;
;;
;; Keri Hoadley
;; Outline - 24 January 2013
;; Assumptions + Approximations - 11 April 2013
;; Parameters and Limits of Parameters - October/November 2013
;;
;; Revision 2 - 19 September 2014
;;
;
; Purpose
; -------
; We want to create models that characterize the radial/rotational disk structure
; of the molecular gas in protoplanetary disks.
; Therefore, we will have a set of free parameters (some more free than others) that
; should be adequate to modeling a disk given the physical conditions in the disk.
;
;
; Outline
; -------
; Based on modeling done by Rosenfeld et al 2012 (Section 3).
;
;
;
;
;
; Assumptions and Approximations
; ------------------------------
; 1. LTE is okay:
;     - H2 emission energy levels are populated according to LTE
;       --> low enough energy and high optical depth associated with H2
;     - in the inner regions of protoplanetary disks, the optical depths are very large
;     - free electrons (from photodissociation of atoms/dust grains) allows assumption to hold
; 2. Keplerian rotation of the disk material
; 3. Stellar mass is treated as a point mass of M_star
; 4. Geometrically thin disk at all radii (or at least radii of H2 probed)
; 5. Gas structure:
;     - vertically isothermal (no large temperature gradients in the vertical axis of the disk)
;     - hydrostatic equilibrium (dP/dr = -rho*g*height)
; 6. The H2 emission energy levels are populated according to LTE
;
;
;
; Parameters & Limits of Parameters
; ---------------------------------
; Parameters:
;   1. M_H2   = total PARTICIPATING (observed) mass of the H2 in the disk
;   2. M_star = total stellar mass, treated as a point source
;   3. R_char = characteristic radium of the H2 <-- observed radius of the molecular emission?
;   4. T_R10  = temperature of the disk at R = 10 AU
;   5. q_T    = temperature gradient of T_R10
;   6. v_turb = turbulence in the disk, which is treated as constant at a given R
;   7. i_disk = inclination angle of the disk with respect to the observer
;
; Limitations:
; 1. M_H2   ~ [0.001, 0.016, 0.0001 steps] M_solar
; 2. M_star ~ [0.1, 1.6, 0.1 steps] <-- take range from ranges of my target stars +/- errors in the masses
; 3. R_char ~ [0.03, 15, 0.01 steps?] AU
; 4. T_R10  ~ [1500, 3000, 1 step?]K <-- find source (b/c Kevin told me, but I should find it)
; 5. q_T    ~ [-5, 5, 0.1 steps?] <-- Rosenfeld paper?
; 6. v_turb ~ [0, 15, 0.1 steps?]km/s   <-- again, find literature for this assumption (sam reason as 4.)
; 7. i_disk ~ [0, 90, 1 steps?] degrees worst case scenario, but maybe we can narrow this range down?
;
;
;
; Revision 2
; ----------
; 1. Revise cross section calculation
;       - take out all line profile values
;       - take out exponential function - this gets considered in the flux calculation only
;       - only define cross section = intrinsic cross section
;       - in tau: take out line profile function
;                 include Matt's new tau_eff values
;
; 2. Revise flux calculation
;       F(H2) ~ F(LyA)(1 - e(-tau))*Bul*(d^2)
;
; 3. Calculate all progressions at one time:
;       - make extra dimension in arrays to hold for progression absorption values
;       - make extra dimension in final flux array for all 12 emission lines considered (by target)
;
;
;




;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SYNTHETIC H2 SPECTRAL MODELS - LINE PROFILE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;
;
; Take a range of velocities, from the first and last velocities in the arrays made from the H2 data,
; this should range [-200,+200] km s-1.
; What I will need to give the model:
;   - speed of light --> c
;   - Doppler velocity --> maybe an array of values around the max velocities of each of the disks? [v_doppler]
;   - observational velocity --> describes the range of velocities explored in velocity space of our H2 data [v_obs]
;   - molecular weight of H2, ~ 2*mH [m_H2]
;   - Boltzmann constant [k_b]
;   - stellar mass [M_star]
;   - turbulent velocity [v_turb]
;   - rest wavelength [lambda_0]
;   - radius of the emitting H2 [r_obs]
;   - temperature gradient at T = 1 AU [q]
;   - temperature that will describe the molecular temperature at 1 AU as constant [T_1AU]
;
;;;;;;;;;;;;;;; *********** ALL TAKEN CARE OF IN grid_disks.pro ************** ;;;;;;;;;;;;;;;;;;;
;;
;;
; Call:
; grid_disks, distance, restwave, M_star, i_disk, T_1AU, grad, v_turbulent, v_Doppler, radius_H2, phi_disk, sight_line, v_obs, T_rH2, H2_line_profile
;
; Inputs:
;   distance, restwave, M_star, i_disk, T_1AU, grad, v_turbulent, v_Doppler
; Outputs:
;   radius_H2, phi_disk, sight_line, v_obs, T_rH2, H2_line_profile
;




;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; SYNTHETIC H2 SPECTRA MODELS - OPACITY FUNCTION (kappa) ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;
;;
;
;
;
; Now that we have complete (Doppler) line profiles are created, we need to use this as a means of determining
; the opacity of the H2
;
;


; A function used in both parts of defining the opacity of the emission line
; This is the exponential expression described by the wavelength of the line and
; the temperature over the gas (as a funtion of radius on the disk).
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - exp_temp ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;FUNCTION exp_temp, T_H2, lambda_0
;  ;
;  ; Input parameters:
;  ;     T_1AU - 2d array, temperature profile of the disk as a function of radius
;  ;     lambda_0 - double, the rest wavelength of the emission H2 line
;  ;
;  ; Ouput:
;  ;     exp_T - 2d array, array of the exponential function, dependent on emission line energy and
;  ;                temperature in the disk.
;  ;
;
;  ; Define any constants to be used
;  k_B = 1.38E-16    ; Boltzmann constant [erg K-1]
;  c = 3.0E10        ; speed of light [cm s-1]
;  h = 6.626E-27     ; Planck constant [erg s]
;
;  ; Conversions:
;  ANGtoCM = 1.0E-8   ; [Ang --> cm] for wavelength unit conversion
;
;  exp_T = DBLARR( N_ELEMENTS(T_H2[*,0]),N_ELEMENTS(T_H2[0,*]) )
;
;  ; Calculate the exponential relation for the energy of the emission line as a function of temperature
;  FOR i = 0, N_ELEMENTS( T_H2[*,0] ) - 1 DO BEGIN
;    exp_T[i,*] = EXP( -(h*c)/((lambda_0*ANGtoCM)*k_B*T_H2[i,*]) )
;  ENDFOR
;
;  RETURN, exp_T
;END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - exp_temp ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
FUNCTION trans_temp, T_H2, trans_energy
  ;
  ; Input parameters:
  ;     T_1AU - 1d array, temperature profile of the disk as a function of radius
  ;     trans_energy - 1D array, the ground state energy of available absorbing H2
  ;
  ; Ouput:
  ;     exp_T - 2d array, array of the exponential function, dependent on emission line energy and
  ;                temperature in the disk.
  ;

  ; Define any constants to be used
  k_B = 1.38E-16    ; Boltzmann constant [erg K-1]
  c = 3.0E10        ; speed of light [cm s-1]
  h = 6.626E-27     ; Planck constant [erg s]

  ; Conversions:
  ANGtoCM = 1.0E-8   ; [Ang --> cm] for wavelength unit conversion

  ; Calculate the exponential relation for the energy of the emission line as a function of temperature
;  exp_trans = EXP( -(h*c)/((lambda_LyA*ANGtoCM)*k_B*T_H2) )
  exp_trans = DBLARR( N_ELEMENTS(T_H2), N_ELEMENTS(trans_energy) )

  FOR i = 0, N_ELEMENTS(trans_energy) - 1 DO BEGIN
      exp_trans[*,i] = EXP( (-1)*trans_energy[i] / (k_B*T_H2) )
  ENDFOR

  RETURN, exp_trans
END




;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - PARTITION FUNCTION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; With the energy distribution of H2 in the disk characterized at each grid point,
; we construct the partition function as the sum of all the energy distribution values
; in the disk. This helps normalize the energy distribution at a given grid point
; to understand how the number density of H2 is laid out in the disk.
;
; THIS REPRESENTS THE H2 POPULATION BEFORE LyA PUMPING!!!!
;
FUNCTION partition_function, T_H2, J_rot, T_H2_range, Z_TH2
  ;
  ; INPUTS:
  ;   exp_temp - 2d array, the Boltzmann energy distribution of gas at each grid point
  ;   ang_mom - int, the angular momentum of the H2 species, described by the final transition state
  ;             for the emission feature.
  ;
  ; OUTPUT:
  ;   part_fnct - 2D array [TH2, progH2], describes the sum of all the energy distribution values at each grid point.
  ;

    ; Calculate the "partition function", or H2 population distribtuion, for the line, given T_rH2
  part_fnct = DBLARR( N_ELEMENTS(T_H2),N_ELEMENTS(J_rot) )

  print,SYSTIME()

  ; Index the appropriate J level
 ; j_index = WHERE(J_rot EQ rotation_level)

  ; For loop: Loop through T_H2 to find where T_H2_range = T_H2[i,k]
  ; The store Z_TH2 into part_fnct[i,k]
  FOR i = 0, N_ELEMENTS(T_H2) - 1 DO BEGIN
;    FOR k = 0, N_ELEMENTS(T_H2[0,*]) - 1 DO BEGIN
      IF T_H2[i] EQ 0.0 THEN BEGIN
        part_fnct[i] = 0.0
      ENDIF ELSE BEGIN
        t_index = WHERE(T_H2_range EQ UINT(T_H2[i]))
        FOR j = 0, N_ELEMENTS(J_rot) - 1 DO BEGIN
 ;         j_index = J_rot[j]
          part_fnct[i,j] = Z_TH2[t_index,j]
        ENDFOR
      ENDELSE
;    ENDFOR
  ENDFOR

  print,SYSTIME()
;  plot,T_H2,part_fnct
;  stop


  RETURN, part_fnct

END





;;
;; A) Cross section of H2 as a function of line profile and internal properties of H2
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - sigma_0 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
; This inital cross section depends only on when the wavelength and thus the intrinsic
; properties of the H2 transitions from[v,J]' to [v,J]".
; Therefore, this will be called in the lambda_0 FOR loop because this will not change
; when we things like M_star, r_H2, or an of the other parameters.
;
FUNCTION sigma_0, lambda_0, j1, j2, A_ul_total
  ;
  ; Parameter Inputs:
  ;     lambda_0 - 1D array [4 elements], LyA pumping wavelength
  ;     j1, j2 - 1D array [4 elements], ground and excited state J values
  ;     A_ul_total - 1D array [4 elements], sum of all Einstein coefficients from a specific progression
  ; Output:
  ;     sig0 - 1D array [4 elements], internal cross section of the absorbing H2
  ;

  ; Conversion of untis
  ANGtoCM = 1.0E-8      ; [Ang --> cm] conversion for wavelength
  c = 3.E10

  ; Define the g_l and g_l+1 values
  g_l0 = DBLARR(N_ELEMENTS(j1))
  g_l1 = DBLARR(N_ELEMENTS(j1))
  FOR i = 0, N_ELEMENTS(j1) -1 DO BEGIN
   IF j1[i] MOD 2 EQ 0 THEN g_l0[i] = (2.*j1[i] + 1.) ELSE g_l0[i] = 3.*(2.*j1[i] + 1.)
   IF j2[i] MOD 2 EQ 0 THEN g_l1[i] = (2.*j2[i] + 1.) ELSE g_l1[i] = 3.*(2.*j2[i] + 1.)
  ENDFOR

  ; Define sig0
  sig0 = ( (lambda_0*ANGtoCM)^3 / (8.0 * !PI * c) )*( g_l1 / g_l0 )*A_ul_total

  ; Return the initial cross section value
  RETURN, sig0

END




;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - Disk Surface Density ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
FUNCTION surface_density, r_H2, r_c, density_grad
  ;
  ; Parameter Inputs:
  ;       r_c - double (variable), characteristic radial location of the H2 gas
  ;       r_H2 - 1d array, radial location of H2 in the disk
  ;       density_grad - double (variable), surface density gradient that describes the
  ;                      power-law and exponential decay functions
  ; Output:
  ;       surface_density_H2 - 1d array, describes the density of H2 as a function of the radius in the disk
  ;

  ; Calculate the power-law relation between r_c and r_H2
  power_law = (r_H2 / r_c)^(-1.0*density_grad)

  ; Calculate the exponential relation of the surface density
  exp_law = EXP( -1.0*(r_H2/r_c)^(2.0 - density_grad) )

  surface_density_H2 = power_law*exp_law

  ; Calculate the approximate surface density (no constants)
 ; surface_density_H2 = DBLARR( N_ELEMENTS(r_H2) )


;  FOR i = 0, N_ELEMENTS(r_H2[*,0]) - 1 DO BEGIN
;    FOR j = 0, N_ELEMENTS(r_H2[0,*]) - 1 DO BEGIN
;  FOR i = 0, N_ELEMENTS(r_H2) - 1 DO BEGIN
;    ; Calculate the power-law relation between r_c and r_H2
;    power_law = (r_H2[i] / r_c)^(-1*density_grad)
;
;    ; Calculate the exponential relation of the surface density
;    exp_law = EXP( -(r_H2[i]/r_c)^(2 - density_grad) )
;
;    surface_density_H2[i] = power_law*exp_law
;  ENDFOR

;  stop

  ; Return the surface density relation to find the number density of H2 later.
  RETURN, surface_density_H2

END




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - Scale Height ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
FUNCTION scale_height, T_AU, r_H2, M_star, phi
  ;
  ; Parameters Input:
  ;     T_AU  - 1d array, based on the position in the disk r_H2
  ;     r_H2   - 1d array, position of H2 in the disk
  ;     M_star - double, stellar mass
  ; Return:
  ;     H_z - 1d array
  ;

  ; Define the scale height function as an array
;  H_z = DBLARR( N_ELEMENTS(r_H2) )

  ; Define constants used to define the scale height
  k_B = 1.38E-16          ; [erg K-1]
  m_H2 = 2.33 * 1.67E-24   ; [g] - mu = 2.33 = mean mole. weight for circumstellar material (Ruden & Pollack 1991)
  G = 6.67259E-8          ; [cm+3 g-1 s-2]

  ; Define conversions
  AUtoCM = 1.496E+13       ; [AU --> cm] - conversion for r_H2 from cm to AU
  MSUNtoG = 1.99E+33       ; [M_solar --> g] - conversion for M_star

  thermal = k_B*T_AU / m_H2
  kinetic = ((r_H2) / (G*M_star)) * (AUtoCM/MSUNtoG)
  H_z = SQRT( thermal * kinetic )*r_H2
;  stop
  ; Loop thru the input array lengths to calculate the scale height as a function of radius
;  FOR i = 0, N_ELEMENTS(r_H2[*,0]) - 1 DO BEGIN
;    FOR j = 0, N_ELEMENTS(r_H2[0,*]) - 1 DO BEGIN

;  FOR i = 0, N_ELEMENTS(r_H2) - 1 DO BEGIN
;    thermal = k_B*T_AU[i] / m_H2
;;    kinetic = ((r_H2[i])^3)*((AUtoCM)^3) / (G*M_star*MSUNtoG)
;    kinetic = ((r_H2[i]) / (G*M_star)) * (AUtoCM/MSUNtoG)
;;    H_z[i] = SQRT( (k_B*T_AU[i]*((r_H2[i]*AUtoCM)*(r_H2[i]*AUtoCM)*(r_H2[i]*AUtoCM)))/(m_H2*G*M_star*MSUNtoG) )/AUtoCM  ; Scale height is in cm
;    H_z[i] = SQRT( thermal * kinetic )*r_H2[i];/AUtoCM  ; Scale height is in cm
;  ENDFOR

  ; Finally, convert H_z from cm to AU
  ;H_z = H_z/AUtoCM
;  stop
  ; Return the scale height to use to determine the number density of H2 in the disk.
  RETURN, H_z

END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - Exponential Disk Density ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Describes the disk density of H2 as a function of height in the disk, z
; We will define z within the function, for a given Hz value.
; Since Hz will be an array of the same size as r_H2, z will not be saved (per say), and we will
; define the z vector for an individual Hz value.
; Then, we will create a disk density profile that will be a 2D array, with the rows representing the
; radius of H2 in the disk and the columns representing the disk height position, z.
;
FUNCTION exp_disk_density, Hz, disk_height
  ;
  ; Parameter inputs:
  ;   Hz - 1d array, scale height of the disk at a given value of r_H2
  ;   disk_height - int, describes where the outer-most emitting H2 is located in z
  ;
  ; Output parameter:
  ;   rho_disk - 3d array, describes the density of H2 in the disk at a given position in
  ;              the radial and vertical location of the protoplanetary disk.
  ;

  ; First order of business: convert AU to cm
  AUtoCM = 1.496E+13       ; [AU --> cm] - conversion for r_H2 from cm to AU
  z_step = 10              ; number of steps to take in z until we reach 1 scale height
                           ; below the outer-most disk atmosphere


  ; Define the density 2d array that will be N_ELEMENTS(Hz) x N_ELEMENTS(Hz)
  rho_disk = DBLARR( z_step+1 )

  ; First, define z by the same number of elements as Hz, but split up evenly between 0 --> +6 Hz.
;  FOR i = 0, N_ELEMENTS(Hz) - 1 DO BEGIN
    ; Starting from the upper-most layer (z=0), calculate the exp. drop-off in density
    ; as we descend into the disk (-z/z_step)*Hz
    FOR z = 0, z_step DO BEGIN
      height = ( disk_height - (FLOAT(z)/z_step) )^2   ; Actually: (d_h*Hz - (z/zstep)*Hz / Hz)^2, but Hz terms cancel
      rho_disk[z] = EXP( -0.5* height )
    ENDFOR
;  ENDFOR

  ; We solved for the disk density of H2, so now we can use it in our program.
  ; Returns only rho[z] (1D) because relationship between all z/Hz at all radii remains the same.
  RETURN, rho_disk


END


;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - CROSS SECTION of H2 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Use the line profile array, sigma_0 of the intrinsci properties of the H2, and the thermal distribution of
; energy of the H2 line to determine the cross section of the H2 in the disk to then to determine the
; opacity of the H2 in the disk as a function of radius and vertical position in the disk.
;
FUNCTION cross_section, sigma_0
  ;
  ; Parameters Called:
  ;   line_profile - 1d array, describes the doppler broadening of the H2 emission line (as a function of r_H2)
  ;   sigma_0 - double, the intrinsic cross section of the H2 species, dependent only on properties of the
  ;             H2 transition
  ;   thermal_distr - 1d array, the energy distribution of the H2 as a function of r_H2
  ;
  ; Output parameter:
  ;   xsection - 1d array, the thermal cross section of the H2 emission line as a function of r_H2
  ;

  ; Define the x-section array
  ;xsection = DBLARR(N_ELEMENTS(line_profile[*,0]),N_ELEMENTS(line_profile[0,*]))
  xsection = sigma_0

  ; Loop through each array value of line_profile and thermal_distr and calculate the xsection
;  FOR i = 0, N_ELEMENTS( thermal_distr[*,0] ) - 1 DO BEGIN
;    FOR j = 0, N_ELEMENTS( thermal_distr[0,*] ) - 1 DO BEGIN
;    ;  xsection[i,j] = sigma_0*(1.0 - thermal_distr[i,j]);*line_profile[i,j]
;      xsection[i,j] = sigma_0;*line_profile[i]
;    ENDFOR
;  ENDFOR
;
  ; THAT WAS EASY!
  RETURN, xsection

END



;;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - H2 Number Density ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Last large step before calculating the opacity of H2 in the disk as a function of radius. Above, we
; calculate the thermal cross sections of the H2 in the disk, so next we need to know the number density
; of H2 at our given transition to understand how the two components together create the opacity
; of gas in the disk.
;
FUNCTION H2_number_density, M_H2, Hz, r_c, j_rot, rho_surface, rho_disk, density_grad, thermal_distr, partition_fnct
  ;
  ; Parameters Input:
  ;   M_H2 - double, estimate of the total mass of H2 emitting in the disk. [units: M_solar]
  ;   H_z - 1d array, describes the scale height of the disk at different radii [units: AU]
  ;   r_c - double, our estimated characteristic radius of H2 in the disk [units: AU]
  ;   l_0 - double, the angular momentum of the final(?) state of the H2 transition
  ;   total_protons - double, essentially 'Z'
  ;   rho_surface - 1d array surface density at given radii in the disk of H2 gas [unitless]
  ;   rho_disk - 2d array, the density of gas in the disk as a function of r_H2 and vertical height in the disk
  ;              [unitless]
  ;   thermal_distr - 1d array, the themal energy distribution of H2 gas in the disk, a function of r_H2
  ;                   [unitless]
  ;
  ; Output:
  ;   n_H2 - 2d array, describing the number density of emitting H2 gas in the disk, and depends on both
  ;          the radius of the H2 from the star and the vertical location in the disk of the H2 gas.
  ;          [cm-3]
  ;

  ; Conversions
  MSOLtoG = 1.99E+33      ; [M_solar --> g] conversion for M_H2
  AUtoCM = 1.496E+13      ; [AU --> cm] - conversion for r_H2 from cm to AU

  ; Define the constants and define g_l
  g_j = DBLARR( N_ELEMENTS(j_rot) )
  FOR j = 0, N_ELEMENTS(j_rot) - 1 DO BEGIN
  IF j_rot[j] MOD 2 EQ 0 THEN g_j[j] = (2*j_rot[j] + 1) ELSE g_j[j] = 3*(2*j_rot[j] + 1)    ; degeneracy of the rotation state of H2
  ENDFOR
  mu_H2 = 2.0*1.67E-24    ; mean molecular weight of H2, = 2* m_H [g]


  num_density_constant = (M_H2*MSOLtoG*(2 - density_grad))/(mu_H2*SQRT(8.0*(!PI)^3)*(r_c*AUtoCM)^2)

  ; Set up a FOR loop to take radius elements from the 1d arrays and multiply them by the vertical
  ; array (of height in the disk) in the density profile of the H2.
  ; n_H2 will also be a 2d array.
  n_H2 = DBLARR( N_ELEMENTS(rho_surface),N_ELEMENTS(rho_disk),N_ELEMENTS(j_rot) )

  FOR i = 0, N_ELEMENTS( rho_surface ) - 1 DO BEGIN
    FOR z = 0, N_ELEMENTS( rho_disk ) - 1 DO BEGIN
      ; Calculate n_H2
      IF partition_fnct[i] EQ 0.0 THEN BEGIN
        n_H2[i,z] = 0.0
      ENDIF ELSE BEGIN
        FOR k = 0, N_ELEMENTS(j_rot) - 1 DO BEGIN
        n_H2[i,z,k] = rho_disk[z]*(num_density_constant*rho_surface[i]*thermal_distr[i])/(partition_fnct[i]*Hz[i]*AUtoCM)*g_j[k]
        ENDFOR
      ENDELSE
    ENDFOR
  ENDFOR

;  stop

  ; Return the number density calculated
  RETURN, n_H2

END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - OPACITY OF H2 GAS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Once we have the number density of emitting H2 in the disk as a function of radius and the
; thermal cross section of said H2 gas, we can find an expression for the opacity of the gas
; in the disk as a function of radius and vertical height above the midplane of the disk.
;
FUNCTION opacity, n_H2, sigma
  ;
  ; Parameters Input:
  ;   n_H2 - 2d array, with rows = radius in the disk, columns = vertical height in the plane.
  ;          Represents the number density of H2 in the protoplanetary disk. [cm-3]
  ;   sigma - 1d array, the thermal cross section of H2 gas in the disk [cm2]?
  ;
  ; Output:
  ;   kappa_H2 - 2d array, opacity of the emitting H2 gas in the LTE disk. [cm+2 g-1]?
  ;

  ; Define the kappa 2d array
  kappa_H2 = DBLARR( N_ELEMENTS(n_H2[*,0,0]),N_ELEMENTS(n_H2[0,*,0]),N_ELEMENTS(n_H2[0,0,*]) )

  ; In a FOR loop, calculate kappa_H2
  FOR i = 0, N_ELEMENTS(n_H2[*,0]) - 1 DO BEGIN
    FOR z = 0, N_ELEMENTS(n_H2[0,*]) - 1 DO BEGIN
    ; Mulitply the values of sigma by the columns of n_H2
      kappa_H2[i,z,*] = sigma*n_H2[i,z]
    ENDFOR
  ENDFOR
;  stop

  RETURN, kappa_H2

END






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; PART 3 - RADIATIVE TRANSFER EQUATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;;
;;
;; Here, we apply all the above calculations to determine the intensity of the emission feature from the H2 in the disk.
;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - EFFECTIVE TAU OF H2 EMISSION PROGRESSION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
FUNCTION effective_tau, NH2, T_H2, column_H2, T_r, taueff
  ;
  ; We have read in the effective tau save files from Matt McJ, so
  ; now we have to match up column_H2 and T_r to calculate an effective
  ; optical depth.

  ; Check that the optical depth values of H2 make sense?
  ; Convert NH2 to alog10 units to compare to column_H2
  colH2 = alog10(NH2)

  ; Round colH2 and T_H2 by the appropriate step size of T_r and column_H2
  step_T = 100
  step_N = 0.2
  colH2 = ROUND(colH2/step_N)*step_N
;  TH2 = ROUND(T_H2/step_T)*step_T
  TH2 = T_H2

  tau_eff = DBLARR( N_ELEMENTS(colH2) )

  ; IF statements: Make sure NH2 and T_H2 are not out of bounds
  FOR k = 0, N_ELEMENTS(colH2[0,0,*]) - 1 DO BEGIN
  IF (TH2 LT min(T_r[*,k])) OR (TH2 GT max(T_r[*,k])) OR (colH2[0,0,k] LT min(column_H2[*,k])) OR (colH2[0,0,k] GT max(column_H2[*,k])) THEN BEGIN
    tau_eff[k] = 0.0
  ENDIF ELSE BEGIN
;    tau_index = WHERE( (colH2 EQ column_H2) AND (TH2 EQ T_r) )
    tau_index = WHERE( (colH2[k] LE column_H2[*,k]+0.1 AND colH2[k] GE column_H2[*,k]-0.1) AND (TH2 LE T_r[*,k]+50.0 AND TH2 GE T_r[*,k]-50.0) )
    IF tau_index NE -1 THEN BEGIN
      tau_eff[k] = taueff[tau_index]
    ENDIF ELSE tau_eff[k] = 0.0;print,'Something went wrong! :('

    ;    tau_eff = taueff[WHERE( (colH2 EQ column_H2) AND (TH2 EQ T_r) )]
  ENDELSE
  ENDFOR
;    stop

  RETURN, tau_eff

END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - Optical Depth Profile ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
FUNCTION optical_depth, n_H2, sigma, grid_dist, line_profile, disk_height, T_H2, column_H2, T_r, tau_eff ;grid_AU
  ;
  ; grid_dist = Hp (scale height) --> good
  ;

  ; Conversion for pc to cm
  PARSECtoAU = 206265.0    ; [pc --> AU]
  AUtoCM = 1.49E+13   ; [AU --> cm], have to convert the sightline to cm

  kappa_H2 = opacity(n_H2, sigma)

  ; Calculate tau
;  tau = DBLARR( N_ELEMENTS(line_profile[*,0]),N_ELEMENTS(line_profile[0,*]) )
  tau = DBLARR( N_ELEMENTS(kappa_H2[*,0,0]), N_ELEMENTS(n_H2[0,*,0]), N_ELEMENTS(n_H2[0,0,*]) )
  NH2 = DBLARR( N_ELEMENTS(kappa_H2[*,0,0]), N_ELEMENTS(n_H2[0,*,0]), N_ELEMENTS(n_H2[0,0,*]) )
  t_eff = DBLARR( N_ELEMENTS(kappa_H2[*,0,0]), N_ELEMENTS(n_H2[0,*,0]), N_ELEMENTS(n_H2[0,0,*]) )

;  stop
  FOR i = 0, N_ELEMENTS(kappa_H2[*,0,0]) - 1 DO BEGIN
    ; Keep a count of tau in z-direction for a given i
    tau_count = 0.0

    FOR z = 0, N_ELEMENTS(kappa_H2[0,*,0]) - 1 DO BEGIN
      tau_count += kappa_H2[i,z,*]*grid_dist[i]*AUtoCM*(FLOAT(z+1)/N_ELEMENTS(kappa_H2[0,*,0]));*line_profile[i];*AUtoCM;*PARSECtoAU;grid_AU*AUtoCM
      tau[i,z,*] = tau_count

      ; Determine tau_eff from Matt's table, using N(H2) and T(r)
      ; Calculate N(H2)
      NH2[i,z,*] = (FLOAT(z+1)/N_ELEMENTS(kappa_H2[0,*,0]))*n_H2[i,z,*]*grid_dist[i]*AUtoCM   ; units: cm-2

      ; Check that we match up the column density, temperature info of the radial position in the disk
      ; and pull from that the associated tau_eff value
      t_eff[i,z,*] = effective_tau(NH2[i,z,*],T_H2[i],column_H2,T_r,tau_eff)

      FOR k = 0, N_ELEMENTS(t_eff[0,0,*]) - 1 DO BEGIN
      IF t_eff[i,z,k] EQ 0.0 THEN BEGIN
        tau[i,z,k] = tau_count[k]
      ENDIF ELSE BEGIN
        tau[i,z,k] = tau_count[k]*t_eff[i,z,k]
      ENDELSE
      ENDFOR
    ENDFOR
  ENDFOR
;  stop

  RETURN, tau

END


;;
;
;;;;;;;;;;;;;;;;;;;;;; FUNCTION - CORRECT LyA FLUX ;;;;;;;;;;;;;;;;;;;;;;;;;
;
; To accurately represent the flux the H2 absorbs to the excited state and emits back down,
; we have to convert the given from Kevin's target flux list of LyA values observed from
; Earth, with the associated distance to the target and A_v value between the target
; and the observer.
;
; Once the flux is corrected to some nominal position in the disk (say, from some distance from
; the star, to represent LyA photons NOT from the stellar surface but from the accretion shock
; front that is exciting a large flux of LyA photons), we can assume that this is the flux the
; H2 molecules see in the disk and
;
FUNCTION unred_flux_to_H2, start_flux, wave, A_v, r_H2, dist, r_LyA
  ;
  ; INPUT:
  ;   start_flux - the flux at the given LyA wavelength from the target observed at Earth
  ;   LyAwave - the pumping LyA wavelength for the H2 excitation [Ang]
  ;   A_v - V-mag extinction of the target from the observer
  ;   r_H2 - 2D matrix, all the radial values of where the H2 is located in our model disk
  ;   dist - distance from the target to the observer [pc]
  ;   r_LyA - location from the target star that the LyA originates (assumed constant)
  ;
  ; OUTPUT:
  ;   flux_at_H2 - 2d array, unreddened flux the H2 sees in the disk
  ;

  ; Define relevent constants
  PARSECtoAU = 206265.0    ; [pc --> AU]
;  AUtoCM = 1.49E+13        ; [AU --> cm], have to convert the sightline to cm
  ANGtoCM = 1.0E-8         ; [Ang --> cm]

  ; Define Rv and therefore Ebv (color excess)
  Rv = 3.1    ; typical ISM value
  Ebv = A_v / Rv

  ; Use ccm_unred.pro to calculate the unreddened LyA flux as observed from Earth
  ; Calls:
  ;   1. wavelegnth - LyA pumping wavelength of interest (Angstroms)
  ;   2. start_flux - associated flux with the LyA pumping wavelength
  ;   3. E(B-V) - color excess, found by relation R_v = A_v/E(B-V)
  ;   5. R_v (optional) - constant to describe ISM distribution of dust
  ; Return:
  ;   4. f_unred - Corrected flux of LyA pumping
  ;
  flux_unreddened = DBLARR( N_ELEMENTS(wave) )
  FOR wa = 0, N_ELEMENTS(wave) - 1 DO BEGIN
    ccm_unred, wave[wa], start_flux[wa], Ebv, funred ;, R_V = Rv
    flux_unreddened[wa] = funred
  ENDFOR
  ; Calculate the LyA flux at the target
  flux_LyA_target = flux_unreddened * ( (dist*PARSECtoAU)^2 / (r_LyA)^2 )  ; both distances are in AU

;  flux_at_H2 = DBLARR( N_ELEMENTS(r_H2) )
  flux_at_H2 = DBLARR( N_ELEMENTS(r_H2[*,0]),N_ELEMENTS(wave) )

  ; Calculate flux at H2, given a radius where the H2 is located
  FOR i = 0, N_ELEMENTS( r_H2 ) - 1 DO BEGIN
;    FOR j = 0, N_ELEMENTS( r_H2[0,*] ) - 1 DO BEGIN
      flux_at_H2[i,*] = flux_LyA_target * ( (r_LyA)^2 / (r_H2[i])^2 )                  ; both distances are in AU
;    ENDFOR
  ENDFOR
;  stop
  RETURN, flux_at_H2    ; units = [ergs cm-2 s-1 Ang-1]

END







;;
;
;;;;;;;;;;;;;;;;;;;;;;;; FUNCTION - TOTAL INTENSITY ;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
FUNCTION total_intensity, tau, sight_line, LyA_flux, T_H2, v_obs, v_binned, H2_grid_size, H2_radius, H_z, distance, A_v, starmass, inclination, turbulence, characteristic_radius, totalH2mass, temp_grad, dens_grad, height
  ;
  ; INPUTS:
  ;
  ; OUTPUT:
  ;   total_intensity - 1d array, of same size as binned v_obs, bins and sums up the
  ;       total intensity from the emitting H2 gas, to give us our final emission line
  ;       feature!!
  ;

  ; Define Rv and therefore Ebv (color excess)
  Rv = 3.1    ; typical ISM value
  Ebv = -1*A_v / Rv     ; -1 indicated reddening the flux

  ; Define grid
  H2_grid = grid_radial_gas(H2_radius)

  ; Set a fraction of H2 that emits from the progression
  ; From Herczeg et al 2004, eta = 0.25
  eta = 0.25

  ; Create a 2d and binned intensity array, to represent intensities of H2 emission at a given v_obs
;  intens_H2 = DBLARR( N_ELEMENTS(kappa[*,0]),N_ELEMENTS(kappa[0,*]) )
  intens_H2 = DBLARR( N_ELEMENTS(v_obs[*,0]),N_ELEMENTS(v_obs[0,*]),N_ELEMENTS(tau[0,0,*]) )

  ; Determine conversion from pc to cm
  PARSECtoAU = 206265.0    ; [pc --> AU]
  AUtoCM = 1.49E+13   ; [AU --> cm], have to convert the sightline to cm

  ; Calculate the intensity of the H2 emission at each grid point
  FOR i = 0, N_ELEMENTS(tau[*,0,0]) - 1 DO BEGIN
    ; Start a counter for the flux estimate from the H2 emission
    flux_H2_disk = DBLARR(N_ELEMENTS(tau[0,0,*]));*0.0

    FOR z = 0, N_ELEMENTS(tau[0,*,0]) - 1 DO BEGIN
      flux_emission = LyA_flux[i,*]*(1.0 - EXP(-1.*tau[i,z,*]))*eta;*tau[i,z]
      flux_H2_disk += flux_emission
    ENDFOR

    ; For r,phi locations (i,j)
    FOR j = 0, N_ELEMENTS(sight_line[0,*]) - 1 DO BEGIN
      ; Calculate radiative transfer in the disk at absorbing H2
;      flux_H2_disk = kappa[i,j]*LyA_flux[i,j]*EXP(-1.*tau[i,j])*H2_grid[i,j]*AUtoCM


      ; Calculate F_H2_obs, the observer frame flux from the H2
;      flux_H2_obs = flux_H2_disk*( (H2_grid[i,j])^2 / (sight_line[i,j]*PARSECtoAU)^2 )
      flux_H2_obs = flux_H2_disk*( H2_grid[i]*H2_grid[i] / (sight_line[i,j]*PARSECtoAU)^2 ) * COS(inclination*!pi/180)*COS(inclination*!pi/180)

      ; Redden the flux as it travels through the ISM for the given A_v value
 ;     ccm_unred, H2_wave, flux_H2_obs, Ebv, flux_H2_red    ; Not working right now but don't know why...
      intens_H2[i,j,*] = flux_H2_obs    ; flux_H2_red
    ENDFOR
  ENDFOR


;  stop
;  plot,v_obs,intens_H2,xr=[-50,50],psym=3
;  stop

  ;;;;; OPTIONAL - INNER v. OUTER RADIUS H2 CONTRIBUTIONS ;;;;;
  ; Make a dummy variable to create a plot of inner radius and outer radius emission
  ; line contributions
;  dummy_var = inner_outer_radii_flux_distribution(H2_radius, intens_H2, v_obs, v_binned, LyA_flux, sight_line, H2_grid_size, T_H2, starmass, inclination, turbulence, characteristic_radius, totalH2mass, temp_grad, dens_grad, height)

;  print,SYSTIME()
;  ; Collapse down the 2d intensity to a 1d binned intensity by v_obs points
;  bin_factor = 1.0    ; km s-1
;  intensity_binned = define_collapsed_line_profile(v_binned, bin_factor, intens_H2, v_obs)

;  plot,v_binned,intensity_binned/N_ELEMENTS(intens_H2[0,*]),xr=[-150,150]
;  stop

  intens_H2 = intens_H2*2.
  RETURN, intens_H2;intensity_binned

END




;;;;;;;;;;;;;;;;;; NEW FUNCTION: v_bin by emission wavelength ;;;;;;;;;;;;;;;;;;;
;
; The previous function creates a 3D array that describes the radial,phi locations of
; resulting H2 flux for a given progression considered for LyA-pumped H2 to some excited
; state.
; Now, we need to bin the intensity array into a 2D array of v_obs by emission transition
; and include the Branching Ration probabilities of the transition occuring.
;
FUNCTION all_emission_line_intensities, A_v, vobs, v_binned, intens_h2, emiss_lines, b_mn, j_emission, j_prog
  ;
  ; Inputs:
  ;   A_v - double, extinction coefficient for reddening
  ;   vobs - 2D array, observed velocity at each radial,phi grid point in disk
  ;   v_binned - 1D array, holds the v_obs values of interest to plot again intensity observed
  ;   intens_h2 - 3D array (radius,phi,progression), intensity array by disk location and
  ;               H2 absoprtion progression
  ;   emiss_lines - 1D array, all the emission lines of interest
  ;   b_mn - 1D array, branching ratios for all transitions (emission lines) from a given progression
  ;   j_emission - 1D array, all j locations of the starting progression by transition
  ;   j_prog - 1D array, starting j values of each given progression. Use to compare against
  ;            j_emission to determine where each transition originates
  ;
  ; Output:
  ;   intensity - 2D array (vobs,emission line), all the final (and reddened) flux
  ;               outputs for each H2 emission line
  ;

  ; Define Rv and therefore Ebv (color excess)
  Rv = 3.1    ; typical ISM value
  Ebv = -1*A_v / Rv     ; -1 indicated reddening the flux

  ; Define the intensity array
  intensity = DBLARR( N_ELEMENTS(v_binned),N_ELEMENTS(emiss_lines) )

  ; Define the intensity binned for each progression
  print,SYSTIME()
  ; Collapse down the 2d intensity to a 1d binned intensity by v_obs points
  bin_factor = 1.0    ; km s-1
  intensity_binned = define_collapsed_line_profile(v_binned, bin_factor, intens_H2, vobs)

  ; Make a FOR LOOP to:
  ; 1. Locate the correct progression by comparing j_emission to j_prog
  ; 2. Multiply intensity by the correct branching ratio (dependent on emiss_line)
  ; 3. Redden by the correct wavelength (emiss_line)
  FOR i = 0, N_ELEMENTS(emiss_lines) - 1 DO BEGIN
    ; 1. Find correct index for the progression emission line arises from
    prog_index = WHERE(j_emission[i] EQ j_prog)

    ; 2. Append intensity_binned of correct progression intensity into intensity array,
    ;    and mulitply by branching ratio
    intens_unred = intensity_binned[*,prog_index]*b_mn[i]

    ; 3. Redden by using ccm_unred
    ccm_unred, emiss_lines[i], intens_unred, Ebv, intens_red
    intensity[*,i] = intens_red

  ENDFOR
;  stop

  RETURN, intensity


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
;  stop
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
;  stop


  RETURN, total_line_profile

END


;;
;
;;;;;;;;;;;; OPTIONAL FUNCTION - INNER v. OUTER RADIUS from H2 FLUX ;;;;;;;;;;;;;;;
;
;
FUNCTION inner_outer_radii_flux_distribution, rH2, intensity, v_obs, v_binned, LyAflux, sight_line, H2_grid_size, T_rH2, starmass, inclination, turbulence, characteristic_radius, totalH2mass, temp_grad, dens_grad, height
  ;
  ; INPUT:
  ;   rH2 - 2d array, radial location of absorbin H2
  ;   intensity - 2d array, emitting intensity of the H2 at a radial location in the disk
  ;   v_obs - 2d array, observed radial velocity of H2
  ;   v_obs_binned - 1d, binned v_obs values
  ;   LyAflux - 2d array, flux of LyA pumping wavelength as seen by the H2 at a given radius
  ;   sight_line - 2d array, sightline distance to each grid point on the disk
  ;   H2_grid_size - double, grid size used, in AU
  ;
  ; OUTPUT:
  ;   Plots the inner and outer radial contributions to the H2 emission line
  ;   dummy_var - some placeholder to call this function in the main program
  ;


  ; Determine conversion from pc to cm
  PARSECtoAU = 206265.0    ; [pc --> AU]
  AUtoCM = 1.49E+13   ; [AU --> cm], have to convert the sightline to cm

  ; Make arrays for the inner, outer intensity distribution values, and the
  ; corresponding velocities
  inner_vel_H2 = DBLARR( N_ELEMENTS(v_obs[*,0]), N_ELEMENTS(v_obs[0,*]) )     ; r < 1.0 AU
  inner_intens_H2 = DBLARR( N_ELEMENTS(intensity[*,0]), N_ELEMENTS(intensity[0,*]) )     ; r < 1.0 AU
  outer_vel_H2 = DBLARR( N_ELEMENTS(v_obs[*,0]), N_ELEMENTS(v_obs[0,*]) )     ; r > 1.0 AU
  outer_intens_H2 = DBLARR( N_ELEMENTS(intensity[*,0]), N_ELEMENTS(intensity[0,*]) )

  ; Set boundary between inner and outer disk
  inner_outer_boundary = 1.0    ; AU

  ; Loop through the radius array, finding where certain radial points sampled
  ; in the disk fall
  FOR i = 0, N_ELEMENTS(v_obs[*,0]) - 1 DO BEGIN
    FOR j = 0, N_ELEMENTS(v_obs[0,*]) - 1 DO BEGIN
      ; Test to see where rH2[i,j] will contribute
;      IF rH2[i,j] LE inner_outer_boundary THEN BEGIN
      IF rH2[i] LE inner_outer_boundary THEN BEGIN
        inner_vel_H2[i,j] = v_obs[i,j]
        inner_intens_H2[i,j] = intensity[i,j]
      ENDIF ELSE BEGIN
        outer_vel_H2[i,j] = v_obs[i,j]
        outer_intens_H2[i,j] = intensity[i,j]
      ENDELSE
    ENDFOR
  ENDFOR

  ; Now, bin the intensities like we did in the Function above.
  bin_factor = 1.0    ; km s-1
  intensity_bin_inner = define_collapsed_line_profile(v_binned, bin_factor, inner_intens_H2, inner_vel_H2)
  intensity_bin_outer = define_collapsed_line_profile(v_binned, bin_factor, outer_intens_H2, outer_vel_H2)

  ; Plot both to see how each contributes to the total intensity
  ang = string("305b)
;  set_plot, 'ps'
;  device, /color, /times, /isolatin1, /landscape,$
;        filename='C:\Users\keho8439\Documents\Comps2\model_comparisons\model_height_6.ps'
;  plot,v_binned,intensity_bin_outer,linestyle=0,thick=2.5,xr=[-100,100],color=!black,xtitle='observed velocity [km s-1]',ytitle='flux [ergs cm-2 s-1 '+ang+'-1]'
;  oplot,v_binned,intensity_bin_inner,linestyle=0,thick=2.5,color=!blue
;  legend_astron,['Outer radius, r > 1 AU','Inner radius, r < 1 AU'], linestyle=[0,0],thick=[2.5,2.5],color=[!black,!blue]
;  legend_astron,['Stellar mass = '+STRTRIM(starmass,2)+' M_sol',$
;                 'Disk inclination = '+STRTRIM(inclination,2)+' degrees',$
;                 'Turbulent velocity = '+STRTRIM(turbulence,2)+' km s-1',$
;                 'Temperature of gas = '+STRTRIM(T_rH2,2)+' K',$
;                 'Temperature gradient of gas = '+STRTRIM(temp_grad,2),$
;                 'Char. radius = '+STRTRIM(characteristic_radius,2)+' AU',$
;                 'Density gradient of gas = '+STRTRIM(dens_grad,2),$
;                 'Mass of emitting H2 = '+STRTRIM(totalH2mass,2)+' M_sol',$
;                 'Gas height = '+STRTRIM(height,2)+' Hz'],/right_legend
;  device, /CLOSE_FILE
;
;  set_plot,'WIN'

  plot,v_binned,intensity_bin_outer,linestyle=0,thick=1.5,xr=[-150,150],xtitle='observed velocity [km s-1]',ytitle='flux [ergs cm-2 s-1 '+ang+'-1]'
  oplot,v_binned,intensity_bin_inner,linestyle=0,thick=1.5,color=!purple
  legend_astron,['Outer radius, r > 1 AU','Inner radius, r < 1 AU'], linestyle=[0,0],thick=[1.5,1.5],color=[!white,!purple]

  stop

  RETURN, intensity_bin_inner

END








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; START MAIN PROGRAM - CREATING THE MODELS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro disk_profile_v2, args=args
;
; args - points to the directory you are storing the created SAV files. Will take the numbers (indices/parm values) in the last file created
;        to continue the model creation where the program left off, should the program be pre-maturely terminated
;

; Set plot commands early
device,decompose=0
setplotcolors

;;; Constrain our variables:
;;    v_obs    (-200,200) km s-1 <-- USED FOR ALL CALCULATIONS, taken from .sav files (cosvelocity.pro)
;     labwave  (list of specific rest wavelengths picked out for each star) Ang, also taken from .sav files (same as above)
;;    r_H2     (0.03,14.0) AU
;;    M_H2     (0.08, 2.0) * 10^-2   M_star
;;    T_1AU    (1500, 3000) K
;;    v_turb   (0,10) km s-1
;;    q        (-3.0,3.0)
;;    i        (0,90) degrees
;;    M_star   (0.08, 2.0) M_star
;;   [phi      (0; -30,30) degrees ]
;;
;; Variables to calculate:
;;    r_obs = (v_obs/cos[phi]*sin[i])^-2 * G * M_star
;;    T_H2  = T_1AU * (r_obs/1AU)^-q
;;    line_profile
;;
;
;
; Define our reference path to pull files and to store model info
; Restore the save file with the velocity information from a target.
ref_path_obs = 'C:\Users\keho8439\Documents\Comps2\obs\'
ref_path_models = 'C:\Users\keho8439\Documents\Comps2\models\'


;;;;;;;;;;;;;;;;;;;;;;;;; PARAMETER SPACE VARIABLES ;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  Define the variables to be used
;

Hp_height = [4]      ; Outer-most scale height the emitting disk is located. units: (z/r)
gamma = [0.75]        ; density gradient (gamma)
T_1AU = [2000.0]     ; Temperature of the gas at 1AU. units: K
q = [0.25]            ; Temperature gradient (for a power law). units: unitless
r_char = [3.0]       ; Characteristic radius of the H2 emission. units: AU
M_H2 =[5.0E-12]      ; Total mass of H2 contributing to the emission feature. units: M_solar



stop

; Specify the sample size of the grid
AUsample = 0.005      ; sampling in the radial component of the disk desired
AUsize_inner = 5000
AUsize_outer = 100
phi_density = 180.


doppler = 0.0  ; Doppler velocity of the target. units: km s-1
turbulence = 0.5     ; Turbulent velocity of the gas in the disk. units: km s-1

;
; We need to input the Einstein A coefficients by hand for a given transition line for the H2 emission.
; Below is the list of wavelengths that the transition line occurs.
; Matching the correct correct wavelength index location to the correct EinsteinA value, input the proper
; Einstein A coefficient by looking them up in Abgrall et al 1993 (either Lyman or Werner bands).
;
; After we define which wavelength we will work in, we can see which reference wavelength matches up to the wavelength being
; used for the target and then match the correct array index to that of EinsteinA to use that correct EinsteinA later on.
;

;;;;; State all states
J_initial = [0, 1, 5, 6]
J_progression = [1, 2, 4, 7]
lower_energy = [1.00, 1.02, 1.20, 1.27]*1.60217657E-12     ; energy of the lower state of the transition, in ergs
LyA_pumping = [1217.21, 1217.64, 1216.07, 1215.73]
Aul_total = [1.8,1.8,1.7,1.7]*1.0E9                     ; total Aul for each progression

;;;;; State all transition important quanities
lambda_ref = [1338.56,1342.26,1393.96,1398.95,1402.65,1446.12,1460.17,1463.83,1467.08,1489.57,1500.45,1504.76,1521.59,1524.65,1580.67]      ;[Ang]
EinsteinA = [3.1E+8,2.8E+8,1.6E+8,2.6E+8,2.3E+8,1.4E+8,1.5E+8,1.4E+8,1.3E+8,1.6E+8,1.7E+8,2.0E+8,6.0E+7,1.9E+8,1.1E+8]   ;[s-1]
v1 = 2   ; ground vibration state, before LyA pumping
J1 = [0,1,1,0,1,5,0,1,6,5,6,5,0,6,6]   ; ground rotation state, before LyA pumping
;LyA_pumping = [1217.21,1217.64,1217.64,1217.21,1217.64,1216.07,1217.21,1217.64,1215.73,1216.07,1215.73,1216.07,1217.21,1215.73,1215.73]
;lower_energy = [1.00,1.02,1.02,1.00,1.02,1.20,1.00,1.02,1.27,1.20,1.27,1.20,1.00,1.27,1.27]*1.60217657E-12     ; energy of the lower state of the transition, in ergs

J2 = [1,2,2,1,2,4,1,2,7,4,7,4,1,7,7]   ; excited rotation state, right after LyA pumping
J3 = [2,3,1,2,3,5,2,3,8,3,6,5,2,8,8]   ; ground state after emission from excited H2 molecule

Aul_temp = DBLARR(N_ELEMENTS(J2))
for i = 0, N_ELEMENTS(J2)-1 DO Aul_temp[i] = Aul_total[WHERE(J2[i] eq J_progression)]
B_ul = EinsteinA/Aul_temp

; Read in fluormod_trans file for total H2 absorption coefficient values
;idl_dir = 'C:\Users\keho8439\Documents\IDL code\kf_cos_fit\'
;fluormod = 'fluormod_trans.idlsav'
;
;restore,idl_dir+fluormod
;
;atot = h2_trans.atotal      ; Summed Aul values by progression
;v_up = h2_trans.vu          ; Excited state v
;j_up = h2_trans.ju          ; Excited stat j
;
;; Find where v = 0, J = [1,2], and v = 1, J = [4,7] to find correct atot
;vj_01 = where(v_up eq 0 and j_up eq 1)
;vj_01 = vj_01[0]
;vj_02 = where(v_up eq 0 and j_up eq 2)
;vj_02 = vj_02[0]
;vj_14 = where(v_up eq 1 and j_up eq 4)
;vj_14 = vj_14[0]
;vj_17 = where(v_up eq 1 and j_up eq 7)
;vj_17 = vj_17[0]
;
;Aul_total = [atot[vj_01],atot[vj_02],atot[vj_14],atot[vj_17]]
;Aul_temp = DBLARR(N_ELEMENTS(J2))
;for i = 0, N_ELEMENTS(J2)-1 DO Aul_temp[i] = Aul_total[WHERE(J2[i] eq J_progression)]
;B_ul = EinsteinA/Aul_temp

;stop

;stop


; Define the angular momentum array (or integer if it's the same for all the lines).
; I need to check this; I don't know how l (ang. mom.) is defined in the H2, so once I find out
; I can include this. For now, it will act as a dummy
;


;;;;;; FULL DIRECTORIES + TARGET LISTS ;;;;;;;;;;;;
; For the models, we will make them specific to each target because
filedir = ['AATau\','BPTau\','CSCha\','DFTau\','DMTau\','GMAur\','LKCA15\','UXTau\','V4046SGR\',$
  'HNTau\','RECX11\','RECX15\','SUAur\','TWHya\'] ; linux
filename = ['AATau_011911kf.sav','BPTau_091111.sav','CSCHA_060711.sav','dftau_reproc052410.sav',$
  'DMTAU_codd.sav','GMAUR_022011kf.sav','LKCA15_coaddJUL10.sav','UXTAU_010511kf.sav',$
  'V4046SGR_071910kf.sav','HNTAU_coaddJUL10.sav','RECX11_021811kf.sav',$
  'RECX15_coaddJUL10.sav','SUAur_040111kf.sav','TWHya_STIS_E140M.sav'];'AATAU_COSall2013.sav'
fileLyAflux = ['AATAU_allFUV_abs_f_053012.sav','BPTAU_allFUV_abs_f_053012.sav',$
  'BPTAU_allFUV_abs_f_053012.sav','DFTAU_allFUV_abs_f_053012.sav',$
  'DMTAU_allFUV_abs_f_053012.sav','GMAUR_allFUV_abs_f_053012.sav',$
  'LKCA15_allFUV_abs_f_053012.sav','UXTAU_allFUV_abs_f_053012.sav',$
  'V4046SGR_allFUV_abs_f_053012.sav','HNTAU_lyman_profiles.sav',$
  'RECX11_lyman_profiles.sav','RECX15_lyman_profiles.sav','SUAUR_lyman_profiles.sav',$
  'TWHYA_lyman_profiles.sav']
targ_distance = [140.,140.,160.,140.,140.,140.,140.,140.,140.,83.,140.,97.,97.,140.,54.]
targ_Av = [0.5,0.5,0.8,0.6,0.0,0.1,0.6,0.2,0.0,0.5,0.0,0.0,0.9,0.0]
targ_r_in = [0.1,0.01,3.0,0.01,0.2,0.5,0.1,0.5,0.2,0.16,0.29,0.21,0.92,0.1]   ; all are under actual values to compensate for anything missed in the inner disk region
targ_r_mid = [1.0,1.0,10.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
targ_r_out = [50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0]
targ_inc = [75.,30.,60.,85.,35.,55.,49.,35.,36.,60.0,70.,60.,62.,4.]
targ_mass = [0.8,0.73,1.05,0.19,0.5,1.2,0.85,1.3,0.85+0.69,0.85,0.80,0.40,2.30,0.60]
targ = 4    ; DM Tau

;filedir = ['AATau\','BPTau\','CSCha\','DFTau\','DMTau\','GMAur\','LKCA15\','UXTau\','V4046SGR\']
;filename = ['AATau_011911kf.sav','BPTau_091111.sav','CSCHA_060711.sav','dftau_reproc052410.sav','DMTAU_codd.sav','GMAUR_022011kf.sav','LKCA15_coaddJUL10.sav','UXTAU_010511kf.sav','V4046SGR_071910kf.sav']
;fileLyAflux = ['AATau_allFUV_abs_f_053012.sav','BPTau_allFUV_abs_f_053012.sav','BPTau_allFUV_abs_f_053012.sav','DFTAU_allFUV_abs_f_053012.sav','DMTAU_allFUV_abs_f_053012.sav','GMAUR_allFUV_abs_f_053012.sav','LKCA15_allFUV_abs_f_053012.sav','UXTAU_allFUV_abs_f_053012.sav','V4046SGR_allFUV_abs_f_053012.sav']
;targ_distance = [140.,140.,160.,140.,140.,140.,140.,140.,140.,83.]
;targ_Av = [0.5,0.5,0.8,0.6,0.0,0.1,0.6,0.2,0.0]
;targ_r_in = [0.1,0.01,3.0,0.01,0.2,0.5,0.1,0.5,0.2]   ; all are under actual values to compensate for anything missed in the inner disk region
;targ_r_mid = [1.0,1.0,10.0,1.0,1.0,1.0,1.0,1.0,1.0]
;targ_r_out = [50., 50., 50., 50., 50., 50., 50., 50., 50.]
;targ_inc = [75., 30., 60., 85., 35., 55., 49., 35., 36.]
;targ_mass = [0.8, 0.73, 1.05, 0.19, 0.5, 1.2, 0.85, 1.3, 0.85+0.69]

;; TEST TARGET DIRECTORY + LIST
;filedir = ['AATau\']
;filename = ['AATau_011911kf.sav']
;fileLyAflux = ['AATau_allFUV_abs_f_053012.sav']
;targ_distance = [140.]
;targ_Av = [0.5]
;targ_r_in = [0.1]
;targ_r_mid = [1.0]
;targ_r_out = [100.0]
;targ_inc = [75.]      ; Disk inclination, defined per target from sub-mm lit.
;targ_mass = [0.8]     ; Stellar mass of the target


; Read in partition function info
partition_save = 'disk_partition_values.sav'
restore,ref_path_models+partition_save      ; Arrays: J, T_H2_range, Z_TH2

; Read in the effective tau tables, given to me by Matt McJ.
taueff_saves = ref_path_models+'tau_eff_'+['1217.21','1217.64','1216.07','1215.73']+'.sav'



;;;;;;;;;;;;;;;;; PRELIMINARY LOOP: TARGET INFO ;;;;;;;;;;;;;;;;;;;;;;
;; This pulls the following infomation for the models:
;; Target name --> opens velocity, wavelength SAVE files
;;    v_obs (2D arr)
;;    labwave (1D arr)
;;    v_exc (1D arr)
;;    v_ground (1D arr)
;;    j_exc (1D arr)
;;    j_ground (1D arr)


FOR target = targ, targ DO BEGIN
;FOR target = 0, N_ELEMENTS(filename) - 1 DO BEGIN
; Open the files to snag the target name, LyA flux values, etc.
restore,ref_path_obs+filedir[target]+filename[target]
restore,ref_path_obs+filedir[target]+fileLyAflux[target]

starname = targname

; Print where we are in the program
print,'Target: '+starname

; Make a directory to store all the SAVE files
FILE_MKDIR,ref_path_models+filedir[target]+'disk_profiles_v2'
print,'directory '+ref_path_models+filedir[target]+'disk_profiles_v2'+' successfully created!'

; Save the parameters in a save file, to open later to compare parameter space with the data.
parm_file = ref_path_models+filedir[target]+'disk_profiles_v2\'+'parameter_list.sav'
save, T_1AU, q, gamma, r_char, M_H2, Hp_height, file=parm_file

;stop

distance = targ_distance[target]
inclination = targ_inc[target]
starmass = targ_mass[target]
Av = targ_Av[target]
r_in = targ_r_in[target]
r_mid = targ_r_mid[target]
r_out = targ_r_out[target]

r_LyA = 0.01      ; AU, a constant for now, drives where LyA photons originate in the disk/accretion front


; Define the disk in (r,phi)
;grid_from_radius, inclination, distance, r_in, r_mid, r_out, AUsize_inner, AUsize_outer, phi_density, xgrid, ygrid, r_grid, phi_grid, r_H2, phi_H2
;r_grid = r_grid*distance
;phi_H2 = phi_H2*180/!pi
;phi_grid = phi_grid*180/!pi

; Restore the wavelength and velocity SAVE files to use in each model for that target.
restore,ref_path_obs+filedir[target]+starname+'_H2wavelengths.sav'   ; labwave, v_exc, v_ground, j_exc, j_ground

print,SYSTIME()


;;;;;;;;;;;;;;;;;;;; FOR-LOOP TIIIIIIIME! ;;;;;;;;;;;;;;;;;;;;;;;
; Loop #1: lambda_0 --> so we can get v_obs once
; THIS IS NOT CORRECT, MAKE v_obs OUTSIDE OF LOOP!

; SO, since we are now calculating all transitions simultaneously in v2.0,
; get rid of restwave dependence of the models, and instead make a directory
; for the target.

;FOR i = 0, 0 DO BEGIN
  ; When we get to a new wavelength, make a new directory for that
  ; wavelength in the model directory
;  wavename = STRTRIM(labwave[i],2)
;  FILE_MKDIR,ref_path_models+filedir[target]+wavename
;  print,'directory '+wavename+' successfully created!'


;  print,'Beginning new wavelength: '+STRTRIM(labwave[i],2)
;  print,'Iterating loop '+STRTRIM(i+1,2)+' of '+STRTRIM(N_ELEMENTS(labwave),2)

  ; Define this currect labwave as a variable name
;  wave_index = WHERE(lambda_ref EQ labwave[i])    ; matching up our list to the testing wavelength of the target
;  restwave = lambda_ref[wave_index]
;  absorption_energy = lower_energy[wave_index]
;  transition_wave = LyA_pumping[wave_index]
  ; Make a 2D array to hold the column, T_r, and tau_eff values for all progressions
  restore,taueff_saves[0]      ; Gives back: column_H2, T_r, tau_eff
  column_tau = DBLARR( N_ELEMENTS(column_H2),N_ELEMENTS(LyA_pumping) )
  T_tau = DBLARR( N_ELEMENTS(column_H2),N_ELEMENTS(LyA_pumping) )
  eff_tau = DBLARR( N_ELEMENTS(column_H2),N_ELEMENTS(LyA_pumping) )

  FOR pr = 0, N_ELEMENTS(LyA_pumping) - 1 DO BEGIN
    restore,taueff_saves[pr]
    column_tau[*,pr] = column_H2
    T_tau[*,pr] = T_r
    eff_tau[*,pr] = tau_eff
  ENDFOR

;  stop
  ; Locate where the appropriate start flux is located from the flux file
  ; wave = LyAwave
  ; flux = LyAflux_red
;  flux_index = 0
;  FOR index = 0, N_ELEMENTS( LyAwave ) - 1 DO BEGIN
;    IF LyAwave[index] LE transition_wave THEN BEGIN
;      flux_index = index
;    ENDIF ELSE BREAK
;  ENDFOR
;  flux_index = WHERE( LyAwave EQ transition_wave )

  start_flux_red = DBLARR( N_ELEMENTS(LyA_pumping) )
  FOR fl = 0, N_ELEMENTS(LyA_pumping) - 1 DO BEGIN
    flux_index = WHERE( LyA_pumping[fl] LE LyAwave+0.01 AND LyA_pumping[fl] GE LyAwave-0.01)
    flux_index = flux_index[0]
    start_flux_red[fl] = LyAflux_red[flux_index]
  ENDFOR
;  stop

;  print,'Verifying wavelength = '+STRTRIM(restwave,2)+' Ang'
;  print,'Ground state energy of the H2 absorption species: '+STRTRIM(absorption_energy,2)+' ergs'
;  print,'LyA pumping wavelength = '+STRTRIM(transition_wave,2)+' Ang'
;
;
  ;;;;;;;;;;;;;;;;;;;;;;; OPTIONAL - Indexing current file list of parm space explored to continue ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Check to see if a directory argument was added to the program.
  if (not keyword_set(args)) then begin
    ; Default the indices to 0
    i_t1au = 0
    i_q = 0
    i_gamma = 0
    i_hp = 0
;    i_dop = 0
;    i_turb = 0
    i_rch = 0
    i_mh2 = 0

  endif else begin
    ; Search through all the files, making an array of all the names
    list_files = FILE_SEARCH(args+'*')

    ; Open the very last file
;    restore,list_files[N_ELEMENTS(list_files)-1]
    filelist = list_files[N_ELEMENTS(list_files)-1]

    ; Use STRMID to determine index value
    i_hp = FIX( STRMID(filelist,65,1) )
    i_gamma = FIX( STRMID(filelist,67,1) )
    i_t1au = FIX( STRMID(filelist,69,1) )
    i_q = FIX( STRMID(filelist,71,1) )
;    i_dop = FIX( STRMID(filelist,73,1) )
;    i_turb = FIX( STRMID(filelist,75,1) )
    i_rch = FIX( STRMID(filelist,77,1) )
    i_mh2 = FIX( STRMID(filelist,79,1) )
  endelse


;  stop





  ; No, don't do this! Above, define a generic v_obs over [-500,500] km s-1, sampling 1 km s-1
  ; Grab v_obs from the velocity file above
;;  v_obs = velocity[*,i]           ; velocity = 2d array from H2 velocity SAV file above

  ; Make a v_binned array for all the v_obs values we expect to measure in the
  ; line of sight.
  v_obs_bin = [ FINDGEN(300) - 300.0, FINDGEN(300) + 1.0 ]


  ;;;;;;;;;;;;;;;;;;;;;;;;;; DEFINE Intrinsic Cross Section of H2 (FUNCTION) ;;;;;;;;;;;;;;;;;;;;;;;
  ; For later, define sig0 (the initial cross section of the H2 at rest) because it only depends
  ; on the rest wavelength.
  sig0 = sigma_0(LyA_pumping, J_initial, J_progression, Aul_total)
  print, sig0
;  stop

  ; Loop #2: disk height
  FOR j = i_hp, N_ELEMENTS(Hp_height) - 1 DO BEGIN
    ; Define as a variable
    disk_height = Hp_height[j]
;    stop


    ; Loop #3: i_disk
;    FOR k = 0, N_ELEMENTS(i_disk) - 1 DO BEGIN
;      ; Define as a variable
;      inclination = i_disk[k]

      ; Loop #4: density gradient, gamma
      FOR l = i_gamma, N_ELEMENTS(gamma) - 1 DO BEGIN
        ; Define as a variable
        density_gamma = gamma[l]


        ; Loops #5: T_1AU
        FOR m = i_t1au, N_ELEMENTS(T_1AU) - 1 DO BEGIN
          ; Define as a variable
          temperature = T_1AU[m]


          ; Loop #6: q
          FOR n = i_q, N_ELEMENTS(q) - 1 DO BEGIN
            ; Define as a variable
            gradient = q[n]


;            ; Loop #7: v_doppler
;            FOR p = i_dop, N_ELEMENTS(v_doppler) - 1 DO BEGIN
;              ; Define as a variable
;              doppler = v_doppler[p]
;
;
;              ; Loop #8: v_turb
;              FOR r = i_turb, N_ELEMENTS(v_turb) - 1 DO BEGIN
;                ; Define as a variable
;                turbulence = v_turb[r]

                ;;;;;;;;;;;;;;;;;;;;;;;;;;;; DEFINE Line Profile (FUNCTION) ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ; Now, we have all the parameters necessary to call the line profile function.
                ; The line profile will depend on v_obs and T_rH2 arrays, as well as other parameters.
                ;
                ;;;;;; Define rH2, phi, sighltine, v_obs, T_rH2, and line_profile through GRID_DISKS ;;;;;;
                grid_disks_v2, r_in, r_mid, r_out, AUsize_inner, AUsize_outer, phi_density, distance, restwave, transition_wave, starmass, inclination, temperature, gradient, turbulence, doppler, r_grid, phi_grid, sightline, v_obs, T_rH2;, H2_line_profile
;                grid_disks, r_grid, r_H2, phi_grid, phi_H2, distance, restwave, transition_wave, starmass, inclination, temperature, gradient, turbulence, doppler, sightline, v_obs, T_rH2, H2_line_profile
                ;grid_disks, AUsample, r_in, distance, restwave, transition_wave, starmass, inclination, temperature, gradient, turbulence, doppler, r_H2, phi, sightline, v_obs, T_rH2, H2_line_profile
                ;LyA_flux_unred = unred_flux_to_H2(start_flux_red, transition_wave, Av, r_H2, distance, r_LyA)

                LyA_flux_unred = unred_flux_to_H2(start_flux_red, LyA_pumping, Av, r_grid, distance, r_LyA)


;                H2_line_profile = H2_line_profile[*,0]    ; only the v=0 component of the line profile is necessary



                ;;;;;;;;;;;;;;;;;;;; DEFINE Exponential Temperature [Partition Function] (FUNCTION) ;;;;;;;;;;;;;;;;;;;;;;;
                ; We can define the exponential temperature relation with the H2 temperature profile
 ;               energy_distribution = exp_temp(T_rH2, restwave)
                transition_energy = trans_temp(T_rH2,lower_energy)
;                stop


                ;;;;;;;;;;;;;;;;;;;; DEFINE Scale Height [Hp] (FUNCTION) ;;;;;;;;;;;;;;;;;;;;;;;
                ; We can define the scale height as a function of r_H2, stellar mass, and temperature at z = 0.
                ;Hp = scale_height(T_rH2, r_H2, starmass, phi)
                Hp = scale_height(T_rH2, r_grid, starmass, phi_grid)

                ;;;;;;;;;;;;;;;;;;;; DEFINE Exponitial Density Distribution (FUNCTION) ;;;;;;;;;;;;;;;;;;;;;;;
                ; We can define the exponential density distribution relation with the scale height at each r_H2
                ; NOTE: FOR NOW, WILL ASSUME z=0 AT ALL LOCATIONS ITHE DISK, SO IGNORE!
                ; (Because rho = 1 for exp( 0 ) = 1)
                rho_disk = exp_disk_density( Hp, disk_height )
       ;         rho_disk = EXP(-0.5*(5)^2)      ; Since I'm using scale height as a means to define tau, the contributing H2 should be a few scale heights above the plane
;                rho_disk = 1.0

                ;;;;;;;;;;;;;;;;;;;;; DEFINE Cross Section (FUNCTION) ;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ; Once we have the line profile created, we can calculate the cross section of the H2
                xsection = cross_section(sig0)
;                plot,v_obs,xsection,xr=[-50,50],xs=1,xtitle='v_obs',ytitle='X-section of emitting H2'
;                stop


                ;;;;;;;;;;;;;;;;;; DEFINE partition function ;;;;;;;;;;;;;;;
                ; Use Kevin's h2pop code to create an array of partition functions at each
                ; location in the disk.
                partition_fnct = partition_function(T_rH2,J_rot_state,T_H2_range,Z_TH2)




                ; Loop #9: r_char
                FOR s = i_rch, N_ELEMENTS(r_char) - 1 DO BEGIN
                  ; Define as a variable
                  characteristic_radius = r_char[s]

                  ; Define the surface density of the disk
                  ;rho_surface = surface_density(r_H2, characteristic_radius, phi)
                  rho_surface = surface_density(r_grid, characteristic_radius, density_gamma)
;                  stop

                  ; Loop #10: Total H2 Mass contributing to the emission line
                  FOR u = i_mh2, N_ELEMENTS(M_H2) - 1 DO BEGIN
                    ; Define as a variable
                    totalH2mass = M_H2[u]

                    ; We can now determine the rest of our intrinsic/thermal variable to then
                    ; determine the radiative transfer behavior of the line (i.e. get the intensity)

                    ; First, get the number density of H2 in the disk (with the transition of restwave)
                    n_H2 = H2_number_density(totalH2mass, Hp, characteristic_radius, J_initial, rho_surface, rho_disk, density_gamma, transition_energy, partition_fnct)
                   ; plot,r_grid,n_H2,xs=1,xr=[0.5,2],xtitle='r [au]',ytitle='number density of emitting H2'
;                    stop

                    ; Now we have all the information to calculate the opacity of the H2
                    tau_z = optical_depth( n_H2, xsection, Hp, H2_line_profile, disk_height, T_rH2, column_tau, T_tau, eff_tau )
;                    plot,r_grid,tau_z,xs=1,xr=[0.5,2],xtitle='r [au]',ytitle='opacity of the emitting H2'
;                    stop

                    ; Calculate intensity
                    ;;;;; CHANGED - return the intensity bin for ALL
                    intensity = total_intensity(tau_z,sightline,LyA_flux_unred,temperature,v_obs,v_obs_bin,AUsample,r_grid,Hp,distance,Av,starmass,inclination,turbulence,characteristic_radius,totalH2mass,gradient,density_gamma,disk_height)

                    ; Calculate intensity of each vobs bin
                    intensity_binned = all_emission_line_intensities(Av, v_obs, v_obs_bin, intensity, lambda_ref, B_ul, J2, J_progression)
           ;         stop




                    ; LAST PART!!
                    ; Writing the files:
                    ; Order: Directory --> model\target\wavelength
                    jname = STRTRIM(j,2) ;STRTRIM(j,2)    ; disk height
;                    kname = STRTRIM(k,2)
                    lname = STRTRIM(l,2) ;STRTRIM(l,2)    ; gamma
                    mname = STRTRIM(m,2) ;STRTRIM(m,2)    ; T_1au
                    nname = STRTRIM(n,2) ;STRTRIM(n,2)    ; q
    ;                pname = STRTRIM(0,2) ;STRTRIM(p,2)    ; v_doppler
    ;                rname = STRTRIM(1,2) ;STRTRIM(r,2)    ; v_turb
                    sname = STRTRIM(s,2) ;STRTRIM(s,2)    ; r_char
                    uname = STRTRIM(u,2) ;STRTRIM(u,2)    ; M_H2
;                    savefile=ref_path_models+filedir[target]+wavename+'\'+starname+'-'+jname+'-'+lname+'-'+mname+'-'+nname+'-'+pname+'-'+rname+'-'+sname+'-'+uname+'.sav'
                    savefile=ref_path_models+filedir[target]+'disk_profiles_v2'+'\'+starname+'-'+jname+'-'+lname+'-'+mname+'-'+nname+'-'+sname+'-'+uname+'.sav'
                    print,'Saving to file = '+savefile
          ;          SAVE, starname, lambda_ref, LyA_pumping, temperature, gradient, density_gamma, disk_height, characteristic_radius, totalH2mass, v_obs_bin, intensity_binned, filename=savefile ;r_grid, intensity,
                    print,SYSTIME()

                    intensity = intensity*phi_density

                    ; plot a few things, including flux distribution of H2 by radius, and resulting line profiles
                    ; find max values in both intensity arrays, plot with ymax = max
                    mx_int_r = max(intensity)
                    mx_int_v = max(intensity_binned[0:598,*])

                    ; make intensity v. radius in one window
                    WINDOW,0
                    th = 6
                    ln = 0
           ;         set_plot, 'ps'
           ;         device, /color, /times, /isolatin1, /landscape, XOFFSET=0, YOFFSET=30, XSIZE=30.0, YSIZE=19.0,$
           ;               filename=ref_path_models+'zzz_disk_profiles_v2_parm_effects\'+$
           ;               'T_H2='+STRTRIM(FIX(temperature),2)+$
           ;               '-q='+STRTRIM(FIX(gradient),2)+$
           ;               '-r_c='+STRTRIM(FIX(characteristic_radius),2)+$
           ;               '-gamma='+STRTRIM(FIX(density_gamma+0.1),2)+$
           ;               '-M_H2='+STRTRIM(FIX(alog10(totalH2mass)),2)+$
           ;               '-z='+STRTRIM(FIX(disk_height),2)+'.ps'


                    plot,r_grid,intensity[*,0,0],xr=[0.01,100],yr=[0,mx_int_r*1.2],ys=1,font=0,xtitle='radius (AU)',$
                      ytitle='flux distribution by H2 absorption progression',thick=th,linestyle=ln,/xlog
                    oplot,r_grid,intensity[*,0,1],color=!red,thick=th,linestyle=ln
                    oplot,r_grid,intensity[*,0,2],color=!dyellow,thick=th,linestyle=ln
                    oplot,r_grid,intensity[*,0,3],color=!green,thick=th,linestyle=ln
                    legend_astron,['T_H2 = '+STRTRIM(temperature,2)+' K',$
                      'q = '+STRTRIM(gradient,2),$
                      'r_c = '+STRTRIM(characteristic_radius,2)+' AU',$
                      'gamma = '+STRTRIM(density_gamma,2),$
                      'M_H2 = '+STRTRIM(totalH2mass,2)+' M_sol',$
                      'z/r = '+STRTRIM(disk_height,2)],/right_legend
                    legend_astron,['[0,1]','[0,2]','[1,4]','[1,7]'],color=[!white,!red,!dyellow,!green],$
                      linestyle=[ln,ln,ln,ln],thick=[th,th,th,th],/left_legend


                    ; make intenstiy_binned v. v_obs_bin in another window
                    WINDOW,1
                    J2 = [1,2,2,1,2,4,1,2,7,4,7,4,1,7,7]
                    plot,v_obs_bin,intensity_binned[*,0],psym=10,xr=[-50,50],yr=[0,mx_int_v],ys=1,font=0,$
                      xtitle='velocity (km s-1)',ytitle='flux (ergs cm-2 s-1 Ang-1)',thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,1],psym=10,color=!red,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,2],psym=10,color=!dred,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,3],psym=10,color=!orange,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,4],psym=10,color=!dorange,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,5],psym=10,color=!yellow,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,6],psym=10,color=!lgreen,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,7],psym=10,color=!green,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,8],psym=10,color=!lblue,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,9],psym=10,color=!blue,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,10],psym=10,color=!cyan,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,11],psym=10,color=!dcyan,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,12],psym=10,color=!purple,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,13],psym=10,color=!pink,thick=th,linestyle=ln
                    oplot,v_obs_bin,intensity_binned[*,14],psym=10,color=!dpink,thick=th,linestyle=ln
                    legend_astron,['T_H2 = '+STRTRIM(temperature,2)+' K',$
                      'q = '+STRTRIM(gradient,2),$
                      'r_c = '+STRTRIM(characteristic_radius,2)+' AU',$
                      'gamma = '+STRTRIM(density_gamma,2),$
                      'M_H2 = '+STRTRIM(totalH2mass,2)+' M_sol',$
                      'z/r = '+STRTRIM(disk_height,2)],/right_legend

                    ; try fitting a gaussian to the intensities, and return fwhm
                    mx = DBLARR( 15 )
                    FOR i = 0, 14 DO mx[i] = max(intensity_binned[0:598,i])
                    fwhm = mx / 2.
                    i_fwhm = INTARR( 15 )
                    FOR i = 0, 14 DO BEGIN
                      i_temp = WHERE(fwhm[i] LE intensity_binned[0:300,i]+(fwhm[i]*0.1) AND fwhm[i] GE intensity_binned[0:300,i]-(fwhm[i]*0.1))
                      i_fwhm[i] = i_temp[0]
                    ENDFOR
                    vline,v_obs_bin[i_fwhm],linestyle=1,thick=th

                    WINDOW,3
                    plot,r_grid,n_H2[*,0,0],linestyle=ln,thick=th,xr=[0.01,100],xs=1,yr=[1.E-2,1.E+13],ys=1,$
                      /ylog,/xlog,$
                      xtitle='radius (AU)',ytitle='number density (cm-3)'
                    oplot,r_grid,n_H2[*,0,1],linestyle=ln,thick=th,color=!dred
                    oplot,r_grid,n_H2[*,0,2],linestyle=ln,thick=th,color=!orange
                    oplot,r_grid,n_H2[*,0,3],linestyle=ln,thick=th,color=!lgreen
                    legend_astron,['[0,1]','[0,2]','[1,4]','[1,7]'],color=[!white,!dred,!orange,!lgreen],$
                      linestyle=[ln,ln,ln,ln],thick=[th,th,th,th],/left_legend
  ;                  device, /CLOSE_FILE
  ;                  set_plot,'WIN'

                    stop
                  ENDFOR   ; ends M_H2 FOR-loop
                ENDFOR   ; ends r_char FOR-loop
;              ENDFOR   ; ends v_turbulent FOR-loop
;            ENDFOR   ; ends v_doppler FOR-loops
          ENDFOR   ; ends temperature gradients (q) FOR-loop
        ENDFOR   ; ends temp at 1 AU (T_1AU) FOR-loop
      ENDFOR   ; ends gamma (surface density gradient) FOR-loop
;    ENDFOR   ; ends disk inclination (i_disk) FOR-loop
  ENDFOR   ; ends scale height (Hp_height) FOR-loop
; ENDFOR   ; ends wavelength FOR-loop

ENDFOR   ; end TARGET FOR-loop

print, 'End of program'




END
