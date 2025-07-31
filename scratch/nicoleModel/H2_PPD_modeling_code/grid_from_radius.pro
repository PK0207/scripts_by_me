; Try radial logarithmic scale + linear phi scale, then figure out [x,y] points for
; each [r,phi] location
; 

;pro grid_from_radius
pro grid_from_radius, inclination, distance, rmin, rmid, rmax, rin_ndensity, rout_ndensity, phi_ndensity, xgrid, ygrid, r_grid, phi_grid, r_H2, phi_H2

;inclination = 0.0 
;distance = 100.0    ; pc
;
;rmin = 0.1
;rmid = 1.0
;rmax = 10.0
;
;rin_ndensity = 10
;rout_ndensity = 10
;phi_ndensity = 360

; Call the scale_vector function from coyote
rinner = scale_vector(FINDGEN(rin_ndensity),rmin,rmax)
router = scale_vector(FINDGEN(rout_ndensity),rmid,rmax)

rinout = 10.^(0.02*FINDGEN(200) - 1.4)

; Define r_grid from the above scale vector quantities
r_grid = rinout / distance

; Define phi at every 1 degree increment
phi_grid = scale_vector(FINDGEN(phi_ndensity),0.0,2*!pi)

; Make x,y grid empty arrays of the size of rows = radius, columns = phi
; and fill in by solving for x, y grid spaces
xgrid = DBLARR( N_ELEMENTS(r_grid),N_ELEMENTS(phi_grid) )
ygrid = DBLARR( N_ELEMENTS(r_grid),N_ELEMENTS(phi_grid) )

; Solve for tan(phi)
tan_phi = TAN(phi_grid)

; FOR loops- loop through radial and azimuthal positions on the disk,
; then solve for x, y. y will be corrected by a factor of cos(inclination) before
; being stored in ygrid
FOR radial = 0, N_ELEMENTS(r_grid) - 1 DO BEGIN
  FOR azimuth = 0, N_ELEMENTS(phi_grid) - 1 DO BEGIN
    ; 1. x-grid
    IF (phi_grid[azimuth] GE 90.0*!pi/180.) AND (phi_grid[azimuth] LT 270.*!pi/180.) THEN BEGIN
      x_dummy = -1.0*SQRT( (r_grid[radial]*r_grid[radial]) / (1.0 + tan_phi[azimuth]*tan_phi[azimuth]) ) 
    ENDIF ELSE BEGIN
      x_dummy = SQRT( (r_grid[radial]*r_grid[radial]) / (1.0 + tan_phi[azimuth]*tan_phi[azimuth]) )
    ENDELSE
    xgrid[radial,azimuth] = x_dummy / COS(inclination*!pi/180.)
    
    ; 2. y, then y-grid
    IF phi_grid[azimuth] GE (180.0*!pi/180.) THEN BEGIN
      y_dummy = -1.0*SQRT( (r_grid[radial]*r_grid[radial]) / (1.0 + (1.0/(tan_phi[azimuth]*tan_phi[azimuth]))) )
    ENDIF ELSE BEGIN
      y_dummy = SQRT( (r_grid[radial]*r_grid[radial]) / (1.0 + (1.0/(tan_phi[azimuth]*tan_phi[azimuth]))) )
    ENDELSE
    ygrid[radial,azimuth] = y_dummy ;/ COS(inclination*!pi/180.)
    
  ENDFOR
ENDFOR

; Now recalculate radial and azimuthal positions of the gas disk
r_H2 = DBLARR( N_ELEMENTS(r_grid),N_ELEMENTS(phi_grid) )
phi_H2 = DBLARR( N_ELEMENTS(r_grid),N_ELEMENTS(phi_grid) )

FOR x = 0, N_ELEMENTS(r_grid) - 1 DO BEGIN
  FOR y = 0, N_ELEMENTS(phi_grid) - 1 DO BEGIN
    ; 1. radial position
    r_dummy = SQRT(xgrid[x,y]*xgrid[x,y] + ygrid[x,y]*ygrid[x,y]) * distance
    IF r_dummy LT 40.0 THEN r_H2[x,y] = r_dummy ELSE r_H2[x,y] = r_dummy;0.0
    
    ; 2. azimuthal position
    phi_dummy= ATAN(ygrid[x,y] / xgrid[x,y])   ; NEED TO FIX FOR CERTAIN REGIONS OF x,y
    
    ; Test to see where where we are in the disk, and how phi needs to be
    ; shifted to fit our idea of the model (phi = 0, 180 at y = 0
    IF ygrid[x,y] LE 0.0 THEN BEGIN
      IF xgrid[x,y] LE 0.0 THEN BEGIN
        phi_H2[x,y] = phi_dummy + !pi
      ENDIF ELSE BEGIN
        phi_H2[x,y] = phi_dummy + 2.0*!pi
      ENDELSE
    ENDIF ELSE BEGIN
      IF xgrid[x,y] Le 0.0 THEN BEGIN
        phi_H2[x,y] = phi_dummy + !pi
      ENDIF ELSE BEGIN
        phi_H2[x,y] = phi_dummy
      ENDELSE
    ENDELSE
    
  ENDFOR
ENDFOR



RETURN

END