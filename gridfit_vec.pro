;***********************************************************
; Bob Zavala me fecit
; email: robert.t.zavala.civ@mail.mil
;***********************************************************

;;DESCRIPTION
; A procedure to test for the presence of a binary in NPOI
; data exported from an OYSTER chameleon file. Incorporates 
; a file selection widget and a GUI for adjusting parameters
; used in the gridfit search.
;
; For each delta-mag a
; minimum chi-squared is found, and the global
; minimum for all the delta-mags is found. This is
; a vectorized version of cubegrid.pro and the vectorization
; owes a debt of gratitude to my Marist High School junior
; year physics teacher Mr. Stankevitz.
;
;;SOME HISTORY
; Added a GUI to make grid parameter changing easier.
; In the past I had to edit the code and re-compile. 
;
; Added a file selection widget to eliminate entering the 
; .dat filename at the command line.
;
; Modification of cubegrid to handle files created with the
; sData.pro procedure. This will grab the diameters from the
; saved structure instead of requiring them to be input
; by hand. The diameters are obtained in doublesearch
; procedure below, and d2 is set to 0.5*d1.
; 
; The structure is also expanded to hold some additional star
; information.
; 
; Note that if riflag returns True than the diameter was not
; in diameter.bsc and the diameter was computed using the R-I
; color as a OYSTER procedure call within sData.pro

;;RELATED PROGRAMS
; sData.pro: wrote out the visibility data used here
;
; doubleSearch: procedure below does the grid search.
;
; alldata.pro: Now used to write out the visibility data
;            from the Chameleon file. Replaced sdata.pro
;***********************************************************
PRO gridfit_vec

; Procedure to conduct a grid search test for a binary in NPOI 
; visibility data. Data are searched in a vectorized manner 
; that spreads work across multiple CPU's. Limitation is then 
; memory, and an available memory check is made. Widegts are 
; used to select the input data file and to set parameters 
; necessary for the grid search and for output plots that 
; are requested by the user. 

COMPILE_OPT IDL2, STRICTARRSUBS

; Common statements

COMMON VISDATA, datafile

; Use file selection widget to select input visibility
; data file.

datafile = DIALOG_PICKFILE(/READ, FILTER='*.dat')

; Call the gridfit parameter setting widget. This also 
; contains a button to start the gridfit code. 

tab_grid_widget

END

;***********************************************************

function ww_buttoncall, ev

; Button handling function for the GRIDFIT widget created 
; by tab_grid_widget EXCEPT GRIDFIT as that is handled 
; by an EVENT_PRO call. 

COMMON VISDATA

WIDGET_CONTROL, ev.TOP, GET_UVALUE=stash

; CASE statements to select based on event based by button

CASE ev.ID OF
  
  ; Select a new/next input .dat file
  stash.bDatFile: datafile = DIALOG_PICKFILE(/READ, FILTER='*.dat')
    
  ; Show the Help file
  stash.bHelp: XDISPLAYFILE,'gridfit_help.txt',TITLE='GRIDFIT_VEC HELP'   
    
  ; User is DONE
  stash.bDone: WIDGET_CONTROL, ev.TOP, /DESTROY    

ENDCASE

end

;***********************************************************

pro cubegridm_vec, event

;;HOW TO USE
; 0) variable names enclosed in ()
; 1) The input datafile is provided via a COMMON statement
; 2) A widget and COMMON statement is sued to set parameters
;    used in the grid search.
; 3) The same widget in 2) is used to control which output
;    lots are made for the user.
; 4) A test is built into this procedure so that the requested
;    search area will not be allowed to exceed available memory.

;;OUTPUT
; You will get a postscript file with the name (outprefix)+.ps
; and an idl .dat file containing a structure (gridout) which
; has the gridding output results for each delta-mag step
; plus soime general inofrmation on the star.
;
;;PERFORMANCE
; Checking against cubergidm with same search area in x and y
; and same pixelSize and deltamag intervals gives a factor of
; 2 improvement in time of cubegridm_vec versus cubegridm.
; test conducted "back in the day" on alex circa 2010.
; 
;;MEMORY USAGE
; Memory required is measured against the data cube that 
; is sent via vectorization for the GRIDFIT. The cube axes 
; are:
;     relative R.A.  (mas) pixels set by pixelSize
;     relative Decl. (mas) pixels set by pixelSize
;     V^2 points sampled by number of u-v samples
;     
; Knowing these 3 axes, you have the number of data points 
; as a volume. The memory required is:
;                80.09 Bytes/(data point)
;
;************************************************************

COMMON VISDATA
COMMON GRIDPARAMS, xmin, xmax, ymin, ymax, pixSize, $
                numdmag, makeplots_q, makedMagPlots_q
COMMON SHAREIT, diam, sep, pa, theMin, theMax, sym_Delta, $ 
                sym_ChiSq, starname, npoiid, riflag, $
                x_min, x_max, y_min, y_max, pixelSize, $
                n_x, n_y, n_total, num_dmag, outprefix, $
                obsDate, makeplots, makedMagPlots

COMPILE_OPT IDL2, STRICTARRSUBS

; Initialize datafile, number of deltamag iterations, output strings
; and whether to make plots

; Now set the dmag steps. Note that as the increment 
; is 0.1 as indicated in the next FOR loop the 
; largest delta-mag we examine is num_dmag*0.1

; Error trap here to make sure datafile is a file and not just 
; a path. A user who does not select a file via DIALOG_PICKFILE 
; will have a valid path for datafile but no actual file. This is 
; a valid result from DIALOG_PICKFILE but will cause a crash in 
; pro doublesearch 

fileInform = FILE_INFO(datafile)
fileTest = fileInform.regular

IF fileTest THEN BEGIN
    print,'In cubegridm_vec datafile is: ',datafile
ENDIF ELSE BEGIN
    errMessage=['Required input .dat file is missing',$ 
               'Please click on Datfile file to select an input file'] 
    z=DIALOG_MESSAGE(errMessage, /ERROR, /CENTER)
    return
ENDELSE

; Restore the input file so we can get the number of 
; data points to use in memory estimate

; RESTORE, FILENAME=datafile

; I will use the IDL_Savefile object for the restore 
; process for practise

fObj = OBJ_NEW('IDL_Savefile',datafile)
fObj->Restore,'sData'

num_dmag=numdmag*10
outprefix = STRMID(datafile,0,STRLEN(datafile)-4)
obsDate = STRMID(datafile,strlen(datafile)-1-13,10)
makeplots = makeplots_q
PRINT,'Will cubegridm make plots? ___',makeplots,'___'

; Determine memory available and maximum size of an area to search
; for this computer. This returns a string array with the required 
; information in freeOutput[1]

SPAWN,'free',freeOutput
memoryLine = STRSPLIT(freeOutput[1],/EXTRACT)
; Now get maximum memory in kilobytes, convert to bytes
availableMemory=LONG(memoryLine[6])
PRINT,' '
PRINT,'Available memory is: '+STRING(availableMemory)+' KiloBytes' 
PRINT,' '

; Select the range in mas of x and y space to search
; Searching the inner 60 mas is ~2x resolution of 4meter
; speckle.
; We search only one half of the u-v plane.
; Now set the following via the COMMON statement transferred 
; variables.  
x_min = xmin
x_max = xmax
y_min = ymin
y_max = ymax
pixelSize = pixSize

n_x = FLOOR((x_max - x_min)/pixelSize) + 1
n_y = FLOOR((y_max - y_min)/pixelSize) + 1
n_total = n_x * n_y

; Check that requested region is not too big
; and on a 64 bit machine we need 80.09 bytes/pixel/data point
num_data_points = N_ELEMENTS(sData.v2cd)
maxMemoryNeeded = 80.09*n_total*num_data_points

IF maxMemoryNeeded/1000 GE availableMemory THEN BEGIN
    print,'*****************************************'
    print,maxMemoryNeeded/1000,FORMAT='("Maximum memory needed is: ",5x,I12," KiloBytes")'
    print,'and exceeds available memory of '+memoryLine[6]+'000 Bytes.'
    print,'Stopping, recommend you resize the gridding'
    print,'region or obtain more memory.'
    print,'*****************************************' 
    print,' '
    RETURN
ENDIF ELSE BEGIN 
    print,'*****************************************'
    print,maxMemoryNeeded/1000,FORMAT='("Memory needed is: ",3x,I12," KiloBytes")'
    print,'*****************************************' 
    WAIT, 15
ENDELSE    


; We are finished with the savefile datafile object and the sdata copy
fObj->Cleanup 
DELVAR,sdata

;
; Define structure to contain deltamag, sep, pa, min and max ChiSq 
; and minimum flag for each deltamag call to doubleSearch
; Note that the first 4 entries are like a header and are single
; elements, not arrays.

gridout = {grid, starname:' ', npoiid:' ', date:'', diam:0.0, riflag:' ', $ 
           deltamag:FLTARR(num_dmag+1), sep:FLTARR(num_dmag+1), $
           pa:FLTARR(num_dmag+1), theMin:FLTARR(num_dmag+1), $ 
           theMax:FLTARR(num_dmag+1)}

; Initialize start time and loop calls to doubleSearch
; Note we are looping over 0.1 mag steps.

start_time = SYSTIME(1)
FOR i=0,num_dmag DO BEGIN
    deltamag=0.1*i
    doubleSearch,datafile,deltamag,makedMagPlots_q
    gridout.deltamag[i]=deltamag
    gridout.sep[i]=sep
    gridout.pa[i]=pa
    gridout.theMin[i]=theMin
    gridout.theMax[i]=theMax
ENDFOR
stepps = STRCOMPRESS(STRING(num_dmag+1))


; Find minimum, report it later with other areas

globalMin=MIN(gridout.theMin,minLoc)


; Only make plot if we really want it
; If we do make a single page with three plots on it

IF makeplots EQ 'yes' THEN BEGIN
    USERSYM, COS(FINDGEN(17)*(!PI*2/16.0)), SIN(FINDGEN(17)*(!PI*2/16.0))
    set_plot,'PS'
    !P.MULTI=[0,1,3]
    DEVICE,filename=outprefix+'.ps',xsize=7.25,ysize=10.0,xoffset=0.5, $ 
           yoffset=0.5,/inches,FONT_SIZE=20
    PLOT,gridout.deltamag,gridout.theMin,XTITLE=sym_Delta+'!6mag',YTITLE=sym_ChiSq, $
      TITLE='!6'+datafile,CHARTHICK=2.0,THICK=3.0,XTHICK=3.0,YTHICK=3.0, $
      PSYM=-8, SYMSIZE=0.5, XRANGE=[-0.1,num_dmag/10], XSTYLE=1, $ 
      SUBTITLE='High res: Global minimum at '+sym_Delta+'mag ='+STRING(gridout.deltamag[minLoc],FORMAT='(F5.2)')$
      +' '+sym_Chisq+' = '+STRING(globalMin,FORMAT='(F7.2)')
    OPLOT,[gridout.deltamag[minLoc],gridout.deltamag[minLoc]],[0,max(gridout.theMax)],THICK=3.0, $
          LINESTYLE=2
    USERSYM, COS(FINDGEN(17)*(!PI*2/16.0)), SIN(FINDGEN(17)*(!PI*2/16.0)),/FILL
    PLOT,gridout.deltamag,gridout.sep,XTITLE=sym_Delta+'mag',YTITLE='Separation (mas)', $
      CHARTHICK=2.0,THICK=3.0,XTHICK=3.0,YTHICK=3.0, $
      PSYM=8, SYMSIZE=0.5, XRANGE=[-0.1,num_dmag/10], XSTYLE=1, $ 
      SUBTITLE='Global minimum at sep (mas) ='+STRING(gridout.sep[minLoc],FORMAT='(F6.2)')
    OPLOT,[gridout.deltamag[minLoc],gridout.deltamag[minLoc]], $ 
          [0,max(gridout.sep)+0.2*max(gridout.sep)], $ 
          THICK=3.0, LINESTYLE=2
    PLOT,gridout.deltamag,gridout.pa,XTITLE=sym_Delta+'mag',YTITLE='PA (deg)', $
      CHARTHICK=2.0,THICK=3.0,XTHICK=3.0,YTHICK=3.0, $
      PSYM=8, SYMSIZE=0.5, XRANGE=[-0.1,num_dmag/10], XSTYLE=1, $ 
      SUBTITLE='Global minimum at PA (deg) ='+STRING(gridout.pa[minLoc],FORMAT='(F6.2)')
    OPLOT,[gridout.deltamag[minLoc],gridout.deltamag[minLoc]],[0,200],THICK=3.0, $
          LINESTYLE=2
    DEVICE,/close_file
    set_plot,'X'    
ENDIF

; Set the returned header-like information in output structure
gridout.starname = starname
gridout.npoiid = npoiid
gridout.date = obsDate
gridout.diam = diam
gridout.riflag = riflag

; save the grid output file
save,gridout,    filename=outprefix+'.gridout_vec.sav'

;get ending time
totaltime = STRCOMPRESS(STRING(SYSTIME(1)-start_time))

;Tell user how long this grid test took.
PRINT, 'Total time to for grid test = '+totaltime+' sec'

print,'********************'
print,'For star ',starname,' initial binary grid search finds:'
print,'Minimum reduced chi-squared of ',globalMin,' found at deltaMag of ',gridout.deltamag[minLoc]
print,'********************'

print,'**********************************'
print,'Grid search for a binary completed.'
print,'Data saved to an idl files called '+outprefix+'.*.gridout.sav'
print,'Plots saved as '+outprefix+'.ps '
print,'**********************************'
;
end

;********************************************************************

PRO doubleSearch,datafile,deltamag,makedMagPlots_q

COMPILE_OPT DEFINT32, STRICTARR, STRICTARRSUBS

COMMON SHAREIT

; Modification of Chris Tycner's binarySearch.pro so that the 
; channel widths are properly included in the binary star 
; visibility model equation.

; I also changed it so that a FOR loop is vectorized.

; RELATED PROGRAM: sData.pro, cubegridm (above)
; sData procedure reads in the data from OYSTER CHA files and 
; puts it into a struture appropriate for this and other tasks. 
; This exists as a stand-alone procedure and is a modification 
; of sdata

; cubegridm (above) is the controlling program for the grid 
; search

; CHANNEL WIDTHS problem: Some cha files from the TPF survey have 
; a zero channel width. This creates NaN entries when we calculate
; the sinc_factor. I will make a test for zero channel widths 
; and impose a model channel width based on a fit to the 
; the channel widths. A 2nd order polynomial fit gives:
; chanwidths =
; 3.33149e-08-2.35728e-09*(channel_number)+6.24985e-11*(channel_number)^2
; I probably could have hard coded the channel widths but I 
; got fancy with this fit.

; PLOTTING VARIABLES/FLAGS
; psOutput controls whether postscript or X window plots of 
; individual delta-mag results are made. 
;
; makedMagPlots controls whether we make those plots at all.
; You can probably skip this unless you want to visualize the 
; chi-squared space.

pi = !DPI
; make sure the multi flag from cubegridm doesn't impact us here.
!P.MULTI=0

!EXCEPT=0
; Report to user the line number of any math errors. You might get one
; at line containing sincfactor if you pass through the origin. 
; We allow this error and don't fix it to aviod IF loops and tests 
; which will slow the code. Also, if we eventually vectorize entire
; code for speed, we won't be able to use any IF statements or other tests. 

; Define plot symbol Delta,ChiSq

sym_Delta='!7D!X!N' 
sym_ChiSq='!7v!S!E!72!R!B!I!7m!N!X'

LOADCT, 3, /SILENT

; psOutput = 0 makes plots to screen, 1 makes ps plots
psOutput = 1

; Now set a flag for whether to make plots or not
; Set yes or no in the widget and passed in the 
; input parameter
makedMagPlots=makedMagPlots_q

;spawn,'pwd',here
;here=here+'/'
;here = here+'/NPOI/TPF/'
;file=here+datafile
RESTORE, FILENAME=datafile

pick_ch = 30 ; select channel to exclude
index = WHERE(sData.channel NE pick_ch)
PRINT, 'Excluding ' + STRING(N_ELEMENTS(sData.v2cd)-N_ELEMENTS(index)) + ' data points.'
PRINT, 'Keeping ' + STRING(N_ELEMENTS(index)) + ' data points.'

num_datapts=N_ELEMENTS(index)

print,'num_datapts is: ',num_datapts

; put the u, v coordinates, calibrated v^2 values into vectors.
; make a vector of the uv radii, channel widths and 
; channel effective wavelengths
u = sData.u[index]
v = sData.v[index]
v2 = sData.v2cd[index]
v2_err = sData.v2cd_err[index]
uvrad = SQRT(u*u+v*v)
chanw = sData.chan_width
lam_eff = sData.uv_wave

; Check for any zero channel widths 
zeroChanInd = WHERE(chanw EQ 0, zerocount)

; Replace zero channel widths with the parabolic fit estimate
IF zerocount NE 0 THEN BEGIN
    chanw[zeroChanInd] = 3.33149e-08 $ 
                         - 2.35728e-9*sData.channel[zeroChanInd] $ 
                         + 6.24985e-11*sData.channel[zeroChanInd]^2
ENDIF

; calculate the PA of each u-v point (measured E from N)
basePA = 90.0d - 180.0d/pi*ATAN(v,u) ; in deg

stringy = STRMID(STRCOMPRESS(STRING(pixelSize)),0,4)
PRINT, 'Requested region is ' + STRING(n_x) + ' by ' + STRING(n_y)
PRINT, 'I will increment by '+stringy+' mas.'
PRINT, 'with total number of elements = ' + STRING(n_total)

; set up the x and y vectors
x = DBLARR(n_total)
y = DBLARR(n_total)
counter = 0L
FOR i = 0, n_x-1 DO BEGIN
    FOR j = 0, n_y-1 DO BEGIN
        x[counter] = x_min + DOUBLE(i)*pixelSize
        y[counter] = y_min + DOUBLE(j)*pixelSize
        counter += 1
    END
END

; calculate the PA of x-y point
; Note: x is E and y is N
xyPA = 90.0d - 180.0d/pi*ATAN(y,x)   ; in deg

; calculate separations for each x-y point in radians
xySep = SQRT(x^2 + y^2)/1000d/206265d

; chose binary parameters
; ratio = brightness ratio, d1 and d2 are diameters of
; stars 1 and 2 in mas. mas_to_rad is conversion factor
; mas -> radians
ratio = 10^(-1*deltamag/2.5)
mas_to_rad = 1/(206265*1000.)
diam = sData.diam
starname = sData.name
riflag = sData.riflag
npoiid = sData.npoiid 
d1=sData.diam
d2=d1*0.5
d1rad = d1*mas_to_rad
d2rad = d2*mas_to_rad
; Calculate V^2 vector for the two stars
arg_1 = pi*uvrad*d1rad
arg_2 = pi*uvrad*d2rad
vstar1 = 2*beselj(arg_1,1)/arg_1
vstar2 = 2*beselj(arg_2,1)/arg_2

; Start gridding

PRINT, 'Gridding chi^2 space...'

; calculate ch^2 for each x-y point
chi2 = DBLARR(n_total)

;
; TO VECTORIZE THIS K LOOP
; first_term = xySep*cos(!pi*xyPA/180d)#cos(!pi*basePA/180d)
; second_term = xySep*sin(!pi*xyPA/180d)#sin(!pi*basePA/180d)
; proj_sep = first_term + second_term 
; problem occurs with sinc_arg. 
; sinc_arg_right = !pi*chanw*uvrad/lam_eff gives a [31] and not a [1,31] dimension
; fix that reform(array, 1, n_elements(chanw)) and get the sinc_arg
; with correct dimensions? 

; convert position angles to radians
basePA = pi*basePA/180d
xyPA = pi*xyPA/180d
; create array of proj_sep using cos(a-b) identity
first_term = cos(xyPA)#cos(basePA)
second_term = sin(xyPA)#sin(basePA)
; resize the xySep 
xySep = xySep#(make_array(num_datapts,VALUE=1.0)) 
proj_sep = xySep*(first_term + second_term)
; clear up some memory by setting first_term and 
; second_term to a single integer. Note delvar doesn't work here
first_term = 1
second_term = 1
; make an array of 1.0
ones=xyPA/xyPA
sinc_arg_1 = ones#(pi*chanw*uvrad/lam_eff)
sinc_arg = sinc_arg_1*proj_sep
; get some back memory again, sinc_arg_1
sinc_arg_1 = 1
sinc_factor = sin(sinc_arg)/sinc_arg
; get some memory back again, again: sinc_arg
sinc_arg = 1
s = ones#uvrad
vstar1Sqrd = vstar1^2
vstar1     = ones#vstar1
vstar2     = ones#vstar2
vstar1Sqrd = ones#vstar1Sqrd
v2         = ones#v2
v2_errSqrd = ones#v2_err^2
v2binary_model = $
       (vstar1Sqrd + (sinc_factor*ratio*vstar2)^2 + $
       2d*sinc_factor*ratio*vstar1*vstar2*COS(2d*pi*s*proj_sep))/((1d)+ratio)^2

chi2 = TOTAL((((v2 - v2binary_model)^2)/(v2_errSqrd))/num_datapts, 2, /DOUBLE)

; The use of # above VECTORIZES the k loop below
;FOR k = 0, n_total-1 DO BEGIN
;    IF k MOD 400000 EQ 0 THEN PRINT,'k = ',k
;    ;proj_sep = xySep[k]*COS(pi*(xyPA[k]-basePA)/180d)
;    proj_sep = proj_sep_big[k,0:30]
;    sinc_arg = pi*chanw*uvrad*proj_sep/lam_eff
;    sinc_factor = sin(sinc_arg)/sinc_arg
;    v2binary_model = $
;      (vstar1^2 + (sinc_factor*ratio*vstar2)^2 + $ 
;                   2d*sinc_factor*ratio*vstar1*vstar2*COS(2d*pi*s*proj_sep))/((1d)+ratio)^2
;    chi2[k] = TOTAL((((v2 - v2binary_model)^2)/(v2_err^2))/num_datapts, /DOUBLE)
;    PRINT, 'Gridding ' + STRING(FORMAT='(F5.1)', FLOAT(k)/FLOAT(n_total-1)*100.0) + ' complete'
;END

;

index_min = WHERE(chi2 EQ MIN(chi2))
PRINT,'*************************************'
PRINT,'*****Binary grid search for '+starname
PRINT,'*****For a delta-mag of '+STRCOMPRESS(STRING(deltamag))
PRINT,'*************************************'
atx = STRTRIM(STRING(x[index_min]),1)
aty = STRTRIM(STRING(y[index_min]),1)
printMin = STRTRIM(STRING(MIN(chi2)),1)
PRINT, 'Min reduced chi^2 of ' + printMin +' found at x=' + atx +' y=' + aty
PRINT, 'This corresponds to Sep =' + STRING(SQRT(x[index_min]^2 + y[index_min]^2)) + $
  ' and PA of ' + STRING(90.0d - 180.0d/pi*ATAN(y[index_min],x[index_min]))

; Set some variables used for printing strings, and via the 
; COMMON statement returning variable values to cubegrid procedure.

subname=STRING(deltamag,FORMAT='(F3.1)')
theMin=MIN(chi2)
StrtheMin=STRCOMPRESS(theMin)
theMax=MAX(chi2)
StrtheMax=STRCOMPRESS(theMax)
rho_0=SQRT(x[index_min[0]]^2 + y[index_min[0]]^2)
pa_0=90.0d - 180.0d/pi*ATAN(y[index_min[0]],x[index_min[0]])
sep=rho_0
pa=pa_0
IF makedMagPlots EQ 'yes' THEN BEGIN
; set-up variables for determining sep and pa and for plotting
    plot_x = (DINDGEN(n_x)/DOUBLE(n_x-1))*(x_max-x_min)+x_min
    plot_y = (DINDGEN(n_y)/DOUBLE(n_y-1))*(y_max-y_min)+y_min
    plot_chi2 = DBLARR(n_x, n_y)
    counter = 0L
    FOR i = 0, n_x-1 DO BEGIN
        FOR j = 0, n_y-1 DO BEGIN
            plot_chi2[i,j] = chi2[counter]
            counter += 1
        END
    END

    IF ~psOutput THEN WINDOW, 1, XSIZE=700, YSIZE=700, TITLE='Chi^2 Space'
    IF psOutput THEN SET_PLOT, 'ps'
    IF psOutput THEN DEVICE, /COLOR, BITS_PER_PIXEL=8
    IF psOutput THEN DEVICE, FILENAME=datafile+'.'+subname+'.ps', ENCAPSULATED=0
    IF psOutput THEN DEVICE, XSIZE=7, YSIZE=7, XOFFSET=0.5, YOFFSET=3, /INCHES
    IF psOutput THEN pick_thick=4 ELSE pick_thick=1

    CONTOUR, plot_chi2, plot_x, plot_y, XRANGE=[MAX(plot_x),MIN(plot_x)], NLEVELS=10, $
      TITLE=file, SUBTITLE='Min '+sym_ChiSq+' of '+STRMID(StrtheMin,0,5)+' found at sep of '$ 
      +STRMID(STRTRIM(rho_0,2),0,5) + $
      ' and PA of ' + STRMID(STRTRIM(pa_0,2),0,6) + $
      '!C Maximum '+sym_ChiSq+ ' = '+STRMID(StrtheMax,0,5) + $
      '!C !C '+sym_Delta+'mag ='+STRING(deltamag,FORMAT='(F6.2)')+ $ 
      ', D1 ='+STRING(d1,FORMAT='(F6.2)') + $
      ' mas, D2 ='+STRING(d2,FORMAT='(F6.2)')+' mas.', $  
      XSTYLE=1, YSTYLE=1, C_ANNOTATION=''
    OPLOT, x[index_min], y[index_min], PSYM=1, THICK=pick_thick
    IF psOutput THEN DEVICE, /CLOSE_FILE

    SET_PLOT, 'X'

ENDIF

END

;********************************************************************
;+
; :Author: Bob Zavala
;-
; :Description
;
; Event handler for cubegridm_vec input widget
; Also includes event handler for initializing
; all the GRID PARAMS and the PLOTS plotting parameters.
; until a stand-alone event handler is made for this
; event.

PRO tab_grid_widget_event, ev

  COMMON GRIDPARAMS

  COMPILE_OPT hidden

  ; Retrieve the anonymous structure contained in the user value of
  ; the top-level base widget.
  WIDGET_CONTROL, ev.TOP, GET_UVALUE=stash

;  print,' '
;  print,' stash returns: ', stash
;  print,' '

  ;stop

  ; Check that COMMON gets the correct values returned from
  ; ww_init_grid_params after values are changed in the widget

  ;  print,'The COMMON variables before WIDGET_CONTROL are:'
  ;  print,''
  ;  print,'xmax = ',xmax
  ;  print,'xmin = ',xmin
  ;  print,'ymax = ',ymax
  ;  print,'ymin = ',ymin
  ;  print,'pixSize = ',pixSize
  ;  print,'numdmag = ',numdmag
  ;  print,'makeplots_q is: ',makeplots_q,'.'
  ;  print,'And makedMagPlots_q is: ',makedMagPlots_q,'.'
  ;  print,' '

  ; Use CASE statements to watch and act on button clicks
  ; I need to allow for clicks on the tabs GRID PARAMS and
  ; PLOTS and that their event ID is wTab

  ;  print,' '
  ;  print,' ev.ID is: ',ev.ID
  ;  print,' '

  CASE ev.ID OF

    stash.wTab:

    stash.wT1:

    stash.wT2:

    stash.wXmax: BEGIN
      WIDGET_CONTROL, ev.ID, get_value=xmax
      print,'Xmax = '+xmax[0]
      print,' '
      xmax = DOUBLE(xmax[0])
    END

    stash.wXmin: BEGIN
      WIDGET_CONTROL, ev.ID, get_value=xmin
      print,'Xmin = '+xmin[0]
      print,' '
      xmin = DOUBLE(xmin[0])
    END

    stash.wYmax: BEGIN
      WIDGET_CONTROL, ev.ID, get_value=ymax
      print,'Ymax = '+ymax[0]
      print,' '
      ymax = DOUBLE(ymax[0])
    END

    stash.wYmin: BEGIN
      WIDGET_CONTROL, ev.ID, get_value=ymin
      print,'Ymin = '+ymin[0]
      print,' '
      ymin = DOUBLE(ymin[0])
    END

    stash.wPixSize: BEGIN
      WIDGET_CONTROL, stash.wPixSize, get_value=PixSize
      print,'Pixel Size =  '+PixSize[0]
      print,' '
      pixSize = DOUBLE(PixSize[0])
    END

    stash.wDmag: BEGIN
      WIDGET_CONTROL, stash.wDmag, get_value=Dmag
      print,'numdmag set to: '+Dmag[0]
      print,' '
      numdmag = FIX(Dmag[0])
    END

    stash.wBgroup1: BEGIN
      WIDGET_CONTROL, ev.ID, get_value=do_sum_plots
      ;print,'do_sum_plots = ',do_sum_plots
      if do_sum_plots then makeplots_q = 'no' $
      else makeplots_q = 'yes'
      print,'Save Summary Plots = '+makeplots_q
      print,' '
    END

    stash.wBgroup2: BEGIN
      WIDGET_CONTROL, ev.ID, get_value=do_dmag_plots
      if do_dmag_plots then makedMagPlots_q = 'no' $
      else makedMagPlots_q = 'yes'
      print,'Save Delta-Mag Plots = '+makedMagPlots_q
      print,' '      
    END

    ; stash.bDone: WIDGET_CONTROL, ev.TOP, /DESTROY

  ENDCASE

  ; Check that COMMON gets the correct values returned from
  ; ww_init_grid_params after values are changed in the widget

  ;  print,' '
  ;  print,'The COMMON variables after WIDGET_CONTROL are:'
  ;  print,' '
  ;  print,'xmax = ',xmax
  ;  print,'xmin = ',xmin
  ;  print,'ymax = ',ymax
  ;  print,'ymin = ',ymin
  ;  print,'pixSize = ',pixSize
  ;  print,'numdmag = ',numdmag
  ;  print,'makeplots_q is: ',makeplots_q,'.'
  ;  print,'And makedMagPlots_q is: ',makedMagPlots_q,'.'
  ;  print,' '

  ; Use CASE statements to watch and act on button clicks
  ; I need to allow for clicks on the tabs GRID PARAMS and
  ; PLOTS and that their event ID is wTab

  ;  print,' '
  ;  print,' ev.ID is: ',ev.ID
  ;  print,' '


END

;*********************************************************************

; :Description
; Procedure to initialize the gridding and plotting parameters

PRO ww_init_grid_params, ev, stash

  COMMON GRIDPARAMS

  print,' '
  ;print,'ww_init_grid_params is called with event: ',STRTRIM(STRING(ev),1)
  print,'Gridfit parameters initialized as: '
  print,' '

  ; Set the gridfit parameters
  ; Note that these are returned from the widget as a
  ; string and are doubles for fitting and need conversion
  
  WIDGET_CONTROL, stash.wXmax, get_value=xmax
  xmax = DOUBLE(xmax[0])
  print,'  xmax = ',xmax
  WIDGET_CONTROL, stash.wXmin, get_value=xmin
  xmin = DOUBLE(xmin[0])
  print,'  xmin = ',xmin
  WIDGET_CONTROL, stash.wYmax, get_value=ymax
  ymax = DOUBLE(ymax[0])
  print,'  ymax = ',ymax
  WIDGET_CONTROL, stash.wYmin, get_value=ymin
  ymin = DOUBLE(ymin[0])
  print,'  ymin = ',ymin
  WIDGET_CONTROL, stash.wPixSize, get_value=PixSize
  pixSize = DOUBLE(PixSize[0])
  print,'  pixSize = ',pixSize
  WIDGET_CONTROL, stash.wDmag, get_value=Dmag
  numdmag = FIX(Dmag[0])
  print,'  numdmag = ',numdmag


  ;print,' '
  ; Set the plotting parameters
  WIDGET_CONTROL, stash.wBgroup1, get_value=do_sum_plots
  if do_sum_plots then makeplots_q = 'no' $
  else makeplots_q = 'yes'
  WIDGET_CONTROL, stash.wBgroup2, get_value=do_dmag_plots
  if do_dmag_plots then makedMagPlots_q = 'no' $
  else makedMagPlots_q = 'yes'
  print,' '
  print,' Save Summary Plots   = ',makeplots_q
  print,' Save Delta-Mag Plots = ',makedMagPlots_q
  print,' '

END

;*********************************************************************

; :Description
; Widget creation routine for cubegridm_vec input widget
; Note the default font found via inserting these two lines after wT1 = line:
; wTabFont = WIDGET_INFO(wTab, /FONTNAME)
; print,'wTabFont: ',wTabFont
; was:
; wTabFont: -Misc-Fixed-Medium-R-SemiCondensed--13-120-75-75-C-60-ISO8859-1


PRO tab_grid_widget, LOCATION=location

  ; Load available device fonts into a string array simplifying font selection

  device, get_fontnames=fontes, set_font='*'

  ; Create the top-level base (wTLB) and the tab.
  wTLB = WIDGET_BASE(xsize=310, ysize=310, title='GRIDFIT INPUTS', /COLUMN, /BASE_ALIGN_TOP)
  wTab = WIDGET_TAB(wTLB, xsize=290, ysize=240, LOCATION=location)

  ;; Create the first widget tab (wT1) base to hold grid input parameters
  wT1 = WIDGET_BASE(wTab, TITLE='GRID PARAMS', COLUMN = 2)
  wXmaxLabel = WIDGET_LABEL(wT1, VALUE='Enter Xmax (mas)', FONT=fontes[210])
  wXmax= WIDGET_TEXT(wT1, VALUE='+20.0', /EDITABLE);, /ALL_EVENTS comment out /ALL_EVENTS
  ;need Enter to register an event)
  wXminLabel = WIDGET_LABEL(wT1, VALUE='Enter Xmin (mas)', FONT=fontes[210])
  wXmin= WIDGET_TEXT(wT1, VALUE='-20.0', /EDITABLE)
  wPixSizeLabel = WIDGET_LABEL(wT1, VALUE='Pixel size (mas)', FONT=fontes[210])
  wPixSize = WIDGET_TEXT(wT1, VALUE='0.2', /EDITABLE)

  wYmaxLabel = WIDGET_LABEL(wT1, VALUE='Enter Ymax (mas)', FONT=fontes[210])
  wYmax= WIDGET_TEXT(wT1, VALUE='+20.0', /EDITABLE)
  wYminLabel = WIDGET_LABEL(wT1, VALUE='Enter Ymin (mas)', FONT=fontes[210])
  wYmin= WIDGET_TEXT(wT1, VALUE='-20.0', /EDITABLE)
  wDmagLabel = WIDGET_LABEL(wT1, VALUE='Delta mag. (mag)', FONT=fontes[210])
  wDmag = WIDGET_TEXT(wT1, VALUE='5', /EDITABLE)

  ;; Create the 2nd widget tab (wT2) to hold plotting buttons
  wT2 = WIDGET_BASE(wTab, TITLE='PLOTS', /COLUMN)
  ; Set default (SET_VALUE) to save summary plots
  wLabel1 = WIDGET_LABEL(wT2, VALUE='Save Summary Plots?', FONT=fontes[227])
  wBgroup1 = CW_BGROUP(wT2, ['Yes', 'No'], FONT=fontes[227], $
    /ROW, /EXCLUSIVE, SET_VALUE=0, /RETURN_NAME)
  ; Set default (SET_VALUE) not to save summary plots
  wLabel2 = WIDGET_LABEL(wT2, VALUE='Save Delta-Mag Plots?', FONT=fontes[227])
  wBgroup2 = CW_BGROUP(wT2, ['Yes', 'No'], FONT=fontes[227], $
    /ROW, /EXCLUSIVE, SET_VALUE=1, /RETURN_NAME)

  ; Create a base widget to hold the 'GRIDFIT' button, and
  ; the button itself.
  wControl = WIDGET_BASE(wTLB, /ROW)
  bGridFit = WIDGET_BUTTON(wControl, VALUE='GRIDFIT', FONT=fontes[138], EVENT_PRO='cubegridm_vec')
  
  ; Add to the base widget wControl a button to call the widget for 
  ; input data file selection. Incorporate this into a controlling 
  ; procedure for the buttons
  bDatFile = WIDGET_BUTTON(wControl, VALUE='Datfile', FONT=fontes[138], EVENT_FUNC='ww_buttoncall')
  
  ; Add a Help button
  bHelp = WIDGET_BUTTON(wControl, VALUE='Help', FONT=fontes[138], EVENT_FUNC='ww_buttoncall')
  
  ; Add a Done button
  bDone = WIDGET_BUTTON(wControl, VALUE='Done', FONT=fontes[138], EVENT_FUNC='ww_buttoncall')

  ; Create an anonymous structure to hold widget IDs. This
  ; structure becomes the user value of the top-level base
  ; widget.
  stash = { wTab:wTab, wT1:wT1, wT2:wT2, wXmax:wXmax, wXmin:wXmin, wYmax:wYmax, wYmin:wYmin, $
    wPixSize:wPixSize, wDmag:wDmag, wBgroup1:wBgroup1, wBgroup2:wBgroup2, bGridFit:bGridFit, $
    bDatFile:bDatFile, bHelp:bHelp, bDone:bDone }

  ; Realize the widgets, set the user value of the top-level
  ; base to the structure stash, and call XMANAGER to manage everything.
  WIDGET_CONTROL, wTLB, /REALIZE
  WIDGET_CONTROL, wTLB, SET_UVALUE=stash
  XMANAGER, 'tab_grid_widget', wTLB, /NO_BLOCK

  ; Call the parameter initialization procedure in order to
  ; set the defaults, using the event id for the top level base
  ; wTLB as that is logically the one to use here. I tested and
  ; the id's for the tabs GRID PARAMS and PLOTS also work.
  ww_init_grid_params, wTLB, stash


END
