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

  COMMON SHAREIT, x_min, x_max, y_min, y_max, pixelSize, $
                  num_dmag, makeplots, makedMagPlots

  COMPILE_OPT hidden
  
  ; Retrieve the anonymous structure contained in the user value of
  ; the top-level base widget.
  WIDGET_CONTROL, ev.TOP, GET_UVALUE=stash
  
;  print,' '
;  print,' stash returns: ', stash
;  print,' '
    
  ; Check that COMMON gets the correct values returned from 
  ; ww_init_grid_params after values are changed in the widget
  
;  print,'The COMMON variables before WIDGET_CONTROL are:'
;  print,''
;  print,'x_max = ',x_max
;  print,'x_min = ',x_min
;  print,'y_max = ',y_max  
;  print,'y_min = ',y_min  
;  print,'pixelSize = ',pixelSize
;  print,'num_dmag = ',num_dmag
;  print,'makeplots is: ',makeplots,'.'
;  print,'And makedMagPlots is: ',makedMagPlots,'.'
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
        x_max = DOUBLE(xmax)
      END
      
    stash.wXmin: BEGIN
        WIDGET_CONTROL, ev.ID, get_value=xmin
        x_min = DOUBLE(xmin)
      END
      
    stash.wYmax: BEGIN
        WIDGET_CONTROL, ev.ID, get_value=ymax
        y_max = DOUBLE(ymax)
      END
      
    stash.wYmin: BEGIN
        WIDGET_CONTROL, ev.ID, get_value=ymin
        y_min = DOUBLE(ymin)
      END
    
    stash.wPixSize: BEGIN
       WIDGET_CONTROL, stash.wPixSize, get_value=PixSize
       pixelSize = DOUBLE(PixSize)       
      END
      
    stash.wDmag: BEGIN
        WIDGET_CONTROL, stash.wDmag, get_value=Dmag
        num_dmag = FIX(Dmag)
      END
        
    stash.wBgroup1: BEGIN
        WIDGET_CONTROL, ev.ID, get_value=do_sum_plots
        print,'do_sum_plots = ',do_sum_plots
        if do_sum_plots then makeplots = 'no' $
        else makeplots = 'yes' 
      END
    
    stash.wBgroup2: BEGIN 
        WIDGET_CONTROL, ev.ID, get_value=do_dmag_plots
        if do_dmag_plots then makedMagPlots = 'no' $
        else makedMagPlots = 'yes'        
      END      
    
    ; User is DONE
    stash.bDone: WIDGET_CONTROL, ev.TOP, /DESTROY
    
  ENDCASE
  
  ; Check that COMMON gets the correct values returned from
  ; ww_init_grid_params after values are changed in the widget

  print,' '
  print,'The COMMON variables after WIDGET_CONTROL are:'
  print,' '
  print,'x_max = ',x_max
  print,'x_min = ',x_min
  print,'y_max = ',y_max
  print,'y_min = ',y_min
  print,'pixelSize = ',pixelSize
  print,'num_dmag = ',num_dmag
  print,'makeplots is: ',makeplots,'.'
  print,'And makedMagPlots is: ',makedMagPlots,'.'
  print,' '

  ; Use CASE statements to watch and act on button clicks
  ; I need to allow for clicks on the tabs GRID PARAMS and
  ; PLOTS and that their event ID is wTab

  print,' '
  print,' ev.ID is: ',ev.ID
  print,' '
 
  
END

;*********************************************************************

; :Description
; Procedure to initialize the gridding and plotting parameters

PRO ww_init_grid_params, ev, stash

  COMMON SHAREIT, x_min, x_max, y_min, y_max, pixelSize, $
    num_dmag, makeplots, makedMagPlots

  print,' '
  print,'ww_init_grid_params is called with event: ',STRTRIM(STRING(ev),1)
  print,' '
  ; Set the gridfit parameters
  ; Note that these are returned from the widget as a 
  ; string and are doubles for fitting and need conversion
  WIDGET_CONTROL, stash.wXmax, get_value=xmax
  x_max = DOUBLE(xmax)
  print,'  x_max = ',x_max
  WIDGET_CONTROL, stash.wXmin, get_value=xmin
  x_min = DOUBLE(xmin)
  print,'  x_min = ',x_min
  WIDGET_CONTROL, stash.wYmax, get_value=ymax
  y_max = DOUBLE(ymax)
  print,'  y_max = ',y_max
  WIDGET_CONTROL, stash.wYmin, get_value=ymin
  y_min = DOUBLE(ymin)
  print,'  y_min = ',y_min
  WIDGET_CONTROL, stash.wPixSize, get_value=PixSize
  pixelSize = DOUBLE(PixSize)
  print,'  pixelSize = ',pixelSize
  WIDGET_CONTROL, stash.wDmag, get_value=Dmag
  num_dmag = FIX(Dmag)
  print,'  num_dmag = ',num_dmag
  
  
  print,' ' 
  ; Set the plotting parameters
  WIDGET_CONTROL, stash.wBgroup1, get_value=do_sum_plots
  if do_sum_plots then makeplots = 'no' $
  else makeplots = 'yes'
  WIDGET_CONTROL, stash.wBgroup2, get_value=do_dmag_plots
  if do_dmag_plots then makedMagPlots = 'no' $
  else makedMagPlots = 'yes'
  print,' '
  print,'In ww_init_grid_params makeplots = ',makeplots
  print,'In ww_init_grid_params makedMagPlots = ',makedMagPlots

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
    
  ; Create a base widget to hold the 'Done' button, and
  ; the button itself.
  wControl = WIDGET_BASE(wTLB, /ROW)
  bDone = WIDGET_BUTTON(wControl, VALUE='Done', FONT=fontes[138])
  
  ; Create an anonymous structure to hold widget IDs. This
  ; structure becomes the user value of the top-level base
  ; widget.
  stash = { wTab:wTab, wT1:wT1, wT2:wT2, wXmax:wXmax, wXmin:wXmin, wYmax:wYmax, wYmin:wYmin, $ 
          wPixSize:wPixSize, wDmag:wDmag, wBgroup1:wBgroup1, wBgroup2:wBgroup2, bDone:bDone }

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
