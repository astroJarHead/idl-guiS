pro gatorplot_event,ev

  widget_control,ev.id,get_uvalue=uval
  widget_control,ev.top,get_uvalue=t
  
  CASE uval OF
     
    'filelist' : BEGIN
      longName=dialog_pickfile(title='Select one (1) IDL PRO file',filter='*.pro',GET_PATH=thePath)
      afile=STRMID(longName,STRLEN(longName)-STRLEN(thePath)-1,/REVERSE_OFFSET)
      WIDGET_CONTROL, ev.id, set_value=afile
    END
    
    'xlinlog' : 
    
    'titletext' : 
    
    'export' : BEGIN
      WIDGET_CONTROL,t[0],get_value=exportProFile
      WIDGET_CONTROL,t[1],get_value=outFileName
      if outFileName EQ '' THEN BEGIN
        error1=DIALOG_MESSAGE('Oops, you did not provide an Export filename.',/ERROR)
      endif else begin 
        OPENW,LUN,outFileName,/GET_LUN
        PRINTF,LUN,exportProFile
        FREE_LUN, LUN
        aResponse=DIALOG_MESSAGE('PRO file name saved to: '+outFileName,/INFORMATION)
      endelse  
    END
    
    'draw' : SHADE_SURF, DIST(150)
    
    'DONE': WIDGET_CONTROL, ev.TOP, /DESTROY
    
  ENDCASE

end

pro gatorplot

; Inspired by Craig Warner at the University of Florida:
; https://users.astro.ufl.edu/~warner/IDL5220/week6.html

; Call this procedure to create the widget. 

  ; create base widget
  base=widget_base(xsize=800,ysize=640,title='GatorPlot v3.03  by: Bob Zavala')
  ; a title display
  labeltitle=widget_label(base,value='Export Title:',yoffset=162)
  ; an export button
  export=widget_button(base,value='Export',yoffset=136,xoffset=40,uvalue='export')
  ; Two text widgets for lists of files and a title
  filelist=widget_text(base,xsize=18,ysize=5,/scroll,/editable,/all_events, $ 
  yoffset=28,uvalue='filelist')
  titletext=widget_text(base,xsize=10,/editable,/all_events,yoffset=180, $
  uvalue='titletext')
  ; Now let us create a draw_widget to draw a plot
  draw=widget_draw(base,xsize=640,ysize=512,xoffset=160,uvalue='draw',/button_events)
  ; a compound widget example
  xlinlog=cw_bgroup(base,['Linear','Log'],yoffset=580,xoffset=300,/row,/frame, $ 
  /exclusive, uvalue='xlinlog',label_top='X-axis:',set_value=0)
  ; Add a 'Done' button
  doneButton = WIDGET_BUTTON(base, value='Done', UVALUE='DONE', xoffset=40, $ 
  yoffset=220)
 
  ; Set uvalue of base to an array containing the widget id's of many 
  ; widgets within or under base. 
  widget_control,base,set_uvalue=[filelist,titletext,draw,xlinlog]
  
  ; Create the widget!
  WIDGET_CONTROL, base, /REALIZE
  ; Use XMANAGER to register widget and handle events
  XMANAGER, 'gatorplot', base

end