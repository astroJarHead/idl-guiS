; Development Event handler for cubegridm_vec.pro GUI
; 

PRO gpar_widget_event, ev

  Dmag = ''

  WIDGET_CONTROL, ev.TOP, GET_UVALUE=textwid ; textwid = ev.ID of the text box above 'Done' button
                                             ; They increment top to bottom by one
  WIDGET_CONTROL, ev.ID, GET_UVALUE=uval

  CASE uval OF

    'ONE' : WIDGET_CONTROL, textwid, SET_VALUE='Button 1 Pressed '+STRTRIM(STRING(ev.ID),1)+' '+STRTRIM(STRING(textwid),1)

    'TWO' : WIDGET_CONTROL, textwid, SET_VALUE='Button 2 Pressed '+STRTRIM(STRING(ev.ID),1)+' '+STRTRIM(STRING(textwid),1)

    'DONE': WIDGET_CONTROL, ev.TOP, /DESTROY
    
    ELSE: BEGIN 
      WIDGET_CONTROL, textwid, SET_VALUE='Delta mag set '+STRTRIM(STRING(ev.ID),1)
      ; the next line uses the event ID of the text window where the Delta mag
      ; is entered
      WIDGET_CONTROL, ev.ID, GET_VALUE=Dmag
      print,'Dmag recorded as: ',Dmag ; Here's how I can use the Dmag and other variables. 
    END
  ENDCASE

END



PRO gpar_widget
; call this, and the widget appears

common GridParam,Dmag

  ; Initialize delta-mag
  Dmag = 4.0

  base = WIDGET_BASE(/COLUMN,XSIZE=200)

  button1 = WIDGET_BUTTON(base, VALUE='One', UVALUE='ONE')

  button2 = WIDGET_BUTTON(base, VALUE='Two', UVALUE='TWO')
  
  wLabel = WIDGET_LABEL(base, VALUE='Enter Delta mag.')
  
  wText= WIDGET_TEXT(base, /EDITABLE, value=string(Dmag,format='(f4.1)'), $  
         UVALUE='DMAG') ; generate event on EndOfLine

  text = WIDGET_TEXT(base, XSIZE=20)

  button3 = WIDGET_BUTTON(base, value='Done', UVALUE='DONE')

  WIDGET_CONTROL, base, SET_UVALUE=text ; sets textwid in WIDGET_CONTROL, ev.TOP in event handler

  WIDGET_CONTROL, base, /REALIZE

  XMANAGER, 'gpar_widget', base

END

