pro best_fit_1_11, config_file
; BESTFIT version 1.11 by Ivano Baronchelli
; march 2014

;changes:
; 1.0  -> first version
; 1.10 -> extinction is now considered
; 1.11 -> models are read only one only time during the hypercube
;         building.

version=strcompress(string(1))+"."+strcompress(string(11))
print, " "
print, "**********************************************"
print, "**********   BESTFIT version ",version,"   ***********" 
print, "**********        March 2014       ***********" 
print, "***********  by Ivano Baronchelli ************" 
print, "**********************************************"
print, " "

; CONFIGURATION FILE
readcol, config_file, PARAM_KEYW,PARAM,format='a,a',/silent

print, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

;GET_FLUX
IDX_get_flux=where(PARAM_KEYW eq 'GET_FLUX')
GET_FLUX=PARAM[IDX_get_flux[0]]

;WAV_OR_FILTER
IDX_wav_or_filter=where(PARAM_KEYW eq 'WAV_OR_FILTER')
WAV_OR_FILTER=PARAM[IDX_wav_or_filter[0]]

;FILTER_GET_FLUX
IDX_filter_get_flux=where(PARAM_KEYW eq 'FILTER_GET_FLUX')
FILTER_GET_FLUX=PARAM[IDX_filter_get_flux[0]]

;FILTER_CONVOL
IDX_filters_convol=where(PARAM_KEYW eq 'FILTER_CONVOL')
FILTER_CONVOL=PARAM[IDX_filters_convol[0]]
if FILTER_CONVOL eq 'no' then print, "Central wavelengt of the filters considered"
if FILTER_CONVOL eq 'yes' or (GET_FLUX eq 'yes' and WAV_OR_FILTER eq 'filter') then begin
print, "Convolution with filter response considered"
;FILTERS_RESPONSES_FILE
IDX_filters_responses_file=where(PARAM_KEYW eq 'FILTERS_RESPONSES_FILE')
FILTERS_RESPONSES_FILE=PARAM[IDX_filters_responses_file[0]]
print, "Filters responses from ",FILTERS_RESPONSES_FILE
endif

; FILTERS FILE
; FILTER_NAME
; FILTER_WLGHT
IDX_filters_file=where(PARAM_KEYW eq 'FILTERS_FILE')
FILTERS_FILE=PARAM[IDX_filters_file[0]]
if FILTER_CONVOL eq 'no' then begin
readcol,FILTERS_FILE,FILTER_WLGHT, FILTER_NAME,format='d,a',STRINGSKIP='#',/silent
N_FILTERS=strcompress(n_elements(FILTER_NAME))
minL=min(FILTER_WLGHT)/2.5 ;for visualization
maxL=max(FILTER_WLGHT)*2.5 ;for visualization
endif
if FILTER_CONVOL eq 'yes' then begin
readcol,FILTERS_FILE,FILTER_LOG,format='l',STRINGSKIP='#',/silent
N_FILTERS=strcompress(n_elements(FILTER_LOG))
endif
print,"Number of filters:",strcompress(N_FILTERS)

;INPUT_UNITS
IDX_input_units=where(PARAM_KEYW eq 'INPUT_UNITS')
INPUT_UNITS=PARAM[IDX_input_units[0]]

; OUTPUT_PATH
IDX_output_path=where(PARAM_KEYW eq 'OUTPUT_PATH')
OUTPUT_PATH=PARAM[IDX_output_path[0]]
if OUTPUT_PATH[0] eq '#' then OUTPUT_PATH=''

; OUT_BEST_FIT_SED_TXT (yes/no)
IDX_out_best_fit_sed_txt=where(PARAM_KEYW eq 'OUT_BEST_FIT_SED_TXT')
OUT_BEST_FIT_SED_TXT=PARAM[IDX_out_best_fit_sed_txt[0]]

; OUT_BEST_FIT_SED_PS (yes/no)
IDX_out_best_fit_sed_ps=where(PARAM_KEYW eq 'OUT_BEST_FIT_SED_PS')
OUT_BEST_FIT_SED_PS=PARAM[IDX_out_best_fit_sed_ps[0]]

;SED_TYPE
IDX_sed_type=where(PARAM_KEYW eq 'SED_TYPE')
SED_TYPE=PARAM[IDX_sed_type[0]]

; OUTPUT FILE NAME 1
IDX_output_file1=where(PARAM_KEYW eq 'OUTPUT_FILE1')
OUTPUT_FILE1=OUTPUT_PATH+PARAM[IDX_output_file1[0]]

; OUTPUT FILE NAME 2
IDX_output_file2=where(PARAM_KEYW eq 'OUTPUT_FILE2')
OUTPUT_FILE2=OUTPUT_PATH+PARAM[IDX_output_file2[0]]


; EXTINCTION
IDX_extinct=where(PARAM_KEYW eq 'EXTINCT')
EXTINCT=PARAM[IDX_extinct[0]]
; AV_MIN
IDX_Avmin=where(PARAM_KEYW eq 'AV_MIN')
AV_MIN=double(PARAM[IDX_Avmin[0]])
; AV_MAX
IDX_Avmax=where(PARAM_KEYW eq 'AV_MAX')
AV_MAX=double(PARAM[IDX_Avmax[0]])
; AV_STEP
IDX_Av_step=where(PARAM_KEYW eq 'AV_STEP')
AV_STEP=double(PARAM[IDX_Av_step[0]])

IF EXTINCT ne 'yes' then begin
AV_MIN=0.
AV_MAX=0.
AV_STEP=1.
endif

;KUR93_PAR
IDX_KUR93=where(PARAM_KEYW eq 'KUR93')
KUR93=PARAM[IDX_KUR93[0]]
if KUR93 eq 'yes' then print,'Using Kurucz 93 Models'
if KUR93 eq 'no' then print,'Using generic templated models'
if KUR93 ne 'yes' and KUR93 ne 'no' then begin
print,  "wrong KUR93 parameter"
stop
endif 

; (METALLICITY) MODELS
; Models of metallicities if Kurucz models are used.
; SED models if generic models are used
IDX_metallicity_models=where(PARAM_KEYW eq 'METAL_MODELS_FILE')
METAL_MODELS_FILE=PARAM[IDX_metallicity_models[0]]
readcol,METAL_MODELS_FILE,METAL_MODEL,format='a',STRINGSKIP='#',/silent
str0=''
if KUR93 eq 'yes' then str0=' metallicity'
print,"Number of"+str0+" models:",strcompress(n_elements(METAL_MODEL))

; VERBOSE
; The  analized model is printed, togheter with the corresponding chi2 and the best fit chi2
IDX_verbose=where(PARAM_KEYW eq 'VERBOSE')
VERBOSE=PARAM[IDX_verbose[0]]

;VISUALIZE_FIT
IDX_visualize_fit=where(PARAM_KEYW eq 'VISUALIZE_FIT')
VISUALIZE_FIT=PARAM[IDX_visualize_fit[0]]

; REFEREMENT FILTER
; The fluxes will be normalized to the flux in this filter 
IDX_ref_filter=where(PARAM_KEYW eq 'REF_FILTER')
val=PARAM[IDX_ref_filter[0]] ; OPTIONS: -1 or positive integer number
if val ne -1 then begin
REFEREMENT_FILTER=PARAM[IDX_ref_filter[0]]-1 
if FILTER_CONVOL eq 'no' then print, "Referement filter: ",FILTER_NAME[REFEREMENT_FILTER]
if FILTER_CONVOL eq 'yes' then print, "Referement filter: ",REFEREMENT_FILTER+1
endif
if val eq -1 then  REFEREMENT_FILTER=PARAM[IDX_ref_filter[0]]

;-------------------------------------------------------------------
; EXTINCTION
if EXTINCT eq 'yes' then begin
print,'Av min :',strcompress(string(AV_MIN))
print,'Av max :',strcompress(string(AV_MAX))
print,'Av step:',strcompress(string(AV_STEP))
endif

; ONLY FOR KUR 93 MODELS
if KUR93 eq 'yes' then begin

;GRAVITY_MIN - only Kuruck models
IDX_gravity_min=where(PARAM_KEYW eq 'GRAVITY_MIN')
GRAVITY_MIN=float(PARAM[IDX_gravity_min[0]])
print,"Gravity (log_g relative to solar) min:",strcompress(GRAVITY_MIN)

;GRAVITY_MAX - only Kuruck models
IDX_gravity_max=where(PARAM_KEYW eq 'GRAVITY_MAX')
GRAVITY_MAX=float(PARAM[IDX_gravity_max[0]])
print,"Gravity (log_g relative to solar) max:",strcompress(GRAVITY_MAX)

; T_MIN - only Kuruck models
IDX_T_min=where(PARAM_KEYW eq 'T_MIN')
T_MIN=float(PARAM[IDX_T_min[0]])
print,"T min:",string(T_MIN)+" K"

; T_MAX - only Kuruck models
IDX_T_max=where(PARAM_KEYW eq 'T_MAX')
T_MAX=float(PARAM[IDX_T_max[0]])
print,"T max:",string(T_MAX)+" K"

; TEMPERATURE MODELS - only Kuruck models
T_available=[3500.0 , 3750.0 , 4000.0 , 4250.0 , 4500.0 , 4750.0 , 5000.0 , 5250.0 , 5500.0 , 5750.0 , 6000.0 , 6250.0 , 6500.0 , 6750.0 , 7000.0 , 7250.0 , 7500.0 , 7750.0 , 8000.0 , 9250.0 , 9500.0 , 9750.0 , 10000.0 , 10500.0 , 11000.0 , 11500.0 , 12000.0 , 12500.0 , 13000.0 , 14000.0 , 15000.0 , 16000.0 , 17000.0 , 18000.0 , 19000.0 , 20000.0 , 21000.0 , 22000.0 , 23000.0 , 24000.0 , 25000.0 , 26000.0 , 27000.0 , 28000.0 , 29000.0 , 30000.0 , 31000.0 , 32000.0 , 33000.0 , 34000.0 , 35000.0 , 37500.0 , 40000.0 , 42500.0 , 45000.0 , 47500.0, 50000.0]
T_selected_idx=where(T_available ge T_MIN and T_available le T_MAX)
T_names=strcompress(string(long(T_available[T_selected_idx])),/remove_all)
print,"Number of temperature models (available, used):",strcompress(n_elements(T_available)),strcompress(n_elements(T_selected_idx))

; GRAVITY MODELS - only Kuruck models
g_available=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
g_selected_idx=where(g_available ge GRAVITY_MIN and g_available le GRAVITY_MAX)
g_names="G"+strcompress(string(long(10*g_available[g_selected_idx])),/remove_all)
lt1=where(g_available[g_selected_idx] lt 1.0)
if lt1[0] ne -1 then g_names[lt1]="G0"+strcompress(string(long(10*g_available[g_selected_idx[lt1]])),/remove_all)
print,"Number of gravity models (available, used):",strcompress(n_elements(g_available)),strcompress(n_elements(g_selected_idx))

endif;END - ONLY FOR KUR 93 MODELS - END


;GRAPHICAL VALUES GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
;PS_THICK
IDX_ps_thick=where(PARAM_KEYW eq 'PS_THICK')
PS_THICK=PARAM[IDX_ps_thick[0]]
;PS_CHSIZE
IDX_ps_chsize=where(PARAM_KEYW eq 'PS_CHSIZE')
PS_CHSIZE=PARAM[IDX_ps_chsize[0]]
;PS_XSIZE
IDX_ps_xsize=where(PARAM_KEYW eq 'PS_XSIZE')
PS_XSIZE=PARAM[IDX_ps_xsize[0]]
;PS_YSIZE
IDX_ps_ysize=where(PARAM_KEYW eq 'PS_YSIZE')
PS_YSIZE=PARAM[IDX_ps_ysize[0]]
;RIGHT_EXPANSION
IDX_right_expansion=where(PARAM_KEYW eq 'RIGHT_EXPANSION')
RIGHT_EXPANSION=PARAM[IDX_right_expansion[0]]
;LEFT_EXPANSION
IDX_left_expansion=where(PARAM_KEYW eq 'LEFT_EXPANSION')
LEFT_EXPANSION=PARAM[IDX_left_expansion[0]]
;UP_EXPANSION
IDX_up_expansion=where(PARAM_KEYW eq 'UP_EXPANSION')
UP_EXPANSION=PARAM[IDX_up_expansion[0]]
;LOW_EXPANSION
IDX_low_expansion=where(PARAM_KEYW eq 'LOW_EXPANSION')
LOW_EXPANSION=PARAM[IDX_low_expansion[0]]
;GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

; CATALOG NAME
IDX_catalog_name=where(PARAM_KEYW eq 'CATALOG_NAME')
CATALOG_NAME=PARAM[IDX_catalog_name[0]]
print,"Catalog name:",CATALOG_NAME
; ORIGINAL CATALOG READING
original_cat=read_ascii(CATALOG_NAME)
;first_stars_fluxes=cat.field01[*,0]

print, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
; END OF PARAMETER READING.






;OPERATIONS ON THE CATALOG:






; CATALOG MAGNITUDE->FLUX CONVERSION
; TRANSFORMATIONS here below:
;TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
;original_cat -> original_Jy_cat ->normalized_cat
;TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
IDX_input_units=where(PARAM_KEYW eq 'INPUT_UNITS')
INPUT_UNITS=PARAM[IDX_input_units[0]]
if INPUT_UNITS eq "Jy" or INPUT_UNITS eq "JY"  then original_Jy_cat=original_cat
if INPUT_UNITS eq "AB_MAG" then begin
original_Jy_cat=original_cat
original_Jy_cat.field01[1:N_FILTERS,*]=10^(-(original_cat.field01[1:N_FILTERS,*]-8.9)/2.5)
FLUX_SUP=10^(-(original_cat.field01[1:N_FILTERS,*]-original_cat.field01[N_FILTERS+1:N_FILTERS*2,*]-8.9)/2.5)
FLUX_INF=10^(-(original_cat.field01[1:N_FILTERS,*]+original_cat.field01[N_FILTERS+1:N_FILTERS*2,*]-8.9)/2.5)
FLUXERR=(FLUX_SUP-FLUX_INF)/2.
original_Jy_cat.field01[N_FILTERS+1:N_FILTERS*2,*]=FLUXERR
; The following rule is applied only in the case the input units are
; AB magnitudes. Otherwise, only the data with associated -99. errors
; are considered as not valid numbers 
nonumbers=where(original_cat.field01 lt -90. or original_cat.field01 gt 90. or original_Jy_cat.field01 ne original_Jy_cat.field01)
original_Jy_cat.field01[nonumbers]=-99.
endif

normalized_cat=original_Jy_cat

;---------------------------------------------------------
; Referement Source
FLUX_RF_SOURCE=dblarr(n_elements(original_Jy_cat.field01[0,*]))
FLUX_RF_SOURCE_test=dblarr(n_elements(original_Jy_cat.field01[0,*]))
;---------------------------------------------------------

; CATALOG NORMALIZATION (only if a referement filter is indicated)
ss=0l
while ss lt n_elements(original_Jy_cat.field01[0,*]) do begin
; IF A REFEREMENT FILTER IS INDICATED, WE NORMALIZE THE CATALOG HERE,
; OTHERWISE, THE NORMaLIZATION IS PERFORMED DURING THE COMPARISON WITH
; THE MODELS
if REFEREMENT_FILTER ne -1 then begin
FLUX_RF_SOURCE[ss]=original_Jy_cat.field01[REFEREMENT_FILTER+1,ss]
normalized_cat.field01[1:N_FILTERS*2,ss]=original_Jy_cat.field01[1:N_FILTERS*2,ss]/FLUX_RF_SOURCE[ss]
endif
original_ASS_ERR_tempor=original_Jy_cat.field01[N_FILTERS+1:N_FILTERS*2,ss]
normalized_ASS_ERR_tempor=normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,ss]
normalized_cat_tempor=normalized_cat.field01[1:N_FILTERS,ss]
idcg=where(original_ASS_ERR_tempor lt 0 or original_ASS_ERR_tempor ge 90.)
if idcg[0] ne -1 then begin
normalized_ASS_ERR_tempor[idcg]=-99.
normalized_cat_tempor[idcg]=-99.
endif
normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,ss]=normalized_ASS_ERR_tempor
normalized_cat.field01[1:N_FILTERS,ss]=normalized_cat_tempor
ss=ss+1
endwhile

; SOURCES CATALOG NORMALIZATION IF NO REFEREMENT FILTER IS INDICATED
if REFEREMENT_FILTER eq -1 then begin
SRC=0L
while SRC lt n_elements(normalized_cat.field01[0,*]) do begin
temp_err_cat=original_Jy_cat.field01[N_FILTERS+1:N_FILTERS*2,SRC]
temp_flx_cat=original_Jy_cat.field01[1:N_FILTERS,SRC]
ok_idxs=where(temp_err_cat gt 0.)
; SIMPLE MEAN ------------
;FLUX_RF_SOURCE[SRC]=mean(temp_flx_cat[ok_idxs])
; WEIGHTED MEAN!!! -------
FLUX_RF_SOURCE[SRC]=total(temp_flx_cat[ok_idxs]/temp_err_cat[ok_idxs])/total(1./temp_err_cat[ok_idxs])
;-------------------------
temp_flx_cat[ok_idxs]=temp_flx_cat[ok_idxs]/FLUX_RF_SOURCE[SRC]
temp_err_cat[ok_idxs]=temp_err_cat[ok_idxs]/FLUX_RF_SOURCE[SRC]
normalized_cat.field01[1:N_FILTERS,SRC]=temp_flx_cat
normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,SRC]=temp_err_cat
SRC=SRC+1
endwhile
endif

;RIPRISTINATE IDs ; ID Ripristinated as in original input catalog --
original_Jy_cat.field01[0,*]=original_cat.field01[0,*]
 normalized_cat.field01[0,*]=original_cat.field01[0,*]





; OPERATIONS ON THE FILTERS:





; READING FILTERS RESPONSE CURVES
if FILTER_CONVOL eq 'yes' or (GET_FLUX eq 'yes' and WAV_OR_FILTER eq 'filter') then begin
print,"Reading the specified filter.RES "
readcol,FILTERS_RESPONSES_FILE,FILTERS_WN,format='l',/silent
endif
if FILTER_CONVOL eq 'yes' then begin
dim_w=max(FILTERS_WN)
FILTERS_WAVEL=dblarr(N_FILTERS,dim_w)
FILTERS_RESPONSES=dblarr(N_FILTERS,dim_w)
;FILTERS_LAMBDA_EFF=dblarr(N_FILTERS)
FILTER_WLGHT=dblarr(N_FILTERS) ; Now they are the lambda_eff (here is not an input)
FILTERS_RESPONSES[*,*]=0.
FILTERS_WAVEL[*,*]=0.
; FILTERS_RESPONSES[0,*] ; wavelengths of the first filter
; FILTERS_RESPONSES[1,*] ; response of the first filter
; FILTERS_RESPONSES[2,*] ; wavelengths of the second filter
; FILTERS_RESPONSES[3,*] ; response of the second filter
NFILT=0L
while NFILT lt N_FILTERS do begin
Njumps=0l ; number of jumps + 1 = filter identifier in the filter file
jump=1l   ; lines to jump
while Njumps+1 ne  FILTER_LOG[NFILT] do begin
jump=1+jump+FILTERS_WN[jump-1]
Njumps=Njumps+1
endwhile
;read .res again
readcol,FILTERS_RESPONSES_FILE,FILTERS_WN2,FILT_WAV,RESP,skipline=jump,numline=FILTERS_WN[jump-1],format='l',/silent
; Eliminating duplicated wavelengths (needed for interpolation):
FILT_WAV=FILT_WAV[rem_dup(FILT_WAV)]
RESP=RESP[rem_dup(FILT_WAV)]
FILTERS_WAVEL[NFILT,0:n_elements(FILT_WAV)-1]=FILT_WAV
FILTERS_RESPONSES[NFILT,0:n_elements(FILT_WAV)-1]=RESP
; LAMBDA_EFF FILTERS 
FILTER_WLGHT[NFILT]=int_tabulated(FILT_WAV,RESP*FILT_WAV)/int_tabulated(FILT_WAV,RESP)
print, 'Lambda_eff [A] filter',strcompress(string(NFILT+1)),' : ',strcompress(string(FILTER_WLGHT[NFILT]))
NFILT=NFILT+1
endwhile
;FILTERS CURVES & Lambda_eff  VISUALIZATION
gt0F=where(FILTERS_WAVEL gt 0 and FILTERS_RESPONSES ge 0)
minL=min(FILTERS_WAVEL[gt0F])
maxL=max(FILTERS_WAVEL[gt0F])
minR=min(FILTERS_RESPONSES[gt0F])
maxR=max(FILTERS_RESPONSES[gt0F])
plot, FILTERS_WAVEL[0,*],FILTERS_RESPONSES[0,*],xrange=[minL,maxL],yrange=[minR,maxR],xtitle="wavelength [A]",ytitle="Response",/xlog
NFILT2=0L
while NFILT2 lt N_FILTERS do begin
oplot,FILTERS_WAVEL[NFILT2,*],FILTERS_RESPONSES[NFILT2,*]
oplot, [FILTER_WLGHT[NFILT2],FILTER_WLGHT[NFILT2]],[0,1],linestyle=2,color=1000
NFILT2=NFILT2+1
endwhile

endif ; if convolution with filter response option (FILTER_CONVOL) is YES


if GET_FLUX eq "yes" and  WAV_OR_FILTER eq "filter" then begin
Njumps=0l ; number of jumps + 1 = filter identifier in the filter file
jump=1l   ; lines to jump
while Njumps+1 ne FILTER_GET_FLUX do begin
jump=1+jump+FILTERS_WN[jump-1]
Njumps=Njumps+1
endwhile
;read .res (again)
readcol,FILTERS_RESPONSES_FILE,FILTERS_WN3,FILT_WAV3,RESP3,skipline=jump,numline=FILTERS_WN[jump-1],format='l',/silent
; Eliminating duplicated wavelengths (needed for interpolation):
GET_FLUX_FILT_WAV=FILT_WAV3[rem_dup(FILT_WAV3)]  ; Wavelengths
GET_FLUX_FILT_RESP=RESP3[rem_dup(FILT_WAV3)]     ; Response
;LAMBDA_EFF FILTER (the filter in which we want to get the flux)
GET_FLUX_FILT_L_EFF=int_tabulated(GET_FLUX_FILT_WAV,GET_FLUX_FILT_RESP*GET_FLUX_FILT_WAV)/int_tabulated(GET_FLUX_FILT_WAV,GET_FLUX_FILT_RESP)
print, 'Lambda_eff [A] filter in which getting the flux (FILTER_GET_FLUX) : ',strcompress(string(GET_FLUX_FILT_L_EFF))
;VISUALIZATION
if FILTER_CONVOL eq 'yes' and VISUALIZE_FIT eq 'yes' then begin
   oplot,GET_FLUX_FILT_WAV,GET_FLUX_FILT_RESP,color=60000
   oplot, [GET_FLUX_FILT_L_EFF,GET_FLUX_FILT_L_EFF],[0,1],linestyle=2,color=65000
endif ; just for visualization
if FILTER_CONVOL eq 'no'  and VISUALIZE_FIT eq 'yes' then begin
plot,GET_FLUX_FILT_WAV,GET_FLUX_FILT_RESP,xtitle="wavelength [A]",ytitle="Response",/xlog,yrange=[0,max(GET_FLUX_FILT_RESP)]
oplot, [GET_FLUX_FILT_L_EFF,GET_FLUX_FILT_L_EFF],[0,1],linestyle=2,color=65000





endif
endif ;if GET_FLUX eq "yes" and  WAV_OR_FILTER eq "filter

cont='no'
print, 'Return to continue'
read,cont






; FITS OF THE DATA WITH THE MODELS






;QUANTITIES TO BE MEMORIZED DURING THE COMPARISON (model-source)

; name of the best fit model (contains the information about
; Metallicity, Temperature, Gravity) 
best_fit_model=strarr(n_elements(normalized_cat.field01[0,*]))
; chi^2 of the best solution , second and third
FINAL_CHI2=dblarr(n_elements(normalized_cat.field01[0,*]))
; Number of valid filters that can be used for the fit
N_VALID_FILTERS=lonarr(n_elements(normalized_cat.field01[0,*]))
; Output SED
SED_WAVELENGHTS=dblarr(n_elements(normalized_cat.field01[0,*]))
SED_FLUXES=dblarr(n_elements(normalized_cat.field01[0,*]))
; Flux computed in a filter/wavelength (null if not required as output)
OUT_FLUX_FILTER=dblarr(n_elements(normalized_cat.field01[0,*]))
; Referement MODEL Flux (for normalization). One per source
FLUX_RF_MODEL=dblarr(n_elements(normalized_cat.field01[0,*]))
; Sources extinction 
FINAL_Av=dblarr(n_elements(normalized_cat.field01[0,*]))





;HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE 
; HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE 
;  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE 

;BUILDING HYPERCUBE - Version 1.1


print, "building hypercube"


IF KUR93 eq 'yes' then begin

N_T_models=n_elements(T_selected_idx) ; number of temperature models
N_G_models=n_elements(g_selected_idx) ; number of gravity models
N_M_models=n_elements(METAL_MODEL) ; number of metallicity models
Npoints=0L
;Model reading (to get the maximum number of points in the models)
M_mod=0L
while M_mod lt N_M_models do begin
T_mod=0L
while T_mod lt N_T_models do begin
name1=strmid(METAL_MODEL[M_mod],strlen(METAL_MODEL[M_mod])-4,strlen(METAL_MODEL[M_mod]))
name2=METAL_MODEL[M_mod]+'/'+name1+'_'+T_names[T_mod]+'.fits'
existence_check1=FILE_TEST(name2)
if existence_check1 eq 1 then begin
model=mrdfits(name2,1,HD,/silent)
m_wavelengths=model.(0)
endif
;LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
if n_elements(m_wavelengths) gt Npoints then Npoints=n_elements(m_wavelengths)
;LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
T_mod=T_mod+1
endwhile
M_mod=M_mod+1
endwhile
hypercube=dblarr(N_G_models,N_T_models,N_M_models,Npoints,2)
hypercube[*,*,*,*]=0.

; NOW WE FULFILL THE HYPERCUBE JUST BUILT
M_mod2=0L
while M_mod2 lt N_M_models do begin
T_mod2=0L
while T_mod2 lt N_T_models do begin
G_mod2=0L
while G_mod2 lt N_G_models do begin
;----------------------------------
name1=strmid(METAL_MODEL[M_mod2],strlen(METAL_MODEL[M_mod2])-4,strlen(METAL_MODEL[M_mod2]))
name2=METAL_MODEL[M_mod2]+'/'+name1+'_'+T_names[T_mod2]+'.fits'
existence_check1=FILE_TEST(name2)

if existence_check1 eq 1 then begin
model=mrdfits(name2,1,HD,/silent)
m_tag_nms=TAG_NAMES(model) ; get the tag names in the model
tag_nm_idx=where(m_tag_nms eq g_names[G_mod2]) ;gets the tag name in the structure, corresponding to the first gravity model used.
tag_name=m_tag_nms[tag_nm_idx]
; wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
hypercube[G_mod2,T_mod2,M_mod2,*,0]=model.(0) ; wavelength
hypercube[G_mod2,T_mod2,M_mod2,*,1]=model.(tag_nm_idx) ; model value
; wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
endif
;----------------------------------
G_mod2=G_mod2+1
endwhile
T_mod2=T_mod2+1
endwhile
M_mod2=M_mod2+1
endwhile

ENDIF ; IF KUR93 eq 'yes'



IF KUR93 eq 'no' then begin
Npoints=0.
Nmodels=n_elements(METAL_MODEL)
NMOD=0L
while NMOD lt Nmodels do begin
 model_name=METAL_MODEL[NMOD]
 readcol,model_name,WAV,FLX,format='d,d',/silent
 ;LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 if n_elements(WAV) gt Npoints then Npoints=n_elements(WAV)
 ;LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
NMOD=NMOD+1
endwhile

hypercube=dblarr(Nmodels,Npoints,2)
hypercube[*,*,*]=0.

NMOD2=0L
while NMOD2 lt Nmodels do begin
 model_name=METAL_MODEL[NMOD2]
 readcol,model_name,WAV,FLX,format='d,d',/silent
hypercube[NMOD2,*,0]=WAV
hypercube[NMOD2,*,1]=FLX
NMOD2=NMOD2+1
endwhile


ENDIF ;IF KUR93 eq 'no'



;HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE 
; HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE 
;  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE  HYPERCUBE 



; EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT

; EXTINCTION
Allen_ext=dblarr(22,2)

; Allen (1976) table for K(lambda)/Rv
Allen_ext[0,0]=1000.
Allen_ext[1,0]=1110.
Allen_ext[2,0]=1250.
Allen_ext[3,0]=1430.
Allen_ext[4,0]=1670.
Allen_ext[5,0]=2000.
Allen_ext[6,0]=2220.
Allen_ext[7,0]=2500.
Allen_ext[8,0]=2850.
Allen_ext[9,0]=3330.
Allen_ext[10,0]=3650.
Allen_ext[11,0]=4000.
Allen_ext[12,0]=4400.
Allen_ext[13,0]=5000.
Allen_ext[14,0]=5530.
Allen_ext[15,0]=6700.
Allen_ext[16,0]=9000.
Allen_ext[17,0]=10000.    ; 1 um
Allen_ext[18,0]=20000.
Allen_ext[19,0]=100000.    ; 10 um
Allen_ext[20,0]=1000000.   ; 100 um
Allen_ext[21,0]=10000000.  ; 1000 um

Allen_ext[0,1]=4.2
Allen_ext[1,1]=3.7
Allen_ext[2,1]=3.3
Allen_ext[3,1]=3.0
Allen_ext[4,1]=2.7
Allen_ext[5,1]=2.8
Allen_ext[6,1]=2.9
Allen_ext[7,1]=2.3
Allen_ext[8,1]=1.97
Allen_ext[9,1]=1.69
Allen_ext[10,1]=1.58
Allen_ext[11,1]=1.45
Allen_ext[12,1]=1.32
Allen_ext[13,1]=1.13
Allen_ext[14,1]=1.0
Allen_ext[15,1]=0.74
Allen_ext[16,1]=0.46
Allen_ext[17,1]=0.38
Allen_ext[18,1]=0.11
Allen_ext[19,1]=0.0
Allen_ext[20,1]=0.0
Allen_ext[21,1]=0.0

Emodel_number=long((AV_MAX-AV_MIN)/AV_STEP)

; EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT EXT




if KUR93 eq 'yes' then begin
Tmodel_number=n_elements(T_selected_idx)
Gmodel_number=n_elements(g_selected_idx)
endif

if KUR93 eq 'no' then begin
Tmodel_number=1.
Gmodel_number=1.
endif

;---------------------------------------------------------------------
; TEST CHI2 ON SOURCES - TEST CHI2 ON SOURCES - TEST CHI2 ON SOURCES - 
;---------------------------------------------------------------------
SOURCE=0L
while SOURCE lt n_elements(normalized_cat.field01[0,*]) do begin
if VERBOSE eq 'yes' then print, 'Source number: '+strcompress(SOURCE+1)
; READ the MODELS (TOO MANY TO BE MEMORIZED) AND COMPARE WITH
; THE EFFECTIVE SOURCES, FILTER BY FILTER
MET_MOD=0L
while MET_MOD lt n_elements(METAL_MODEL) do begin
TEM_MOD=0L
while TEM_MOD lt Tmodel_number do begin
GRA_MOD=0L
while GRA_MOD lt Gmodel_number do begin




; RECONSTRUCTION OF THE MODEL NAME
if KUR93 eq 'yes' then begin
Model_effective_name=strmid(METAL_MODEL[MET_MOD],strlen(METAL_MODEL[MET_MOD])-4,strlen(METAL_MODEL[MET_MOD]))
model_name=METAL_MODEL[MET_MOD]+'/'+Model_effective_name+'_'+T_names[TEM_MOD]+'.fits'+'  -   GRAVITY='+G_names[GRA_MOD]
endif
if KUR93 eq 'no' then begin
model_name=METAL_MODEL[MET_MOD]
endif
if VERBOSE eq 'yes' then print, model_name

;ttt='ttt'
;read,ttt

; CHECK MODEL EXISTENCE (changed version 1.1) ***
if KUR93 eq 'no' then are_there_num=where(hypercube[MET_MOD,*,1] gt 0)
if KUR93 eq 'yes' then are_there_num=where(hypercube[GRA_MOD,TEM_MOD,MET_MOD,*,1] gt 0)

if are_there_num[0] eq -1 then check_exist=0
if are_there_num[0] ne -1 then check_exist=1
;**************************************************



if check_exist eq 1 then begin; only if the model is present in the



;____________________________________________________________
; GET MODEL FROM HYPERCUBE

if KUR93 eq 'yes' then begin
if GRA_MOD eq 0 then begin 
MODEL_WAVELENGTH=hypercube[GRA_MOD,TEM_MOD,MET_MOD,*,0] ; NEW V1.1
MODEL_FREQUENCY=3./MODEL_WAVELENGTH ; conversion: [A]-->10^18[Hz]
endif;if GRA_MOD eq 0
ORIGINAL_MODEL_FLUX_no_ext=hypercube[GRA_MOD,TEM_MOD,MET_MOD,*,1] ; NEW V1.1
 endif ;if KUR93 eq 'yes'


if KUR93 eq 'no' then begin
MODEL_WAVELENGTH=hypercube[MET_MOD,*,0]
ORIGINAL_MODEL_FLUX_no_ext=hypercube[MET_MOD,*,1]
endif
;_____________________________________________________________


EXT_MOD=0L
AV_VALUE=AV_MIN
while AV_VALUE le AV_MAX do begin
;while EXT_MOD lt Emodel_number do begin


; XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
; EXTINCTION
if EXTINCT ne 'yes' then begin
ORIGINAL_MODEL_FLUX=ORIGINAL_MODEL_FLUX_no_ext
endif
if EXTINCT eq 'yes' then begin
ALLEN_EXT_INTERP=interpol(Allen_ext[*,1],Allen_ext[*,0],MODEL_WAVELENGTH)
ORIGINAL_MODEL_FLUX=ORIGINAL_MODEL_FLUX_no_ext*10^(-0.4*AV_VALUE*ALLEN_EXT_INTERP)
endif

; XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

; temporary:
;print, SOURCE, MET_MOD,TEM_MOD,GRA_MOD,AV_VALUE


; ARE THERE ELEMENTS IN THE ANALIZED MODEL?
valid_numbers=where(ORIGINAL_MODEL_FLUX ne 0 and ORIGINAL_MODEL_FLUX eq ORIGINAL_MODEL_FLUX)
; IF THERE AREN?T VALID NUMBERS, WE JUMP EVERYTHING, HERE...
if valid_numbers[0] ne -1 then begin


; MODEL PHYSICAL UNITS CONVERSION (IF KURUCZ MODELS)
; ORIGINAL_MODEL_FLUX=model.(tag_name_idx) <-- read above
; IF KURUCZ MODELS, these fluxes are espressed in units of ergs
; cm^{-2} s^{-1} A^{-1}. We want to convert them in units of
; ergs cm^{-2} s^{-1} Hz^{-1}. The conversion is explained in
; Allen's Astrophysical Quantities, pag 183.
; I[ergs cm^-2 s^-1 A^-1]=3.336x(10^-19)xF^2[Hz]xI[ergs cm^-2 s^-1 Hz^-1]
; HERE, the frequencies are expressed in units of 10^18 Hz
; because of memory problems.

if KUR93 eq 'yes' then MODEL_FLUX=(ORIGINAL_MODEL_FLUX/3.336)*10/((MODEL_FREQUENCY)^2); Now the fluxes are in units of [ergs cm^-2 s^-1 Hz^-1] ,in units of 10^18 Hz !!!

if KUR93 eq 'no' then MODEL_FLUX=ORIGINAL_MODEL_FLUX



if FILTER_CONVOL eq 'no' then begin
; MODEL NORMALIZATION

if REFEREMENT_FILTER ne -1 then FLUX_RF_MODEL[SOURCE]=interpol(MODEL_FLUX,MODEL_WAVELENGTH,FILTER_WLGHT[REFEREMENT_FILTER])

 if REFEREMENT_FILTER eq -1 then begin
; model fluxes at the filter wavelengths
FLUX_MODEL_IN_F=interpol(MODEL_FLUX,MODEL_WAVELENGTH,FILTER_WLGHT)
; valid values indexes
ok_numbers=where(normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,SOURCE] gt 0)
; fluxes in the catalog (despite the name is not normalized, in this case)
FLXS_CAT=normalized_cat.field01[1:N_FILTERS,SOURCE]
ERRS_CAT=normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,SOURCE]
ok_numbers1=where(ERRS_CAT gt 0.)
; Referement flux to which normalize: MEAN OF THE FLUXES:
;FLUX_RF_MODEL[SOURCE]=mean(FLUX_MODEL_IN_F[ok_numbers]/FLXS_CAT[ok_numbers])
; Referement flux to which normalize: WEIGTHED MEAN OF THE FLUXES:
FLUX_RF_MODEL[SOURCE]=total(FLUX_MODEL_IN_F[ok_numbers]/FLXS_CAT[ok_numbers])/total(1./FLXS_CAT[ok_numbers])
endif

MODEL_FLUX_NORM=MODEL_FLUX/FLUX_RF_MODEL[SOURCE]
; MODEL NORMALIZED FLUXES: VALUES AT EACH FILTER WAVELENGTH
MODEL_FLUX_WFILTER=interpol(MODEL_FLUX_NORM, MODEL_WAVELENGTH, FILTER_WLGHT)
endif ; no convolution with filters




if FILTER_CONVOL eq 'yes' then begin
MODEL_FLUX_WFILTER=dblarr(N_FILTERS)

if REFEREMENT_FILTER ne -1 then begin
; MODEL NORMALIZATION
; Here we normalize with respect to the CONVOLVED flux in the
; referement filter
; CONVOLUTION MODEL-REF_FILTER: LLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
; interpolation of the spectral model in the referement filter
; wavelengths
idno0_REF=where(FILTERS_WAVEL[REFEREMENT_FILTER,*] gt 0.)
MODEL_INTERP_REF_FIL=interpol(MODEL_FLUX,MODEL_WAVELENGTH,FILTERS_WAVEL[REFEREMENT_FILTER,idno0_REF])
MODEL_IN_REF_FILTER=MODEL_INTERP_REF_FIL*FILTERS_RESPONSES[REFEREMENT_FILTER,idno0_REF]
;---
AREA_MODEL_IN_REF_FILTER=int_tabulated(FILTERS_WAVEL[REFEREMENT_FILTER,idno0_REF],MODEL_IN_REF_FILTER)
;----
AREA_RESP_REF_FILTER=int_tabulated(FILTERS_WAVEL[REFEREMENT_FILTER,idno0_REF],FILTERS_RESPONSES[REFEREMENT_FILTER,idno0_REF])
;-----
;LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
; Referement flux (convolved)
FLUX_RF_MODEL[SOURCE]=AREA_MODEL_IN_REF_FILTER/AREA_RESP_REF_FILTER
; NORMALIZED MODEL FLUX -
MODEL_FLUX_NORM=MODEL_FLUX/FLUX_RF_MODEL[SOURCE]

; MODEL NORMALIZED FLUXES: VALUES in EACH FILTER (CONVOLVED) 
NFILT3=0L
while NFILT3 lt N_FILTERS do begin
no0_idx=where(FILTERS_WAVEL[NFILT3,*] gt 0.)

MODEL_INTERP_FIL=interpol(MODEL_FLUX_NORM,MODEL_WAVELENGTH,FILTERS_WAVEL[NFILT3,no0_idx])
MODEL_IN_FILTER=MODEL_INTERP_FIL*FILTERS_RESPONSES[NFILT3,no0_idx]
AREA_MODEL_IN_FILTER=int_tabulated(FILTERS_WAVEL[NFILT3,no0_idx],MODEL_IN_FILTER)
AREA_RESP_FILTER=int_tabulated(FILTERS_WAVEL[NFILT3,no0_idx],FILTERS_RESPONSES[NFILT3,no0_idx])
; Referement flux (convolved)
MODEL_FLUX_WFILTER[NFILT3]=AREA_MODEL_IN_FILTER/AREA_RESP_FILTER
NFILT3=NFILT3+1
endwhile
endif ; (if REF_FILTER ne -1)


if REFEREMENT_FILTER eq -1 then begin
NFILT4=0L
while NFILT4 lt N_FILTERS do begin
no0_idx=where(FILTERS_WAVEL[NFILT4,*] gt 0.)
MODEL_INTERP_FIL=interpol(MODEL_FLUX,MODEL_WAVELENGTH,FILTERS_WAVEL[NFILT4,no0_idx])
MODEL_IN_FILTER=MODEL_INTERP_FIL*FILTERS_RESPONSES[NFILT4,no0_idx]
AREA_MODEL_IN_FILTER=int_tabulated(FILTERS_WAVEL[NFILT4,no0_idx],MODEL_IN_FILTER)
AREA_RESP_FILTER=int_tabulated(FILTERS_WAVEL[NFILT4,no0_idx],FILTERS_RESPONSES[NFILT4,no0_idx])
; Referement flux (convolved)
MODEL_FLUX_WFILTER[NFILT4]=AREA_MODEL_IN_FILTER/AREA_RESP_FILTER
NFILT4=NFILT4+1
endwhile
;LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
; Referement flux (convolved) simple mean
;;;;;;;;;;;
;FLUX_RF_MODEL[SOURCE]=mean(MODEL_FLUX_WFILTER)
;;;;;;;;;;;
; Referement flux (convolved) weigthed  mean
ERRS_CAT2=normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,SOURCE]
ok_numbers2=where(ERRS_CAT2 gt 0.)
FLUX_RF_MODEL[SOURCE]=total(MODEL_FLUX_WFILTER[ok_numbers2]/ERRS_CAT2[ok_numbers2])/total(1./ERRS_CAT2[ok_numbers2])
;;;;;;;;;;;

; NORMALIZED MODEL FLUX -
MODEL_FLUX_NORM=MODEL_FLUX/FLUX_RF_MODEL[SOURCE]
MODEL_FLUX_WFILTER=MODEL_FLUX_WFILTER/FLUX_RF_MODEL[SOURCE]

endif  ; (if REF_FILTER eq -1)

endif ; (if FILTER_CONVOL eq 'yes' then begin)






; TO CONVOLVE FILTERS WITH FLUXES: EXAMPLE. THIS PART OF CODE IS NOT
; USED IN THIS PROGRAM: it is just an example
example='no'
if example eq  'do_it' then begin
;read filter transmission
readcol,'~/filters/24um.txt',x15,y15
;reade spectrum
readcol,'/templates/Sey2_template_norm.sed',x,y
res=interpol(y,x,x15)
f_conv=res*y15
f15l=int_tabulated(x15,f_conv)
flw3=int_tabulated(x15,y15)
y15em=f15l/flw3  ;convolved flux
endif


;CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE
;CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE
; VARIANCE SOURCE-MODEL
  DIFF=fltarr(N_FILTERS)
  SIGMA=fltarr(N_FILTERS)
  DIFF=MODEL_FLUX_WFILTER-normalized_cat.field01[1:N_FILTERS,SOURCE]
  ; Associated uncertainties (from the catalog):
  SIGMA=normalized_cat.field01[N_FILTERS+1:2*N_FILTERS,SOURCE]
  no99_idx=where(normalized_cat.field01[1:N_FILTERS,SOURCE] gt 0)
  CHI2=total((DIFF[no99_idx]/SIGMA[no99_idx])^2)
;CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE
;CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE-CORE
  N_VALID_FILTERS[SOURCE]=n_elements(no99_idx)


; COMPUTING CHI2 (AND MEMORIZING THE BEST ONE)
; first model is considered the best (at least before reading the second)
if GRA_MOD eq 0 and TEM_MOD eq 0 and MET_MOD eq 0 then FINAL_CHI2[SOURCE]=CHI2

if CHI2 eq CHI2 then FINAL_CHI2[SOURCE]=min([FINAL_CHI2[SOURCE],CHI2])
if CHI2 eq CHI2 and FINAL_CHI2[SOURCE] ne FINAL_CHI2[SOURCE] then FINAL_CHI2[SOURCE]=CHI2

if VERBOSE eq 'yes' then begin
print,'Av value ->',strcompress(AV_VALUE)
print,'Chi^2 This model->',strcompress(CHI2)
print,'Chi^2 best model->', strcompress(FINAL_CHI2[SOURCE])
endif
if FINAL_CHI2[SOURCE] eq CHI2 then begin
;EXTINCTION
FINAL_Av[SOURCE]=AV_VALUE
; MEMORIZE MODEL NAME AND CHARACTERISTICS
if KUR93 eq 'yes' then best_fit_model[SOURCE]=Model_effective_name+'_'+T_names[TEM_MOD]+'_'+G_names[GRA_MOD]
if KUR93 eq 'no' then best_fit_model[SOURCE]=MET_MOD
SED_WAVELENGHTS=MODEL_WAVELENGTH
SED_FLUXES=MODEL_FLUX_NORM
; MODEL GRAPHICAL REPRESENTATION
titleID='Source '+strcompress(string(normalized_cat.field01[0,SOURCE]))
vect_provv=normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,SOURCE]
idno=where(vect_provv gt 0.)    ;VVV
min_val=min(vect_provv[idno])/10. ;VVV
max_val=max(vect_provv[idno])*10. ;VVV

if VISUALIZE_FIT eq 'yes' then begin
plot,MODEL_WAVELENGTH,MODEL_FLUX_NORM,/xlog,/ylog,title=titleID,ytitle='NORMALIZED FLUX [ergs cm-2 s-1 Hz-1]',xtitle='LAMBDA[A]',charsize=1.5,xr=[minL/LEFT_EXPANSION,maxL/RIGHT_EXPANSION],yr=[min_val/LOW_EXPANSION,max_val*UP_EXPANSION],/xst,/yst


oplot,FILTER_WLGHT, normalized_cat.field01[1:N_FILTERS,SOURCE],psym=6,color=1000,thick=2.5
FN=0l
while FN lt n_elements(FILTER_WLGHT) do begin
if normalized_cat.field01[n_elements(FILTER_WLGHT)+1+FN,SOURCE] gt 0 then begin
oploterror,[FILTER_WLGHT[FN],FILTER_WLGHT[FN]],[normalized_cat.field01[FN+1,SOURCE],normalized_cat.field01[FN+1,SOURCE]],[0.,0.],[normalized_cat.field01[n_elements(FILTER_WLGHT)+1+FN,SOURCE],normalized_cat.field01[n_elements(FILTER_WLGHT)+1+FN,SOURCE]],errcolor=1000,thick=2.5
endif
FN=FN+1
endwhile
endif ;if VISUALIZE_FIT eq 'yes' 

; GET FLUX IN SELECTED FILTER (OUT_FLUX_FILTER). SAME units as output sed.
if GET_FLUX eq "yes" then begin
 if SED_TYPE eq 'JY' or SED_TYPE eq 'Jy' or SED_TYPE eq 'AB_MAG' then cost1=FLUX_RF_SOURCE[SOURCE]
 if SED_TYPE eq 'NORM' then cost1=1.

; FROM WAVELENGTH
if WAV_OR_FILTER eq "wavelength" then begin
OUT_FLUX_FILTER[SOURCE]=interpol(MODEL_FLUX_NORM*cost1, MODEL_WAVELENGTH, FILTER_GET_FLUX)
visual_extr_flux=OUT_FLUX_FILTER[SOURCE]/cost1
WAV_GET_FILTER=double(FILTER_GET_FLUX)
endif ; if WAV_OR_FILTER eq "wavelength"

;FROM FILTER CONVOLUTION
if WAV_OR_FILTER eq "filter" then begin
; CONVOLUTION MODEL-GET_FLUX_FILTER: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; interpolation of the spectral model in the referement filter
; wavelengths
idno0_GETF=where(GET_FLUX_FILT_WAV gt 0.)
MODEL_INTERP_GET_FIL=interpol(MODEL_FLUX_NORM,MODEL_WAVELENGTH,GET_FLUX_FILT_WAV[idno0_GETF])
MODEL_IN_GET_FILTER=MODEL_INTERP_GET_FIL*GET_FLUX_FILT_RESP[idno0_GETF]
;---
AREA_MODEL_IN_GET_FILTER=int_tabulated(GET_FLUX_FILT_WAV[idno0_GETF],MODEL_IN_GET_FILTER)
;----
AREA_RESP_GET_FILTER=int_tabulated(GET_FLUX_FILT_WAV[idno0_GETF],GET_FLUX_FILT_RESP[idno0_GETF])
;-----
; Flux (convolved) in the filter
FLUX_GETF_MODEL=AREA_MODEL_IN_GET_FILTER/AREA_RESP_GET_FILTER
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
WAV_GET_FILTER=GET_FLUX_FILT_L_EFF ; effective wavelength - for visualization
visual_extr_flux=FLUX_GETF_MODEL   ; convolved flux - for visualization
OUT_FLUX_FILTER[SOURCE]=FLUX_GETF_MODEL*cost1; convolved flux - for output
endif ;if WAV_OR_FILTER eq "filter"

; VISUALIZE THE FLUX IN THE SELECTED FILTER/WAVELENGTH (normalized)
if VISUALIZE_FIT eq 'yes' then oplot,[WAV_GET_FILTER,WAV_GET_FILTER], [visual_extr_flux,visual_extr_flux],psym=6,color=60000,thick=3.,symsize=1.5
endif ;if GET_FLUX eq 'yes' 


; WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
; WRITE OUTPUT SEDs TXT
; WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
if OUT_BEST_FIT_SED_TXT eq 'yes' then begin
SED_NAME=OUTPUT_PATH+strcompress(string(long(normalized_cat.field01[0,SOURCE])),/remove_all)+'.sed'
openw, fo1,SED_NAME, /get_lun
m=0l
while m lt n_elements(MODEL_WAVELENGTH) do begin
if SED_TYPE eq 'Jy' or SED_TYPE eq 'JY' then printf,fo1,MODEL_WAVELENGTH[m],MODEL_FLUX_NORM[m]*FLUX_RF_SOURCE[SOURCE],format='(i10,1x,f16.10)' 
if SED_TYPE eq 'NORM' then printf,fo1,MODEL_WAVELENGTH[m],MODEL_FLUX_NORM[m],format='(i10,1x,f16.7)'
if SED_TYPE eq 'AB_MAG' then begin
out_sed_mag=8.9-2.5*alog10(MODEL_FLUX_NORM[m]*FLUX_RF_SOURCE[SOURCE])
if MODEL_FLUX_NORM[m] eq 0 then out_sed_mag=99.
printf,fo1,MODEL_WAVELENGTH[m],out_sed_mag,format='(i10,1x,f16.8)'
endif
m=m+1
endwhile
free_lun, fo1
close, fo1
endif ; if OUT_BEST_FIT_SED_TXT eq 'yes'
; WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
; WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

; SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
; WRITE OUTPUT SEDs PS
; SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
if OUT_BEST_FIT_SED_PS eq 'yes' then begin

SED_NAME_PS=OUTPUT_PATH+strcompress(string(long(normalized_cat.field01[0,SOURCE])),/remove_all)+'.ps'
rapp=1.
p_x_f=10.*PS_XSIZE
p_y_f=rapp*10.*PS_YSIZE
;loadct,3
set_plot,'ps',/copy
device,/ color
TVLCT,[0,255,0,0],[0,0,255,0],[0,0,0,255]
;nero=0
;red=1
;green=2
;blue=3
device, filename=SED_NAME_PS,/color, Xsize=p_x_f, Ysize=p_y_f, /cm

title1='Source '+strcompress(string(long(normalized_cat.field01[0,SOURCE])),/remove_all)
XX_model=MODEL_WAVELENGTH
XX_TITLE='Wavelength [A]'
vect_provv=normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,SOURCE] ;errors
idno=where(vect_provv gt 0.)    ;VVV

if SED_TYPE eq 'Jy' or SED_TYPE eq 'JY' then begin
YY_model=MODEL_FLUX_NORM*FLUX_RF_SOURCE[SOURCE]
YY_FILTERS=normalized_cat.field01[1:N_FILTERS,SOURCE]*FLUX_RF_SOURCE[SOURCE]
YY_ERRS=original_Jy_cat.field01[N_FILTERS+1:N_FILTERS*2,SOURCE]
YY_title='Flux [Jy]'
min_valps=min(YY_FILTERS[idno])/5 ;VVV
max_valps=max(YY_FILTERS[idno])*5 ;VVV
plot,XX_model,YY_model,title=title1,ytitle=YY_title,xtitle=XX_TITLE,xrange=[minL/LEFT_EXPANSION,maxL*RIGHT_EXPANSION],yrange=[min_valps/LOW_EXPANSION,max_valps*UP_EXPANSION],/xlog,/ylog,charsize=PS_CHSIZE,charthick=PS_THICK,xthick=PS_THICK,ythick=PS_THICK,thick=PS_THICK,/xst,/yst
oplot, FILTER_WLGHT,YY_FILTERS,color=1,psym=6
FN2=0l
while FN2 lt n_elements(FILTER_WLGHT) do begin
if normalized_cat.field01[n_elements(FILTER_WLGHT)+1+FN2,SOURCE] gt 0 then begin
oploterror,[FILTER_WLGHT[FN2],FILTER_WLGHT[FN2]],[YY_FILTERS[FN2],YY_FILTERS[FN2]],[0.,0.],[YY_ERRS[FN2],YY_ERRS[FN2]],errcolor=1,thick=2.5
endif
FN2=FN2+1
endwhile
endif

if SED_TYPE eq 'NORM' then begin
YY_model=MODEL_FLUX_NORM
YY_FILTERS=normalized_cat.field01[1:N_FILTERS,SOURCE]
YY_ERRS=normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,SOURCE]
YY_title="Normalized values"
min_valps=min(YY_FILTERS[idno])/5 ;VVV
max_valps=max(YY_FILTERS[idno])*5 ;VVV
plot,XX_model,YY_model,title=title1,ytitle=YY_title,xtitle=XX_TITLE,xrange=[minL/LEFT_EXPANSION,maxL*RIGHT_EXPANSION],yrange=[min_valps/LOW_EXPANSION,max_valps*UP_EXPANSION],/ylog,/xlog,charsize=PS_CHSIZE,charthick=PS_THICK,xthick=PS_THICK,ythick=PS_THICK,thick=PS_THICK,/xst,/yst
oplot, FILTER_WLGHT,YY_FILTERS,color=1,psym=6
FN2=0l
while FN2 lt n_elements(FILTER_WLGHT) do begin
if normalized_cat.field01[n_elements(FILTER_WLGHT)+1+FN2,SOURCE] gt 0 then begin
oploterror,[FILTER_WLGHT[FN2],FILTER_WLGHT[FN2]],[YY_FILTERS[FN2],YY_FILTERS[FN2]],[0.,0.],[YY_ERRS[FN2],YY_ERRS[FN2]],errcolor=1,thick=2.5
endif
FN2=FN2+1
endwhile
endif

if SED_TYPE eq 'AB_MAG' then begin
YY_model=8.9-2.5*alog10(MODEL_FLUX_NORM*FLUX_RF_SOURCE[SOURCE])
YY_FILTERS=8.9-2.5*alog10(normalized_cat.field01[1:N_FILTERS,SOURCE]*FLUX_RF_SOURCE[SOURCE])
YMAG_SUP=8.9-2.5*alog10((normalized_cat.field01[1:N_FILTERS,SOURCE]-normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,SOURCE])*FLUX_RF_SOURCE[SOURCE])
YMAG_INF=8.9-2.5*alog10((normalized_cat.field01[1:N_FILTERS,SOURCE]+normalized_cat.field01[N_FILTERS+1:N_FILTERS*2,SOURCE])*FLUX_RF_SOURCE[SOURCE])
YY_ERRS=(YMAG_SUP-YMAG_INF)/2
YY_title="mag AB"
min_valps=max(YY_FILTERS[idno])+2 ;VVV
max_valps=min(YY_FILTERS[idno])-2 ;VVV
plot,XX_model,YY_model,title=title1,ytitle=YY_title,xtitle=XX_TITLE,xrange=[minL/LEFT_EXPANSION,maxL*RIGHT_EXPANSION],yrange=[min_valps+(alog10(LOW_EXPANSION)*2.5),max_valps-(alog10(UP_EXPANSION)*2.5)],/xlog,charsize=PS_CHSIZE,charthick=PS_THICK,xthick=PS_THICK,ythick=PS_THICK,thick=PS_THICK,/xst,/yst
oplot, FILTER_WLGHT,YY_FILTERS,color=1,psym=6,thick=2.5
FN2=0l
while FN2 lt n_elements(FILTER_WLGHT) do begin
if normalized_cat.field01[n_elements(FILTER_WLGHT)+1+FN2,SOURCE] gt 0 then begin
oploterror,[FILTER_WLGHT[FN2],FILTER_WLGHT[FN2]],[YY_FILTERS[FN2],YY_FILTERS[FN2]],[0.,0.],[YY_ERRS[FN2],YY_ERRS[FN2]],errcolor=1,thick=2.5
endif
FN2=FN2+1
endwhile
endif

if GET_FLUX eq 'yes' then begin
if SED_TYPE eq 'AB_MAG' then out_filter_get_flux_vis=8.9-2.5*alog10(OUT_FLUX_FILTER[SOURCE])
if SED_TYPE eq 'JY' or SED_TYPE eq 'Jy' then out_filter_get_flux_vis=OUT_FLUX_FILTER[SOURCE]
if SED_TYPE eq 'NORM' then out_filter_get_flux_vis=OUT_FLUX_FILTER[SOURCE];/FLUX_RF_SOURCE[SOURCE]
oplot,[WAV_GET_FILTER,WAV_GET_FILTER], [out_filter_get_flux_vis,out_filter_get_flux_vis],psym=6,color=2,thick=4.,symsize=1.5
endif

device,/close
set_plot,'x',/copy
endif ; if OUT_BEST_FIT_SED_PS eq 'yes'
; SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
; SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

endif ; if FINAL_CHI2[SOURCE] eq CHI2 

endif ; if valid_numbers[0] ne -1 --> the gravity model is not present.



AV_VALUE=AV_VALUE+AV_STEP
EXT_MOD=EXT_MOD+1 ; extinction (from internal table)
endwhile

endif ;if the metalicity-temperature model is not present, jumps here

GRA_MOD=GRA_MOD+1 ; gravity (Kurucz '93 models)
endwhile
TEM_MOD=TEM_MOD+1 ; temperature (Kurucz '93 models)
endwhile
MET_MOD=MET_MOD+1 ; metallicity (Kurucz '93 models)
endwhile
SOURCE=SOURCE+1   ; source number from catalog
endwhile
;---------------------------------------------------------------------
; TEST CHI2 ON SOURCES - TEST CHI2 ON SOURCES - TEST CHI2 ON SOURCES - 
;------------------------------- END ---------------------------------




; TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
; WRITE OUTPUT FILE1 (info on single sources)
; TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

openw, fo1, OUTPUT_FILE1, /get_lun
nn=0l
while nn lt n_elements(normalized_cat.field01[0,*]) do begin

if EXTINCT ne 'yes' then extnct=0.
if EXTINCT eq 'yes' then extnct=FINAL_Av[nn]

if GET_FLUX eq 'no' then begin
out_filter_get_flux=normalized_cat.field01[0,*]
out_filter_get_flux[*]=-99.
endif
if SED_TYPE eq 'AB_MAG' and GET_FLUX eq 'yes' then out_filter_get_flux=8.9-2.5*alog10(OUT_FLUX_FILTER[nn])
if (SED_TYPE eq 'JY' or SED_TYPE eq 'Jy') and GET_FLUX eq 'yes' then out_filter_get_flux=OUT_FLUX_FILTER[nn]
if SED_TYPE eq 'NORM' and GET_FLUX eq 'yes' then out_filter_get_flux=OUT_FLUX_FILTER[nn]/FLUX_RF_SOURCE[nn]

IF KUR93 eq 'yes' then begin
    if strmid(strcompress(best_fit_model[nn],/remove_all),1,1) eq 'm' then sign_mtl=-1.
    if strmid(strcompress(best_fit_model[nn],/remove_all),1,1) eq 'p' then sign_mtl=+1.
    METALLICITY=sign_mtl*0.1*float(strmid(strcompress(best_fit_model[nn],/remove_all),2,2))
    if strmid(strcompress(best_fit_model[nn],/remove_all),9,1) eq '_' then hh=4
    if strmid(strcompress(best_fit_model[nn],/remove_all),10,1) eq '_' then hh=5
    TEMPERATURE=long(strmid(strcompress(best_fit_model[nn],/remove_all),5,hh))
    GRAVITY=0.1*float(strmid(strcompress(best_fit_model[nn],/remove_all),strlen(best_fit_model[nn])-2,2))
    printf,fo1,normalized_cat.field01[0,nn],FINAL_CHI2[nn],best_fit_model[nn],N_VALID_FILTERS[nn],TEMPERATURE,METALLICITY,GRAVITY,out_filter_get_flux,extnct,format='(i8,1x,f12.6,1x,a15,1x,i5,1x,i5,1x,f5.2,1x,f5.2,1x,f16.8,1x,f12.6,1x)'
 ENDIF

IF KUR93 eq 'no' then begin
printf,fo1,normalized_cat.field01[0,nn],FINAL_CHI2[nn],best_fit_model[nn]+1,N_VALID_FILTERS[nn],out_filter_get_flux,extnct,format='(i8,1x,f12.6,1x,i8,1x,i5,1x,f16.8,1x,f12.6,1x)'
ENDIF

nn=nn+1
endwhile
free_lun, fo1
close, fo1
; TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
; TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


; OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
; WRITE OUTPUT FILE2 (general info on the fits performed)
; OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
openw, fo2, OUTPUT_FILE2, /get_lun
printf,fo2,"##########################################################"
printf,fo2, "# This is the output file of BESTFIT - version ",version
printf,fo2,"##########################################################"
printf,fo2,"DATE: ", SYSTIME()
printf,fo2,"Photometric catalogue :"
printf,fo2,CATALOG_NAME
printf,fo2,"# Templates SEDs: "

hh1=0l
while hh1 lt n_elements(METAL_MODEL) do begin
printf,fo2,strcompress(string(hh1+1))+' '+METAL_MODEL[hh1]
hh1=hh1+1
endwhile

printf,fo2,"# Characteristics of filters: "
if FILTER_CONVOL eq 'yes' then printf,fo2,"# n   wl_eff"
if FILTER_CONVOL eq 'no' then printf,fo2,"# Characteristic wavelength"
hh2=0l
while hh2 lt N_FILTERS do begin
 if FILTER_CONVOL eq 'yes' then begin
  if hh2 eq REFEREMENT_FILTER then string0='  <<<'
  if hh2 ne REFEREMENT_FILTER then string0=' '
 string=strcompress(string(FILTER_LOG[hh2]))+' '+strcompress(string(FILTER_WLGHT[hh2]))+string0

 printf,fo2,string
 endif
 if FILTER_CONVOL eq 'no' then begin
  if hh2 eq REFEREMENT_FILTER then string0='  <<<'
  if hh2 ne REFEREMENT_FILTER then string0=' '
  string=strcompress(string(FILTER_WLGHT[hh2]))+string0
  printf,fo2,string
 endif
hh2=hh2+1
endwhile

if REFEREMENT_FILTER ne -1 then printf,fo2,"# * The reference filter highlighted with '<<<' "
if REFEREMENT_FILTER eq -1 then  printf,fo2,"# No reference filter present"

if GET_FLUX eq 'yes' then begin
printf,fo2,"# Filter in which the flux is computed"
 if WAV_OR_FILTER eq 'filter' then begin
  printf,fo2,"# n   wl_eff"
  printf,fo2,strcompress(string(FILTER_GET_FLUX))+' '+strcompress(string(GET_FLUX_FILT_L_EFF))
 endif
 if WAV_OR_FILTER eq 'wavelength' then begin
  printf,fo2,"# Characteristic wavelength"
  printf,fo2,strcompress(string(FILTER_GET_FLUX))
 endif
endif

if KUR93 eq 'yes' then begin
 printf,fo2, '# Gravity (log_g relative to solar) min: '+strcompress(GRAVITY_MIN)
 printf,fo2, '# Gravity (log_g relative to solar) max: '+ strcompress(GRAVITY_MAX)
 printf,fo2, '# Number of gravity models (available, used): '+strcompress(n_elements(g_available))+' '+strcompress(n_elements(g_selected_idx))
 printf,fo2, '# T min: '+string(T_MIN)+" K"
 printf,fo2, '# T max: '+string(T_MAX)+" K"
 printf,fo2, '# Number of temperature models (available, used): '+strcompress(n_elements(T_available))+' '+strcompress(n_elements(T_selected_idx))
endif

if EXTINCT eq 'yes' then begin
printf,fo2,'# Minimum Av: ',strcompress(string(AV_MIN))
printf,fo2,'# Maximum Av: ',strcompress(string(AV_MAX))
printf,fo2,'# Av steps  : ',strcompress(string(AV_STEP))
endif



free_lun, fo2
close, fo2
; OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOMODEL_FREQUENCY=3./MODEL_WAVELENGTH ; conversion: [A]-->10^18[Hz]OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
; OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

stop
end
