module tm_module

use fml_lib
use netcdf

implicit none
contains

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!
subroutine setup_model()

! local variables
integer::n 
character(8)::date
character(10)::time
character(5)::zone
integer,dimension(8)::value

call date_and_time(date,time,zone,value)
! print header text to screen
print*,'**************************'
print*,'*** Fortran Matrix Lab ***'
print*,'**************************'
print*,''
print '(1x,I4,A1,I2,A1,I2)',value(1),'/',value(2),'/',value(3)
print '(1x,I2,A1,I2,A1,I2)',value(5),':',value(6),':',value(7)
print*,''
print*,'*************************'
print*,'Input Files:'
print*,'../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aexp_filename)
print*,'../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aimp_filename)
print*,trim(tm_Aremin_filename)
print*,'../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_seaice_filename)
print*,'../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_vol_filename)
if(bg_PO4restore_select)then
print*,'../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_PO4restore_filename)
else
print*,'../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_PO4uptake_filename)
end if
PRINT*,
if(gen_restart_select)then
print*,'Restart File:'
print*,'../output/'//trim(gen_restart_filename)//'.nc'
PRINT*,
end if
print*,'Output File:'
print*,'../output/'//trim(gen_config_filename)//'.nc'
Print*,
print*,'Parameters:'
print*,'Number of TMM timesteps per year:',tm_n_dt
print*,'Biogeochemistry timestep ratio:',bg_dt_ratio
print*,'Seasonal transport:',tm_seasonal
print*,'Number of euphotic zone depth layers:',bg_n_euphotic_lyrs
PRINT*,'PO4 nutrient restoring scheme:',bg_PO4restore_select
print*,'PO4 uptake timescale:',bg_uptake_tau
PRINT*,'DOC export fraction:',bg_DOC_frac
PRINT*,'DOC remin rate:',bg_DOC_k
print*,'*************************'
print*,

! create indices for calculating seasonal TMs
call calc_seasonal_scaling()
call load_data_saving()

! convert parameters to correct units
tm_dt=1.0/tm_n_dt ! TMM timestep
bg_dt=tm_dt*bg_dt_ratio ! BGC timestep

bg_DOC_rfrac=1.0-bg_DOC_frac ! reciprical of DOC fraction
!bg_DOC_k=1.0/bg_DOC_k ! year-1
bg_uptake_tau=(1.0/bg_uptake_tau)*gen_conv_d_yr ! days to years-1

! for defining iSur
if(bg_n_euphotic_lyrs.eq.1)then
	n_euphotic_boxes=4448
elseif(bg_n_euphotic_lyrs.eq.2)then
	n_euphotic_boxes=8840
end if

n_surface_boxes=4448

! tracer indices
ioPO4=1
ioDOP=2
ioDIC=3
ioALK=4

iK1=1
iK2=2
iK0=3
iKB=4
iKw=5
iKSi=6
iKp1=7
iKp2=8
iKp3=9

ioCO2=1
ioCO3=2
ioH=3
iopCO2=4

iaO2=1
iaCO2=2

n_ATM_tracers=2 ! n.b. 1 does not allow array operations
!if(bg_O_select) n_ATM_tracers=n_ATM_tracers+1
!IF(bg_C_select) n_ATM_tracers=n_ATM_tracers+1

! set-up sparse matrices
Aexp%nnz=tm_Aexp_nnz
Aimp%nnz=tm_Aimp_nnz
Aremin%nnz=tm_Aremin_nnz

allocate(Aexp%val(Aexp%nnz,n_seasonal))
allocate(Aexp%row(tm_nbox+1)) 
allocate(Aexp%col(Aexp%nnz))

allocate(Aimp%val(Aimp%nnz,n_seasonal))
allocate(Aimp%row(tm_nbox+1))
allocate(Aimp%col(Aimp%nnz))

allocate(Aremin%val(Aremin%nnz,1)) 
allocate(Aremin%row(tm_nbox+1))
allocate(Aremin%col(Aremin%nnz))

allocate(tracers(tm_nbox,gen_n_tracers))
allocate(tracers_1(tm_nbox,gen_n_tracers))
allocate(C(tm_nbox,3))
allocate(C_consts(tm_nbox,9))
allocate(ATM(n_ATM_tracers))
allocate(J(tm_nbox,gen_n_tracers))
allocate(particles(tm_nbox,gen_n_tracers))
allocate(export(tm_nbox))
allocate(tracers_int(tm_nbox,gen_n_tracers))
allocate(EXPORT_int(tm_nbox,n_seasonal))
allocate(ATM_int(n_ATM_tracers))

allocate(iSur(n_euphotic_boxes))
allocate(tm_seaice_frac(n_euphotic_boxes,n_seasonal))
allocate(tm_windspeed(n_euphotic_boxes,n_seasonal))
allocate(tm_area(tm_nbox))
allocate(tm_T(n_euphotic_boxes,n_seasonal))
allocate(tm_S(n_euphotic_boxes,n_seasonal))
allocate(tm_silica(n_euphotic_boxes,n_seasonal))
allocate(bg_PO4_obs(n_euphotic_boxes,n_seasonal))
allocate(bg_PO4_uptake(n_euphotic_boxes,n_seasonal))

allocate(seaice_dt(n_euphotic_boxes))
allocate(wind_dt(n_euphotic_boxes))
allocate(T_dt(n_euphotic_boxes))
allocate(S_dt(n_euphotic_boxes))
allocate(silica_dt(n_euphotic_boxes))


allocate(tm_vol(tm_nbox))

! surface indices (are the first 4448 entries to the vector)
do n=1,n_euphotic_boxes
	iSur(n)=n
end do


! initialise tracer array
if(gen_restart_select)then
else
	tracers_1(:,ioPO4)=bg_PO4_init
	tracers_1(:,ioDOP)=bg_DOC_init
	tracers_1(:,ioDIC)=bg_DIC_init
	tracers_1(:,ioALK)=bg_ALK_init
	end if
tracers(:,:)=0.0
J(:,:)=0.0
tracers_int(:,:)=0.0
ATM_int(:)=0.0
EXPORT_int(:,:)=0.0
dt_count=1 ! keep track of how many timesteps have passed in one year

!Sc_coeffs(:,iaO2)=(/1953.4 , 128.00 , 3.9918 , 0.050091/)   ! O2
!SC_coeffs(:,iaCO2)=(/2073.1 , 125.62 , 3.6276 , 0.043219/)   ! CO2

Sc_coeffs(:,iaCO2)=(/2116.8 , -136.25 , 4.7353 , -0.092307 , 0.0007555/)   ! CO2 Wanninkhof (2014), Orr et al., 2017, Table 1
SC_coeffs(:,iaO2)=(/1920.4 , -135.6 , 5.2122 , -0.10939 , 0.00093777/)   ! O2 Wanninkhof (2014), Orr et al., 2017, Table 1
      
Bunsen_coeffs(:,iaO2)=(/-58.3877 , 85.8079 , 23.8439 , -0.034892 , 0.015568 , -0.0019387/) ! O2
Bunsen_coeffs(:,iaCO2)=(/ -60.2409 , 93.4517 , 23.3585 , 0.0023517 , -0.023656 , 0.0047036 /) ! CO2

Sol_Orr(:,iaCO2)=(/-160.7333 , 215.4152 , 89.8920 , -1.47759 , 0.029941 , -0.027455 , 0.0053407 /) ! CO2, Orr et al., 2017, Table 2
!Sol_Orr(:,iaCO2)=Sol_Orr(:,iaCO2)*1000.0 ! mol L-1 atm-1 -> mol m-3 atm-1 (Orr et al., 2017, Table 2)

!Sol_Orr(:,iaCO2)=(/-162.8301 , 218.2968 , 90.9241 , -1.47696 , 0.025695 , -0.025225 , 0.0049867/) ! CO2, Orr et al., 2017, Table 2 as mol kg-1 for checking

ATM_vol=7777.0 * sum(tm_area(1:n_surface_boxes)) ! height (m) * total area (m2)
ATM_mol=1.77e20 ! 

! temporary code
C(:,:)=-1.0
C_consts(:,:)=0.0
!tm_T(:,:)=20.0
!tm_S(:,:)=30.0
!tm_windspeed(:,:)=0.0
tm_area=1.0
ATM(iaCO2)=278.0*1.0e-6

call initialise_output()


end subroutine setup_model

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine calc_seasonal_scaling

! local variables
integer::n,count,nn

if(tm_seasonal)then
	n_seasonal=12
else
	n_seasonal=1
end if

allocate(tm_seasonal_scale(tm_n_dt))
allocate(tm_seasonal_rscale(tm_n_dt))
allocate(tm_seasonal_n1(tm_n_dt))
allocate(tm_seasonal_n2(tm_n_dt))

tm_seasonal_n1=0.0
tm_seasonal_n2=0.0

if(tm_seasonal)then

count=1
do n=1,tm_n_dt,tm_n_dt/12
	do nn=0,(tm_n_dt/12)-1
		tm_seasonal_n1(n+nn)=count
		tm_seasonal_n2(n+nn)=count+1
	end do
	count=count+1
end do
where(tm_seasonal_n2==13) tm_seasonal_n2=1

do n=1,tm_n_dt,tm_n_dt/12
	DO nn=0,(tm_n_dt/12)-1
		tm_seasonal_rscale(nn+n)=real(nn)/(real(tm_n_dt/12)-1)
		tm_seasonal_scale(nn+n)=1.0-tm_seasonal_rscale(nn+n)
	end do
end do
		
else
	count=1
	do n=1,tm_n_dt

		tm_seasonal_n1(n)=count
		tm_seasonal_n2(n)=count
		
		tm_seasonal_scale(n)=1.0
		tm_seasonal_rscale(n)=1.0-tm_seasonal_scale(n)
	
	end do
	
end if

!print*,'n1',tm_seasonal_n1
!print*,'n2',tm_seasonal_n2
!print*,'scale',tm_seasonal_scale
!print*,'rscale',tm_seasonal_rscale


!~ ! local variables
!~ integer::n,count,nn
!~ integer,dimension(13)::jul_cumdays
!~ real,dimension(12)::jul_days
!~ real::n_days,sum_val

!~ n_days=1.0/((1.0/tm_n_dt)*365) ! number of days per timestep
!~ jul_days=(/31,28,31,30,31,30,31,31,30,31,30,31/)*int(n_days) ! julian days

!~ ! create indices
!~ jul_cumdays(1)=1
!~ sum_val=1
!~ do n=2,13
	!~ sum_val=sum_val+jul_days(n-1)
	!~ jul_cumdays(n)=sum_val
!~ end do
	
	
!~ if(tm_seasonal)then
	!~ n_seasonal=12
!~ else
	!~ n_seasonal=1
!~ end if

!~ allocate(tm_seasonal_scale(tm_n_dt))
!~ allocate(tm_seasonal_rscale(tm_n_dt))
!~ allocate(tm_seasonal_n1(tm_n_dt))
!~ allocate(tm_seasonal_n2(tm_n_dt))

!~ tm_seasonal_n1=0.0
!~ tm_seasonal_n2=0.0
!~ if(tm_seasonal)then

	!~ count=1
	!~ do n=1,12
		!~ do nn=jul_cumdays(n),jul_cumdays(n+1)-1
			!~ tm_seasonal_n1(nn)=count
			!~ tm_seasonal_n2(nn)=count+1
		!~ end do
		!~ count=count+1
	!~ end do
	!~ where(tm_seasonal_n2==13) tm_seasonal_n2=1
	
	
	!~ sum_val=0.0
	!~ do n=1,12
		!~ do nn=jul_cumdays(n),jul_cumdays(n+1)-1
			!~ tm_seasonal_rscale(nn)=sum_val
			!~ sum_val=sum_val+(1.0/jul_days(n))
		!~ end do
		!~ sum_val=0.0
	!~ end do
	!~ tm_seasonal_scale=1.0-tm_seasonal_rscale


!~ else

	!~ count=1
	!~ do n=1,tm_n_dt

		!~ tm_seasonal_n1(n)=count
		!~ tm_seasonal_n2(n)=count
		
		!~ tm_seasonal_scale(n)=1.0
		!~ tm_seasonal_rscale(n)=1.0-tm_seasonal_scale(n)
	
	!~ end do
	
!~ end if


end subroutine calc_seasonal_scaling
! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_TM_data()

! load data from files

call load_TM_netcdf('../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aexp_filename),Aexp)
call load_TM_netcdf('../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aimp_filename),Aimp)
call load_TM_netcdf(tm_Aremin_filename,Aremin)
!call load_seaice()
call load_TM_bgc_data()
call load_PO4_restore()
call load_volume()

if(bg_PO4restore_select)then
	call load_PO4_restore()
	else 
	call load_PO4_uptake()
end if

if(gen_restart_select)then
	call load_netcdf_restart()
end if


end subroutine load_TM_data

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_tm(A_filename,A)

character(len=100),intent(inout)::A_filename
type(sparse)::A

integer::n

open(unit=20,file='../data/'//trim(tm_data_fileloc)//'/'//trim(A_filename)//'_val') ! TM values
read(unit=20,fmt="(e22.15)") A%val
close(unit=20)

open(unit=20,file='../data/'//trim(tm_data_fileloc)//'/'//trim(A_filename)//'_col') ! TM column index
read(unit=20,fmt="(I8)") A%col
close(unit=20)
 
open(unit=20,file='../data/'//trim(tm_data_fileloc)//'/'//trim(A_filename)//'_rowptr') ! TM compressed row index
read(unit=20,fmt="(I8)") A%row
close(unit=20)


end subroutine
! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

SUBROUTINE load_TM_netcdf(dum_filename,dum_A)

character(len=100),intent(in)::dum_filename
type(sparse)::dum_A

integer::loc_ncid,loc_varid,status
character(len=100)::loc_lname

!print*,'loading TM data from:','OMFG/data/'//trim(tm_data_fileloc)//'/'//trim(dum_filename)
!print*,'loading TM data from:',dum_filename

! open netcdf file
!status=nf90_open('../data/'//trim(tm_data_fileloc)//'/'//trim(dum_filename) , nf90_nowrite,loc_ncid)
status=nf90_open('../data/'//'/'//dum_filename , nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!print*,dum_filename
! matrix values
status=nf90_inq_varid(loc_ncid,'A_val',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,dum_A%val)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!print*,dum_A%val(1:20,1)
!print*,dum_A%val(1:10,12)


! matrix column indices
status=nf90_inq_varid(loc_ncid,'A_col',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,dum_A%col)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!print*,dum_A%col(1:20)

! matrix row pointer indices
STATUS=nf90_inq_varid(loc_ncid,'A_row',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,dum_A%row)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!print*,dum_A%row(1:20)

! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))


end subroutine load_TM_netcdf
! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine load_seaice()

! local variables
integer::n,status,loc_varid,loc_ncid

!print*,'Reading in seaice data from:','../data/'//trim(tm_data_fileloc)//'/'//trim(tm_seaice_filename)

! open netcdf file
status=nf90_open('../data/'//trim(tm_data_fileloc)//'/'//trim(tm_seaice_filename) , nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! matrix values
status=nf90_inq_varid(loc_ncid,'Fice',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,tm_seaice_frac,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!print*,tm_seaice_frac(1:10,1)
!print*,tm_seaice_frac(4448:4448+10,1)

! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

tm_seaice_frac=1.0-tm_seaice_frac

end subroutine 

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine load_TM_bgc_data()

! local variables
integer::n,status,loc_varid,loc_ncid

!print*,'Reading in seaice data from:','../data/'//trim(tm_data_fileloc)//'/'//trim(tm_bgc_data_filename)

! open netcdf file
status=nf90_open('../data/'//trim(tm_data_fileloc)//'/'//trim(tm_bgc_data_filename) , nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! windspeed
status=nf90_inq_varid(loc_ncid,'windspeed',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,tm_windspeed,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! seaice fraction
status=nf90_inq_varid(loc_ncid,'Fice',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,tm_seaice_frac,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! T
status=nf90_inq_varid(loc_ncid,'Tbc',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,tm_T,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! S
status=nf90_inq_varid(loc_ncid,'Sbc',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,tm_S,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! silica
status=nf90_inq_varid(loc_ncid,'silica',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,tm_silica,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))


! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

!print*,tm_seaice_frac(4453:4460,1)
!print*,tm_windspeed(1:10,1)
!print*,tm_T(1:10,1)
!print*,tm_S(1:10,1)
!print*,tm_silica(1:10,1)

tm_seaice_frac=1.0-tm_seaice_frac

end subroutine 
! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_PO4_restore()

! local variables
integer::n,status,loc_varid,loc_ncid

!print*,'Reading in PO4restore data from:','../data/'//trim(tm_data_fileloc)//'/'//trim(tm_PO4restore_filename)

! open netcdf file
status=nf90_open('../data/'//trim(tm_data_fileloc)//'/'//trim(tm_PO4restore_filename), nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! matrix values
status=nf90_inq_varid(loc_ncid,'PO4_Obs',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,bg_PO4_obs,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!print*,bg_PO4_obs(1:10,1)

! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

end subroutine 

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine load_PO4_uptake()

! local variables
integer::n,status,loc_varid,loc_ncid

!print*,'Reading in PO4uptake data from:','../data/'//trim(tm_data_fileloc)//'/'//trim(tm_PO4uptake_filename)

! open netcdf file
status=nf90_open('../data/'//trim(tm_data_fileloc)//'/'//trim(tm_PO4uptake_filename), nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! matrix values
status=nf90_inq_varid(loc_ncid,'PO4_Uptake',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,bg_PO4_uptake,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!print*,bg_PO4_uptake(1:10,1)

! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

end subroutine 

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine load_data_saving()

! local variables
integer::var_count,n
integer::iostat
logical::exist_file
real::tmp

! open, read file for dimension, and allocate
INQUIRE(FILE='../data/'//trim(gen_save_timeseries_file), EXIST=exist_file)
if(exist_file.eqv..false.) print*,'Timeseries input file does not exist'

open(UNIT=10,FILE='../data/'//trim(gen_save_timeseries_file))
var_count=0
DO
	READ(UNIT=10,IOSTAT=iostat,FMT=*)  tmp
	IF (iostat < 0) THEN
       exit
    ELSE
       var_count=var_count+1
    END IF
END DO
close(UNIT=10)
allocate(tm_timeseries(var_count))

! read in timeseries save points
open(UNIT=10,FILE='../data/'//trim(gen_save_timeseries_file))
read(unit=10,iostat=iostat,fmt='(F10.1)') tm_timeseries
close(unit=10)



end subroutine 

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_volume()

! local variables
integer::n,status,loc_varid,loc_ncid

!~ print*,'Reading in volume data'

!~ open(unit=20,file='../data/'//trim(tm_data_fileloc)//'/'//trim(tm_vol_filename)) ! TM values
!~ read(unit=20,fmt="(e22.15)") tm_vol
!~ close(unit=20)
!~ print*,'Reading in volume data - completed'

!print*,'Reading in volume data from:','../data/'//trim(tm_data_fileloc)//'/'//trim(tm_vol_filename)

! open netcdf file
status=nf90_open('../data/'//trim(tm_data_fileloc)//'/'//trim(tm_vol_filename) , nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! matrix values
status=nf90_inq_varid(loc_ncid,'vol',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(loc_ncid,loc_varid,tm_vol)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!print*,tm_vol(1:10,1)

! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))






end subroutine 

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_netcdf_restart()

integer::ncid,varid,status

!print*,'Reading in restart data from:','../output/'//trim(gen_restart_filename)//'.nc'

! open netcdf file
status=nf90_open('../output/'//trim(gen_restart_filename)//'.nc',nf90_nowrite,ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! PO4
status=nf90_inq_varid(ncid,'PO4',varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(ncid,varid,tracers_1(:,ioPO4),start=(/ 1, 12 /),count=(/ tm_nbox , 1 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!print*,tracers_1(1:10,iPO4)

! DOP
status=nf90_inq_varid(ncid,'DOP',varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_get_var(ncid,varid,tracers_1(:,ioDOP),start=(/ 1, 12 /),count=(/ tm_nbox , 1 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!print*,tracers_1(1:10,iDOP)

end subroutine load_netcdf_restart

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

FUNCTION amul(A,Vector)
! output
REAL,dimension(tm_nbox,gen_n_tracers)::amul
! dummy
type(sparse),intent(in)::A
REAL,INTENT(in),dimension(tm_nbox,gen_n_tracers)::Vector
! local
integer::n,nn,i
real::sum_val
integer,dimension(2)::vector_size
real,dimension(A%nnz)::val_tmp

vector_size=shape(Vector)

! tmp code
!print*,dt_count
!print*,tm_seasonal_scale(dt_count),tm_seasonal_rscale(dt_count)
!print*,tm_seasonal_n1(dt_count),tm_seasonal_n2(dt_count)
val_tmp=(tm_seasonal_scale(dt_count)*A%val(:,tm_seasonal_n1(dt_count)))&
+&
((tm_seasonal_rscale(dt_count))*A%val(:,tm_seasonal_n2(dt_count)))
! tmp code
!PRINT*,val_tmp(1:10)



! sum_val=sum_val+A%val(nn,1)*Vector(A%col(nn),i)
DO i=1,vector_size(2)
do n=1,tm_nbox
	sum_val=0.0
	
	!print*,n
	!PRINT*,A%row(n)
	!PRINT*,A%row(n+1)-1
	!PRINT*,val_tmp(A%row(n))
	!PRINT*,Vector(A%col(A%row(n)),i)
	
	do nn=A%row(n),A%row(n+1)-1
		sum_val=sum_val+val_tmp(nn)*Vector(A%col(nn),i)
	end do
	amul(n,i)=sum_val
end do
end do

!sum_val=sum_val+&
!((A%val(nn)*Vector(A%col(nn),i)*tm_seasonal_scale(1))+&
!(A%val(nn)*Vector(A%col(nn),i)*(1.0-tm_seasonal_scale(1))))

end FUNCTION

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

FUNCTION amul_remin(A,Vector)
! output
REAL,dimension(tm_nbox)::amul_remin
! dummy
type(sparse),intent(in)::A
REAL,INTENT(in),dimension(tm_nbox)::Vector
! local
integer::n,nn,i
real::sum_val


do n=1,tm_nbox
	sum_val=0.0
	
	!print*,n
	!PRINT*,A%row(n)
	!PRINT*,A%row(n+1)-1
	!PRINT*,A%val(A%row(n),1)
	!PRINT*,Vector(A%col(A%row(n)))
	
	do nn=A%row(n),A%row(n+1)-1
		sum_val=sum_val+A%val(nn,1)*Vector(A%col(nn))
	end do
	amul_remin(n)=sum_val
end do


end FUNCTION

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine print_to_screen(dum_t,dum_extra)

integer::dum_t
real::dum_extra

print*,'year',dum_t/tm_n_dt,sum(tracers(:,ioPO4)*tm_vol),sum(tracers(:,ioDOP)*tm_vol),sum(tracers(:,ioDIC)*tm_vol), &
 & sum(tracers(:,ioALK)*tm_vol),ATM(iaCO2)*1.0e6,dum_extra



end subroutine print_to_screen

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine write_output_netcdf()

integer::ncid,status,m_dimid,n_dimid,PO4id,DOPid,EXPORTid,n_season
integer::dimids(gen_n_tracers)

! create file
! status=nf90_create('../output/'//trim(gen_config_filename)//'.nc',nf90_clobber,ncid)
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! define dimensions
! status=nf90_def_dim(ncid,'tm_nbox',tm_nbox,m_dimid)
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! status=nf90_def_dim(ncid,'time',n_seasonal,n_dimid)
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! dimids=(/m_dimid,n_dimid/)
! 
! define variable
! status=nf90_def_var(ncid,'PO4',nf90_float,dimids,PO4id)
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! status=nf90_def_var(ncid,'DOP',nf90_float,dimids,DOPid)
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! status=nf90_def_var(ncid,'EXPORT',nf90_float,dimids,EXPORTid)
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! end definition
! status=nf90_enddef(ncid)
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! write data
! status=nf90_put_var(ncid,PO4id,tracers_PO4_int(:,:))
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! status=nf90_put_var(ncid,DOPid,tracers_DOP_int(:,:))
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! status=nf90_put_var(ncid,EXPORTid,EXPORT_int(:,:))
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! status=nf90_close(ncid)
! if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
! 
! print*,'Output written to: ','../output/'//trim(gen_config_filename)//'.nc'

end subroutine write_output_netcdf

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine tm_vars_at_dt()

! linearly interpolate seaice, windstress, T, S at timestep
! tm_seasonal_scale, tm_seasonal_rscale, tm_seasonal_n1, tm_seasonal, n2: from calc_seasonal_scaling
! dt_count: fml.f90

seaice_dt=(tm_seasonal_scale(dt_count)*tm_seaice_frac(:,tm_seasonal_n1(dt_count))) &
+ &
((tm_seasonal_rscale(dt_count))*tm_seaice_frac(:,tm_seasonal_n2(dt_count)))

wind_dt=(tm_seasonal_scale(dt_count)*tm_windspeed(:,tm_seasonal_n1(dt_count))) &
+ &
((tm_seasonal_rscale(dt_count))*tm_windspeed(:,tm_seasonal_n2(dt_count)))

T_dt=(tm_seasonal_scale(dt_count)*tm_T(:,tm_seasonal_n1(dt_count))) &
+ &
((tm_seasonal_rscale(dt_count))*tm_T(:,tm_seasonal_n2(dt_count)))

S_dt=(tm_seasonal_scale(dt_count)*tm_S(:,tm_seasonal_n1(dt_count))) &
+ &
((tm_seasonal_rscale(dt_count))*tm_S(:,tm_seasonal_n2(dt_count)))

silica_dt=(tm_seasonal_scale(dt_count)*tm_silica(:,tm_seasonal_n1(dt_count))) &
+ &
((tm_seasonal_rscale(dt_count))*tm_silica(:,tm_seasonal_n2(dt_count)))

! convert wind_dt to correct units for gas exchange
! not pre-calculated due to non-linear terms (u^2)
!!!! *** windspeed is m/s? so adjust this line of code *** !!!!
wind_dt=(wind_dt**2)*bg_gastransfer_a*seaice_dt*conv_sec_yr


end subroutine tm_vars_at_dt

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine initialise_output()

logical::exist_dir

! make output directory (if doesn't already exist)
INQUIRE(FILE='../output/'//trim(gen_config_filename), EXIST=exist_dir)
if(exist_dir.eqv..false.)then
	call system('mkdir ../output/'//trim(gen_config_filename))
end if

! initialise timeseries files
call initialise_timeseries_output()

end subroutine initialise_output

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine initialise_timeseries_output()

! atm CO2
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_atm_CO2.dat',status='replace')
write(unit=10,fmt='(A44)') '% / year / atmospheric CO2 (ppmv)' 
close(unit=10)

! PO4
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_PO4.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / PO4 inventory (mol) / [PO4] global (mmol m-3)' 
close(unit=10)

! DOP
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_DOP.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / DOP inventory (mol) / [DOP] global (mmol m-3)' 
close(unit=10)

! DIC
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_DIC.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / DIC inventory (mol) / [DIC] global (mmol m-3)' 
close(unit=10)

! ALK
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_ALK.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / ALK inventory (mol) / [ALK] global (mmol m-3)' 
close(unit=10)



end subroutine initialise_timeseries_output

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine integrate_output(loc_t,save_count)

integer,intent(inOUT)::save_count
integer,intent(in)::loc_t
real::scalar



!tracers_PO4_int(:,save_count)=tracers_PO4_int(:,save_count)+tracers(:,ioPO4)*tm_dt*n_seasonal
!tracers_DOP_int(:,save_count)=tracers_DOP_int(:,save_count)+tracers(:,ioDOP)*tm_dt*n_seasonal
!EXPORT_int(:,save_count)=EXPORT_int(:,save_count)+export(:)*tm_dt*n_seasonal

!if(mod(dt_count,tm_n_dt/n_seasonal).eq.0 .and. tm_seasonal)then
!	save_count=save_count+1
!end if

if(loc_t>=tm_timeseries(timeseries_count)*96-47 .and. loc_t<=tm_timeseries(timeseries_count)*96+48)then

	! integrate 
	scalar=bg_dt*tm_save_intra_freq
	tracers_int(:,:)=tracers_int(:,:)+tracers(:,:)*scalar
	ATM_int(:)=ATM_int(:)+ATM(:)*scalar
	t_int=t_int+((loc_t-1.0)/(1.0/bg_dt))*scalar

	! write when reached end of time period
	! and reset integrating arrays
	if(loc_t==tm_timeseries(timeseries_count)*96+48)then
		call write_timeseries_output()
		timeseries_count=timeseries_count+1
		ATM_int(:)=0.0
		tracers_int(:,:)=0.0
		t_int=0.0
	endif
		
	
endif

end subroutine integrate_output

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine write_timeseries_output()

real::rvol_tot

rvol_tot=1.0/sum(tm_vol)

! CO2
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_atm_CO2.dat',position='append')
write(unit=10,fmt='(f12.1,f12.6)') &
t_int , &
ATM_int(iaCO2)*1.0e6
close(unit=10)

!PO4
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_PO4.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6)') &
t_int , &
sum(tracers_int(:,ioPO4)*tm_vol) , &
sum(tracers_int(:,ioPO4)*tm_vol)*rvol_tot*1.0e3
close(unit=10)

!DOP
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_DOP.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6)') &
t_int , &
sum(tracers_int(:,ioDOP)*tm_vol) , &
sum(tracers_int(:,ioDOP)*tm_vol)*rvol_tot*1.0e3
close(unit=10)

!DIC
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_DIC.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6)') &
t_int , &
sum(tracers_int(:,ioDIC)*tm_vol) , &
sum(tracers_int(:,ioDIC)*tm_vol)*rvol_tot*1.0e3
close(unit=10)

!ALK
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_ALK.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6)') &
t_int , &
sum(tracers_int(:,ioALK)*tm_vol) , &
sum(tracers_int(:,ioALK)*tm_vol)*rvol_tot*1.0e3
close(unit=10)


end subroutine write_timeseries_output

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


end module tm_module