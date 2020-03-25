module io_module
! ---------------------------------------------------------------------------------------!
!! Subroutines to load and write data
! ---------------------------------------------------------------------------------------!

use fml_lib
use netcdf

implicit none
contains

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_TM_data()
! ---------------------------------------------------------------------------------------!
!! OLD: Load transport matrix data from file
! ---------------------------------------------------------------------------------------!

call load_TM_netcdf('data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aexp_filename),Aexp)
call load_TM_netcdf('data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aimp_filename),Aimp)
!call load_TM_netcdf(tm_Aremin_filename,Aremin)
!Aremin%val=Aremin%val_n(:,1) ! tmp

call load_TM_grid_data()
call load_TM_bgc_data()

select case(trim(bg_uptake_function))
case('restore')
	call load_PO4_restore()
case('fixed')
	call load_PO4_uptake()
end select

if(bg_martin_remin_spatial)then
	call load_Martin_b_spatial()
else
	bg_martin_b(:)=bg_martin_remin_b ! set global value to whole grid
endif

end subroutine load_TM_data

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

SUBROUTINE load_TM_netcdf(dum_filename,dum_A)
! ---------------------------------------------------------------------------------------!
!! Loads transport matrix data from netCDF files
! ---------------------------------------------------------------------------------------!

character(len=*)::dum_filename
!! netCDF filename
type(sparse)::dum_A
!! sparse matrix arrays

integer::loc_ncid,loc_varid,status
character(len=100)::loc_lname

print*,'loading TM data from: ',trim(dum_filename)

! open netcdf file
status=nf90_open(trim(dum_filename),nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),dum_filename

! matrix values
status=nf90_inq_varid(loc_ncid,'A_val',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'A_val'

status=nf90_get_var(loc_ncid,loc_varid,dum_A%val_n)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'A_val'
!print*,dum_A%val_n(1:10,1)

! matrix column indices
status=nf90_inq_varid(loc_ncid,'A_col',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'A_col'

status=nf90_get_var(loc_ncid,loc_varid,dum_A%col)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'A_col'
!print*,dum_A%col(1:20)

! matrix row pointer indices
STATUS=nf90_inq_varid(loc_ncid,'A_row',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'A_row'

status=nf90_get_var(loc_ncid,loc_varid,dum_A%row)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'A_row'
!print*,dum_A%row(1:20)

! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))


end subroutine load_TM_netcdf

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

SUBROUTINE load_TM_metadata(dum_filename,dum_A)
! ---------------------------------------------------------------------------------------!
!! Loads transport matrix metadata from netCDF file
! ---------------------------------------------------------------------------------------!

character(len=*)::dum_filename
!! netCDF filename
type(sparse)::dum_A
!! sparse matrix arrays

integer::loc_ncid,loc_varid,status
character(len=100)::loc_lname

! open netcdf file
status=nf90_open(trim(dum_filename),nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),dum_filename

! nonzeros
status=nf90_inq_varid(loc_ncid,'nnz',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'nnz'

status=nf90_get_var(loc_ncid,loc_varid,dum_A%nnz)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'nnz'
!print*,dum_A%nnz

! nonzeros
status=nf90_inq_varid(loc_ncid,'nb',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'nb'

status=nf90_get_var(loc_ncid,loc_varid,dum_A%nb)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'nb'
!print*,dum_A%nb

! nonzeros
status=nf90_inq_varid(loc_ncid,'n_season',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'n_time'

status=nf90_get_var(loc_ncid,loc_varid,dum_A%n_time)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'n_time'
!print*,dum_A%n_time

! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))


end subroutine load_TM_metadata

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_TM_grid_data()
! ---------------------------------------------------------------------------------------!
!! Loads transport matrix grid data from netCDF file
! ---------------------------------------------------------------------------------------!

! local variables
integer::n,status,loc_varid,loc_ncid

! open netcdf file
status=nf90_open('data/'//trim(tm_data_fileloc)//'/'//trim(tm_grid_data_filename) , nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'->'//tm_grid_data_filename

! volume
status=nf90_inq_varid(loc_ncid,'volume',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> volume'

status=nf90_get_var(loc_ncid,loc_varid,tm_vol)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> volume'

! area
status=nf90_inq_varid(loc_ncid,'area',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> area'

status=nf90_get_var(loc_ncid,loc_varid,tm_area)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> area'

! longitude
status=nf90_inq_varid(loc_ncid,'Longitude',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> longitude'

status=nf90_get_var(loc_ncid,loc_varid,tm_lon)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'->longitude'

! latitude
status=nf90_inq_varid(loc_ncid,'Latitude',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> latitude'

status=nf90_get_var(loc_ncid,loc_varid,tm_lat)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> latitude'

! depth
status=nf90_inq_varid(loc_ncid,'Depth',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> depth'

status=nf90_get_var(loc_ncid,loc_varid,tm_depth)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> depth'

! bottom depth
status=nf90_inq_varid(loc_ncid,'Depth_Btm',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> depth_btm'

status=nf90_get_var(loc_ncid,loc_varid,tm_depth_btm)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> depth_btm'

! i
status=nf90_inq_varid(loc_ncid,'i',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> i'

status=nf90_get_var(loc_ncid,loc_varid,tm_i)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> i'

! j
status=nf90_inq_varid(loc_ncid,'j',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> j'

status=nf90_get_var(loc_ncid,loc_varid,tm_j)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> j'

! k
status=nf90_inq_varid(loc_ncid,'k',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> k'

status=nf90_get_var(loc_ncid,loc_varid,tm_k)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> k'


! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))


end subroutine load_TM_grid_data

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine load_TM_bgc_data()
! ---------------------------------------------------------------------------------------!
!! Loads transport matrix biogeochemistry data from netCDF file
! ---------------------------------------------------------------------------------------!

! local variables
integer::n,status,loc_varid,loc_ncid

! open netcdf file
status=nf90_open('data/'//trim(tm_data_fileloc)//'/'//trim(tm_bgc_data_filename) , nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'->'//tm_bgc_data_filename

! windspeed
status=nf90_inq_varid(loc_ncid,'windspeed',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> windspeed'

status=nf90_get_var(loc_ncid,loc_varid,tm_windspeed,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> windspeed'

! seaice fraction
status=nf90_inq_varid(loc_ncid,'Fice',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> Fice'

status=nf90_get_var(loc_ncid,loc_varid,tm_seaice_frac,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> Fice'

! T
status=nf90_inq_varid(loc_ncid,'Tbc',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> T'

status=nf90_get_var(loc_ncid,loc_varid,tm_T,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> T'

! S
status=nf90_inq_varid(loc_ncid,'Sbc',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> S'

status=nf90_get_var(loc_ncid,loc_varid,tm_S,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> S'

! silica
status=nf90_inq_varid(loc_ncid,'silica',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> Si'

status=nf90_get_var(loc_ncid,loc_varid,tm_silica,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))//'-> Si'


! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! convert to fraction of grid-box not covered, so as to calculate once only
tm_seaice_frac=1.0-tm_seaice_frac

end subroutine

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_PO4_restore()
! ---------------------------------------------------------------------------------------!
!! Loads [PO4] observations from netCDF file for nutrient restoring subroutine
! ---------------------------------------------------------------------------------------!

! local variables
integer::n,status,loc_varid,loc_ncid

print*,'loading in PO4restore data from: ','/data/'//trim(tm_PO4restore_filename)

! open netcdf file
status=nf90_open('data/'//trim(tm_PO4restore_filename), nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),tm_PO4restore_filename

! matrix values
status=nf90_inq_varid(loc_ncid,'PO4_Obs',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'PO4_Obs'

status=nf90_get_var(loc_ncid,loc_varid,bg_PO4_obs,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'PO4_Obs'
!print*,bg_PO4_obs(1:10,1)

! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

end subroutine

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine load_Martin_b_spatial()
! ---------------------------------------------------------------------------------------!
!! Loads spatially explicit Martin curve "b" data from netCDF
! ---------------------------------------------------------------------------------------!

! local variables
integer::n,status,loc_varid,loc_ncid

print*,'loading in spatially varying b data from: ','data/'//trim(bg_martin_b_input_filename)

! open netcdf file
status=nf90_open('data/'//'/'//trim(bg_martin_b_input_filename), nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),bg_martin_b_input_filename

! matrix values
status=nf90_inq_varid(loc_ncid,'Martin_b',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'Spatial_Martin_b'

status=nf90_get_var(loc_ncid,loc_varid,bg_martin_b)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'Spatial_Martin_b'
!print*,bg_PO4_uptake(1:10,1)

! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

end subroutine load_Martin_b_spatial

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_data_saving()
! ---------------------------------------------------------------------------------------!
!! Loads data specifying timeseries and timeslice save points
! ---------------------------------------------------------------------------------------!

! local variables
integer::var_count,n
integer::iostat
logical::exist_file
real::tmp

! ******************** timeseries ********************
! open, read file for dimension, and allocate
INQUIRE(FILE='data/'//trim(gen_save_timeseries_file), EXIST=exist_file)
if(exist_file.eqv..false.) print*,'Timeseries input file does not exist'

open(UNIT=10,FILE='data/'//trim(gen_save_timeseries_file))
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
open(UNIT=10,FILE='data/'//trim(gen_save_timeseries_file))
read(unit=10,iostat=iostat,fmt='(F10.1)') tm_timeseries
close(unit=10)

! ******************** timeslice ********************
! open, read file for dimension, and allocate
INQUIRE(FILE='data/'//trim(gen_save_timeslice_file), EXIST=exist_file)
if(exist_file.eqv..false.) print*,'Timeseries input file does not exist'

open(UNIT=10,FILE='data/'//trim(gen_save_timeslice_file))
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
allocate(tm_timeslice(var_count))

! read in timeslice save points
open(UNIT=10,FILE='data/'//trim(gen_save_timeslice_file))
read(unit=10,iostat=iostat,fmt='(F10.1)') tm_timeslice
close(unit=10)

end subroutine load_data_saving

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine initialise_output()
! ---------------------------------------------------------------------------------------!
!! Initialises timeslice and timeseries output files
!!
!! - creates output and restart directories if needed
! ---------------------------------------------------------------------------------------!

logical::exist_dir

! make output directory (if doesn't already exist)
INQUIRE(FILE='output/'//trim(gen_config_filename), EXIST=exist_dir)
if(exist_dir.eqv..false.)then
	call system('mkdir output/'//trim(gen_config_filename))
end if

! make restart directory (if doesn't already exist)
INQUIRE(FILE='output/'//trim(gen_config_filename)//'/restart', EXIST=exist_dir)
if(exist_dir.eqv..false.)then
	call system('mkdir output/'//trim(gen_config_filename)//'/restart')
end if

! initialise timeseries and timeslice output files
call initialise_timeseries_output()
call initialise_timeslice_output()

! write out copy of namelist file
open(unit=20,file='output/'//trim(gen_config_filename)//'/parameter_namelist.txt',status='replace')
write( UNIT=20, NML=fml_namelist)
close(unit=20)

end subroutine initialise_output

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine initialise_timeseries_output()

! ---------------------------------------------------------------------------------------!
!! Initialise timeseries output files
! ---------------------------------------------------------------------------------------!

! atm CO2
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_atm_CO2.dat',status='replace')
write(unit=10,fmt='(A44)') '% / year / atmospheric CO2 (mol) / atmospheric CO2 (ppmv)'
close(unit=10)

! PO4
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_ocn_PO4.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / PO4 inventory (mol) / [PO4] global (mmol m-3) / [PO4] surface (mmol m-3)'
close(unit=10)

! DOP
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_ocn_DOP.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / DOP inventory (mol) / [DOP] global (mmol m-3)'
close(unit=10)

! DIC
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_ocn_DIC.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / DIC inventory (mol) / [DIC] global (mmol m-3) / [DIC] surface (mmol m-3)'
close(unit=10)

! ALK
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_ocn_ALK.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / ALK inventory (mol) / [ALK] global (mmol m-3)'
close(unit=10)

! gas exchange
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_airsea_CO2.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / global mean flux (mol yr-1) / global mean flux (mol m-2 yr-1)'
close(unit=10)

! POM EXPORT
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_POM_export.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / global mean flux (mol yr-1) / global mean flux (mol m-2 yr-1)'
close(unit=10)

! CaCO3 EXPORT
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_CaCO3_export.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / global mean flux (mol yr-1) / global mean flux (mol m-2 yr-1)'
close(unit=10)



end subroutine initialise_timeseries_output

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine write_timeseries_output()
! ---------------------------------------------------------------------------------------!
!! Writes timeseries output to file
! ---------------------------------------------------------------------------------------!

real::rvol_tot,rvol_sur_tot,rvol_deep_tot

rvol_tot=1.0/sum(tm_vol)
rvol_sur_tot=1.0/sum(tm_vol(1:n_surface_boxes))
rvol_deep_tot=1.0/sum(tm_vol(41153:52749))

! CO2
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_atm_CO2.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6)') &
t_int , &
ATM_int(iaCO2)*ATM_vol, &
ATM_int(iaCO2)*1.0e6
close(unit=10)

!PO4
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_ocn_PO4.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6,f12.6,f12.6)') &
t_int , &
sum(tracers_int(:,ioPO4)*tm_vol) , &
sum(tracers_int(:,ioPO4)*tm_vol)*rvol_tot*1.0e3 , &
sum(tracers_int(1:n_surface_boxes,ioPO4)*tm_vol(1:n_surface_boxes))*rvol_sur_tot*1.0e3 , &
sum(tracers_int(41153:52749,ioPO4)*tm_vol(41153:52749))*rvol_deep_tot*1.0e3
close(unit=10)

!DOP
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_ocn_DOP.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6)') &
t_int , &
sum(tracers_int(:,ioDOP)*tm_vol) , &
sum(tracers_int(:,ioDOP)*tm_vol)*rvol_tot*1.0e3
close(unit=10)

!DIC
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_ocn_DIC.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6,f12.6,f12.6)') &
t_int , &
sum(tracers_int(:,ioDIC)*tm_vol) , &
sum(tracers_int(:,ioDIC)*tm_vol)*rvol_tot*1.0e3 , &
sum(tracers_int(1:n_surface_boxes,ioDIC)*tm_vol(1:n_surface_boxes))*rvol_sur_tot*1.0e3 , &
sum(tracers_int(41153:52749,ioDIC)*tm_vol(41153:52749))*rvol_deep_tot*1.0e3
close(unit=10)

!ALK
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_ocn_ALK.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6)') &
t_int , &
sum(tracers_int(:,ioALK)*tm_vol) , &
sum(tracers_int(:,ioALK)*tm_vol)*rvol_tot*1.0e3
close(unit=10)

!gas exchange
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_airsea_CO2.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,e20.12)') &
t_int , &
sum(diag_int(:,2)), &
sum(diag_int(:,2))/sum(tm_area(1:n_surface_boxes)) ! check this
close(unit=10)

!POM export
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_POM_export.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12)') &
t_int , &
sum(diag_int(:,3))
close(unit=10)

!POM export
open(unit=10,file='output/'//trim(gen_config_filename)//'/timeseries_CaCO3_export.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12)') &
t_int , &
sum(diag_int(:,6))
close(unit=10)




end subroutine write_timeseries_output

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine initialise_timeslice_output()
! ---------------------------------------------------------------------------------------!
!! Initialises netCDF output
! ---------------------------------------------------------------------------------------!

	integer::ncid,status,m_dimid,n_dimid,PO4id,DOPid,EXPORTid,n_season,DICid,ALKid,airseaCO2id,popreminid,po4uptakeid
	integer::lon_dimid,lat_dimid,depth_dimid
	!integer::dimids(2)
	integer::dimids(4)

!create file
status=nf90_create('output/'//trim(gen_config_filename)//'/fields_netcdf.nc',nf90_clobber,ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

!define dimensions
!status=nf90_def_dim(ncid,'tm_nbox',tm_nbox,m_dimid)
!if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_dim(ncid,'lon',NX,lon_dimid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_dim(ncid,'lat',NY,lat_dimid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_dim(ncid,'depth',NZ,depth_dimid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_dim(ncid,'time',nf90_unlimited,n_dimid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

dimids=(/lon_dimid,lat_dimid,depth_dimid,n_dimid/)
!dimids=(/m_dimid,n_dimid/)

!define variable
status=nf90_def_var(ncid,'PO4',nf90_float,dimids,PO4id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_var(ncid,'DOP',nf90_float,dimids,DOPid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_var(ncid,'DIC',nf90_float,dimids,DICid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_var(ncid,'ALK',nf90_float,dimids,ALKid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_var(ncid,'airsea_CO2_flux',nf90_float,dimids,airseaCO2id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_var(ncid,'PO4_uptake',nf90_float,dimids,po4uptakeid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_var(ncid,'POP_remin',nf90_float,dimids,popreminid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_var(ncid,'CaCO3_export',nf90_float,dimids,popreminid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_var(ncid,'CaCO3_remin',nf90_float,dimids,popreminid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

!status=nf90_def_var(ncid,'EXPORT',nf90_float,dimids,EXPORTid)
!if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

!end definition
status=nf90_enddef(ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! close file
status=nf90_close(ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))



end subroutine initialise_timeslice_output

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine write_output_netcdf()
! ---------------------------------------------------------------------------------------!
!! writes netCDF output
!!
!! - currently just for the last timestep!!
! ---------------------------------------------------------------------------------------!

integer::ncid,status,var_id
real,dimension(NX,NY,NZ)::field

!open file
status=nf90_open('output/'//trim(gen_config_filename)//'/fields_netcdf.nc',nf90_write,ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! inquire ids and write data

! PO4
status=nf90_inq_varid(ncid,'PO4',var_id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
field=v2f(tracers_int(:,ioPO4))
status=nf90_put_var(ncid,var_id,field,start=(/ 1, 1, 1, timeslice_count /),count=(/NX,NY,NZ,1/))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! ! DOP
status=nf90_inq_varid(ncid,'DOP',var_id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
field=v2f(tracers_int(:,ioDOP))
status=nf90_put_var(ncid,var_id,field,start=(/ 1, 1, 1, timeslice_count /),count=(/NX,NY,NZ,1/))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!
! ! DIC
status=nf90_inq_varid(ncid,'DIC',var_id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
field=v2f(tracers_int(:,ioDIC))
status=nf90_put_var(ncid,var_id,field,start=(/ 1, 1, 1, timeslice_count /),count=(/NX,NY,NZ,1/))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!
! ! ALK
status=nf90_inq_varid(ncid,'ALK',var_id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
field=v2f(tracers_int(:,ioALK))
status=nf90_put_var(ncid,var_id,field,start=(/ 1, 1, 1, timeslice_count /),count=(/NX,NY,NZ,1/))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
!
! ! air-sea CO2 flux
status=nf90_inq_varid(ncid,'airsea_CO2_flux',var_id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
field=v2f(diag_int(:,1))
status=nf90_put_var(ncid,var_id,field,start=(/ 1, 1, 1, timeslice_count /),count=(/NX,NY,NZ,1/))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! ! PO4 uptake
status=nf90_inq_varid(ncid,'PO4_uptake',var_id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
field=v2f(diag_int(:,3))
status=nf90_put_var(ncid,var_id,field,start=(/ 1, 1, 1, timeslice_count /),count=(/NX,NY,NZ,1/))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! ! POP remin
status=nf90_inq_varid(ncid,'POP_remin',var_id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
field=v2f(diag_int(:,4))
status=nf90_put_var(ncid,var_id,field,start=(/ 1, 1, 1, timeslice_count /),count=(/NX,NY,NZ,1/))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! ! POP remin
status=nf90_inq_varid(ncid,'CaCO3_export',var_id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
field=v2f(diag_int(:,6))
status=nf90_put_var(ncid,var_id,field,start=(/ 1, 1, 1, timeslice_count /),count=(/NX,NY,NZ,1/))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! ! POP remin
status=nf90_inq_varid(ncid,'CaCO3_remin',var_id)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))
field=v2f(diag_int(:,5))
status=nf90_put_var(ncid,var_id,field,start=(/ 1, 1, 1, timeslice_count /),count=(/NX,NY,NZ,1/))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

! close file
status=nf90_close(ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

!print*,'Output written to: ','../output/'//trim(gen_config_filename)//'.nc'

end subroutine write_output_netcdf

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine write_restart()
! ---------------------------------------------------------------------------------------!
!! Writes state variables to binary for restart
! ---------------------------------------------------------------------------------------!

integer::ios

open(unit=10,form='unformatted',status='replace',action='write',iostat=ios, &
file='output/'//trim(gen_config_filename)//'/restart/restart.bin')
write(10) tracers_1 , ATM
close(10,iostat=ios)

print*,'Restart saved to: '//'output/'//trim(gen_config_filename)//'/restart/restart.bin'

end subroutine write_restart

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_restart()
! ---------------------------------------------------------------------------------------!
!! Loads state variables from binary from restart experiment
! ---------------------------------------------------------------------------------------!

integer::ios

open(unit=10,status='old',form='unformatted',action='read',IOSTAT=ios, &
file='output/'//trim(gen_restart_filename)//'/restart/restart.bin')
read(10,iostat=ios) tracers_1 , ATM
close(10,iostat=ios)

end subroutine load_restart

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine write_PO4_uptake()
! ---------------------------------------------------------------------------------------!
!! Write PO4 uptake to binary for fixed export subroutine
! ---------------------------------------------------------------------------------------!

integer::ios

if(tm_save_PO4_uptake)then

	open(unit=10,form='unformatted',status='replace',action='write',iostat=ios, &
	file='output/'//trim(gen_config_filename)//'/PO4_uptake.bin')
	write(10,iostat=ios) export_save_int
	close(10,iostat=ios)
	print*,'Fixed PO4 uptake saved to: ','output/'//trim(gen_config_filename)//'/PO4_uptake.bin'
endif

end subroutine write_PO4_uptake

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine load_PO4_uptake()
! ---------------------------------------------------------------------------------------!
!! Load PO4 uptake binary for fixed uptake subroutine
! ---------------------------------------------------------------------------------------!

integer::ios
print*,'loading in fixed PO4 uptake data from: ','data/'//trim(tm_PO4uptake_filename)
open(unit=10,status='old',form='unformatted',action='read',IOSTAT=ios, &
file='data/'//trim(tm_PO4uptake_filename))
read(10,iostat=ios) bg_PO4_uptake
close(10,iostat=ios)
end subroutine load_PO4_uptake

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

function v2f(vector)
! ---------------------------------------------------------------------------------------!
!! reorder vector into a 3D field
! ---------------------------------------------------------------------------------------!

real,dimension(NX,NY,NZ)::v2f
real,dimension(tm_nbox)::vector
integer::n

v2f(:,:,:)=0.0
do n = 1,tm_nbox
	 v2f(tm_i(n),tm_j(n),tm_k(n)) = vector(n);
enddo

end function v2f


! ---------------------------------------------------------------------------------------!
! END OF MODULE
! ---------------------------------------------------------------------------------------!

end module io_module
