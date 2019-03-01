module io_module

use fml_lib
use netcdf

implicit none
contains

! ---------------------------------------------------------------------------------------!
! load_TM_data
! - calls routines to load TM data
! ---------------------------------------------------------------------------------------!

subroutine load_TM_data()

call load_TM_netcdf('../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aexp_filename),Aexp)
call load_TM_netcdf('../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aimp_filename),Aimp)
call load_TM_netcdf(tm_Aremin_filename,Aremin)

call load_TM_grid_data()
call load_TM_bgc_data()

select case(trim(bg_uptake_function))
case('restore')
	call load_PO4_restore()
case('fixed')
	call load_PO4_uptake()
end select

end subroutine load_TM_data

! ---------------------------------------------------------------------------------------!
! load_TM_netcdf
! - loads transport matrix data
! ---------------------------------------------------------------------------------------!

SUBROUTINE load_TM_netcdf(dum_filename,dum_A)

character(len=*)::dum_filename
type(sparse)::dum_A

integer::loc_ncid,loc_varid,status
character(len=100)::loc_lname

print*,'loading TM data from:',trim(dum_filename)

! open netcdf file
status=nf90_open(trim(dum_filename),nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),dum_filename

! matrix values
status=nf90_inq_varid(loc_ncid,'A_val',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'A_val'

status=nf90_get_var(loc_ncid,loc_varid,dum_A%val)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'A_val'

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
! load_TM_grid_data
! - loads transport matrix grid variables
! ---------------------------------------------------------------------------------------!

subroutine load_TM_grid_data()

! local variables
integer::n,status,loc_varid,loc_ncid

! open netcdf file
status=nf90_open('../data/'//trim(tm_data_fileloc)//'/'//trim(tm_grid_data_filename) , nf90_nowrite,loc_ncid)
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
! load_TM_bgc_data
! - loads biogeochemical-relevant data for transport matrix
! ---------------------------------------------------------------------------------------!


subroutine load_TM_bgc_data()

! local variables
integer::n,status,loc_varid,loc_ncid

! open netcdf file
status=nf90_open('../data/'//trim(tm_data_fileloc)//'/'//trim(tm_bgc_data_filename) , nf90_nowrite,loc_ncid)
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
! load_PO4_restore
! - loads PO4 observations
! ---------------------------------------------------------------------------------------!

subroutine load_PO4_restore()

! local variables
integer::n,status,loc_varid,loc_ncid

!print*,'Reading in PO4restore data from:','../data/'//trim(tm_data_fileloc)//'/'//trim(tm_PO4restore_filename)

! open netcdf file
status=nf90_open('../data/'//trim(tm_data_fileloc)//'/'//trim(tm_PO4restore_filename), nf90_nowrite,loc_ncid)
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
! load_PO4_uptake
! - loads PO4 uptake from previous run
! ---------------------------------------------------------------------------------------!


subroutine load_PO4_uptake()

! local variables
integer::n,status,loc_varid,loc_ncid

!print*,'Reading in PO4uptake data from:','../data/'//trim(tm_data_fileloc)//'/'//trim(tm_PO4uptake_filename)

! open netcdf file
status=nf90_open('../data/'//trim(tm_data_fileloc)//'/'//trim(tm_PO4uptake_filename), nf90_nowrite,loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),tm_PO4uptake_filename

! matrix values
status=nf90_inq_varid(loc_ncid,'PO4_Uptake',loc_varid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'PO4uptake'

status=nf90_get_var(loc_ncid,loc_varid,bg_PO4_uptake,start=(/ 1, 1 /),count=(/ n_euphotic_boxes , 12 /))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)),'PO4uptake'
!print*,bg_PO4_uptake(1:10,1)

! close netcdf file
status=nf90_close(loc_ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

end subroutine load_PO4_uptake

! ---------------------------------------------------------------------------------------!
! load_data_saving
! - loads times for saving timeslice and timeseries
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

end subroutine load_data_saving

! ---------------------------------------------------------------------------------------!
! initialise_output
! - creates output and restart directories if needed
! - calls subroutines to initialise timeslice and timeseries output
! ---------------------------------------------------------------------------------------!


subroutine initialise_output()

logical::exist_dir

! make output directory (if doesn't already exist)
INQUIRE(FILE='../output/'//trim(gen_config_filename), EXIST=exist_dir)
if(exist_dir.eqv..false.)then
	call system('mkdir ../output/'//trim(gen_config_filename))
end if

! make restart directory (if doesn't already exist)
INQUIRE(FILE='../output/'//trim(gen_config_filename)//'/restart', EXIST=exist_dir)
if(exist_dir.eqv..false.)then
	call system('mkdir ../output/'//trim(gen_config_filename)//'/restart')
end if

! initialise timeseries files
call initialise_timeseries_output()

end subroutine initialise_output

! ---------------------------------------------------------------------------------------!
! intitialise_timeseries_output
! - write headers for timeseries output
! ---------------------------------------------------------------------------------------!

subroutine initialise_timeseries_output()

! atm CO2
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_atm_CO2.dat',status='replace')
write(unit=10,fmt='(A44)') '% / year / atmospheric CO2 (mol) / atmospheric CO2 (ppmv)'
close(unit=10)

! PO4
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_PO4.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / PO4 inventory (mol) / [PO4] global (mmol m-3) / [PO4] surface (mmol m-3)'
close(unit=10)

! DOP
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_DOP.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / DOP inventory (mol) / [DOP] global (mmol m-3)'
close(unit=10)

! DIC
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_DIC.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / DIC inventory (mol) / [DIC] global (mmol m-3) / [DIC] surface (mmol m-3)'
close(unit=10)

! ALK
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_ALK.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / ALK inventory (mol) / [ALK] global (mmol m-3)'
close(unit=10)

! gas exchange
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_airsea_CO2.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / global mean flux (mol yr-1) / global mean flux (mol m-2 yr-1)'
close(unit=10)

! POM EXPORT
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_POM_export.dat',status='replace')
write(unit=10,fmt='(A100)') '% / year / global mean flux (mol yr-1) / global mean flux (mol m-2 yr-1)'
close(unit=10)



end subroutine initialise_timeseries_output

! ---------------------------------------------------------------------------------------!
! write_timeseries_output
! - writes output to timeseries output files
! ---------------------------------------------------------------------------------------!


subroutine write_timeseries_output()

real::rvol_tot,rvol_sur_tot,rvol_deep_tot

rvol_tot=1.0/sum(tm_vol)
rvol_sur_tot=1.0/sum(tm_vol(1:n_surface_boxes))
rvol_deep_tot=1.0/sum(tm_vol(41153:52749))

! CO2
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_atm_CO2.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6)') &
t_int , &
ATM_int(iaCO2)*ATM_vol, &
ATM_int(iaCO2)*1.0e6
close(unit=10)

!PO4
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_PO4.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6,f12.6,f12.6)') &
t_int , &
sum(tracers_int(:,ioPO4)*tm_vol) , &
sum(tracers_int(:,ioPO4)*tm_vol)*rvol_tot*1.0e3 , &
sum(tracers_int(1:n_surface_boxes,ioPO4)*tm_vol(1:n_surface_boxes))*rvol_sur_tot*1.0e3 , &
sum(tracers_int(41153:52749,ioPO4)*tm_vol(41153:52749))*rvol_deep_tot*1.0e3
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
write(unit=10,fmt='(f12.1,e20.12,f12.6,f12.6,f12.6)') &
t_int , &
sum(tracers_int(:,ioDIC)*tm_vol) , &
sum(tracers_int(:,ioDIC)*tm_vol)*rvol_tot*1.0e3 , &
sum(tracers_int(1:n_surface_boxes,ioDIC)*tm_vol(1:n_surface_boxes))*rvol_sur_tot*1.0e3 , &
sum(tracers_int(41153:52749,ioDIC)*tm_vol(41153:52749))*rvol_deep_tot*1.0e3
close(unit=10)

!ALK
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_ocn_ALK.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,f12.6)') &
t_int , &
sum(tracers_int(:,ioALK)*tm_vol) , &
sum(tracers_int(:,ioALK)*tm_vol)*rvol_tot*1.0e3
close(unit=10)

!gas exchange
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_airsea_CO2.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12,e20.12)') &
t_int , &
sum(diag_int(:,2)), &
sum(diag_int(:,2))/sum(tm_area(1:n_surface_boxes)) ! check this
close(unit=10)

!gas exchange
open(unit=10,file='../output/'//trim(gen_config_filename)//'/timeseries_POM_export.dat',position='append')
write(unit=10,fmt='(f12.1,e20.12)') &
t_int , &
sum(diag_int(:,3))
close(unit=10)




end subroutine write_timeseries_output

! ---------------------------------------------------------------------------------------!
! write_output_netcdf
! - writes output to timeslice output file
! - n.b. currently just last timestep!
! ---------------------------------------------------------------------------------------!


subroutine write_output_netcdf()

integer::ncid,status,m_dimid,n_dimid,PO4id,DOPid,EXPORTid,n_season,DICid,ALKid,airseaCO2id
integer::dimids(2)

!create file
status=nf90_create('../output/'//trim(gen_config_filename)//'/fields_netcdf.nc',nf90_clobber,ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

!define dimensions
status=nf90_def_dim(ncid,'tm_nbox',tm_nbox,m_dimid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_def_dim(ncid,'time',1,n_dimid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

dimids=(/m_dimid,n_dimid/)

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

!status=nf90_def_var(ncid,'EXPORT',nf90_float,dimids,EXPORTid)
!if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

!end definition
status=nf90_enddef(ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

!write data
status=nf90_put_var(ncid,PO4id,tracers(:,ioPO4))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_put_var(ncid,DOPid,tracers(:,ioDOP))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_put_var(ncid,DICid,tracers(:,ioDIC))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_put_var(ncid,ALKid,tracers(:,ioALK))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_put_var(ncid,airseaCO2id,diag(:,1))
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

!status=nf90_put_var(ncid,EXPORTid,EXPORT_int(:,:))
!if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

status=nf90_close(ncid)
if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

print*,'Output written to: ','../output/'//trim(gen_config_filename)//'.nc'

end subroutine write_output_netcdf

! ---------------------------------------------------------------------------------------!
! write_restart
! - binary dump of state variables
! ---------------------------------------------------------------------------------------!

subroutine write_restart()

integer::ios

open(unit=10,form='unformatted',status='replace',action='write',iostat=ios, &
file='../output/'//trim(gen_config_filename)//'/restart/restart.bin')
write(10) tracers_1 , ATM
close(10,iostat=ios)

end subroutine write_restart

! ---------------------------------------------------------------------------------------!
! load_restart
! - load in state variables in binary dump from previous run
! ---------------------------------------------------------------------------------------!

subroutine load_restart()

integer::ios

open(unit=10,status='old',form='unformatted',action='read',IOSTAT=ios, &
file='../output/'//trim(gen_restart_filename)//'/restart/restart.bin')
read(10,iostat=ios) tracers_1 , ATM
close(10,iostat=ios)

end subroutine load_restart

! ---------------------------------------------------------------------------------------!
! END OF MODULE
! ---------------------------------------------------------------------------------------!

end module io_module