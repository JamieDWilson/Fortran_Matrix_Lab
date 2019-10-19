module fml_lib

implicit none
save

! ******************* namelist definitions ***********************!
! transport matrix parameters
CHARACTER(LEN=100)::tm_Aexp_filename,tm_Aimp_filename
namelist /fml_namelist/ tm_Aexp_filename,tm_Aimp_filename
real::tm_native_dt=1200.0
real::tm_dt_scale=1
namelist /fml_namelist/tm_native_dt,tm_dt_scale
integer::tm_n_dt
integer::tm_nbox
namelist /fml_namelist/ tm_n_dt,tm_nbox
character(len=100)::tm_PO4restore_filename,tm_PO4uptake_filename
namelist /fml_namelist/ tm_PO4restore_filename,tm_PO4uptake_filename
character(len=100)::tm_bgc_data_filename='TMM_MIT_data.nc'
character(len=100)::tm_grid_data_filename='TMM_MIT_grid.nc'
namelist /fml_namelist/ tm_bgc_data_filename,tm_grid_data_filename
logical::tm_seasonal=.true.
namelist /fml_namelist/ tm_seasonal
character(len=100)::tm_data_fileloc
namelist /fml_namelist/ tm_data_fileloc
logical::tm_save_PO4_uptake=.false.
namelist /fml_namelist/ tm_save_PO4_uptake

! biogeochemical parameters
real::bg_uptake_tau=30.0
namelist /fml_namelist/ bg_uptake_tau
real::bg_DOC_k=0.5
namelist /fml_namelist/ bg_DOC_k
real::bg_dt_ratio= 1.0
namelist /fml_namelist/ bg_dt_ratio
real::bg_DOC_frac =0.66
namelist / fml_namelist / bg_DOC_frac
integer::bg_n_euphotic_lyrs = 2
namelist / fml_namelist / bg_n_euphotic_lyrs
character(len=100)::bg_uptake_function = 'restore'
namelist / fml_namelist / bg_uptake_function
logical::bg_O_select = .false.
namelist / fml_namelist / bg_O_select
logical::bg_C_select = .false.
namelist / fml_namelist / bg_C_select
logical::bg_restore_atm_CO2=.false.
namelist / fml_namelist / bg_restore_atm_CO2
real::bg_restore_atm_CO2_target=278.0
namelist / fml_namelist / bg_restore_atm_CO2_target
real::bg_gastransfer_a=6.97e-7
namelist / fml_namelist / bg_gastransfer_a
real::bg_rain_ratio,bg_CaCO3_length_scale=2100.0
namelist / fml_namelist / bg_rain_ratio,bg_CaCO3_length_scale
real::bg_martin_remin_b=-0.858
namelist / fml_namelist / bg_martin_remin_b
logical::bg_martin_remin_spatial=.false.
namelist / fml_namelist / bg_martin_remin_spatial
character(len=100)::bg_martin_b_input_filename=''
namelist /fml_namelist / bg_martin_b_input_filename

! general model parameters
integer::gen_n_tracers=2
namelist / fml_namelist / gen_n_tracers
integer::gen_runtime_years=1
namelist / fml_namelist / gen_runtime_years
character(len=100)::gen_save_timeseries_file
character(len=100)::gen_save_timeslice_file
namelist /fml_namelist/ gen_save_timeseries_file,gen_save_timeslice_file


! end of namelist definitions

! define sparse matrix type
type sparse
	real,allocatable,dimension(:,:)::val_n
	real,allocatable,dimension(:)::val
	integer(KIND=4),allocatable,dimension(:)::row
	integer(KIND=4),allocatable,dimension(:)::col
	integer::nnz
	integer::nb
	integer::n_time
end type sparse

! define the sparse matrices
type(sparse)::Aexp
type(sparse)::Aimp
!type(sparse)::Aremin
!type(sparse)::I ! identity matrix
type(sparse)::Aconv ! convert indices


! ******************* allocatable ***********************!
real,dimension(:,:),allocatable::tm_seaice_frac
real,dimension(:,:),allocatable::tm_windspeed
real,dimension(:,:),allocatable::tm_T
real,dimension(:,:),allocatable::tm_S
real,dimension(:,:),allocatable::tm_silica
real,dimension(:),allocatable::tm_area,tm_vol,tm_lon,tm_lat,tm_depth,tm_depth_btm
integer,dimension(:),allocatable::tm_i,tm_j,tm_k,tm_wc


real,dimension(:),allocatable::seaice_dt
real,dimension(:),allocatable::wind_dt
real,dimension(:),allocatable::T_dt
real,dimension(:),allocatable::S_dt
real,dimension(:),allocatable::silica_dt

real,dimension(:),allocatable::bg_martin_b
real,dimension(:,:),allocatable::bg_PO4_obs
real,dimension(:,:),allocatable::bg_PO4_uptake

real,dimension(:,:),allocatable::tracers
real,dimension(:,:),allocatable::tracers_1
real,dimension(:,:),allocatable::C
real,dimension(:,:),allocatable::C_consts
real,dimension(:,:),allocatable::J
real,dimension(:,:),allocatable::Jatm
real,dimension(:,:),allocatable::particles
real,dimension(:),allocatable::ATM
real,dimension(:),allocatable::export_save
real,dimension(:,:),allocatable::export_save_int
real,dimension(:,:),allocatable::diag
real,dimension(:,:),allocatable::tracers_int
real,dimension(:),allocatable::ATM_int
real::t_int=0.0
real,dimension(:,:),allocatable::diag_int

integer,dimension(:),allocatable::iSur
real,dimension(:),allocatable::tm_timeseries,tm_timeslice



! ******************* global variables ***********************!
real::gen_conv_d_yr=365.25
real::bg_DOC_rfrac
REAL::bg_PO4_init=0.00217 ! initial PO4 (mmol m-3 -> mol m-3) (Kriest et al., 2010)
REAL::bg_DOC_init=0.0001*1E-03 ! inital DOP (mmol m-3 -> mol m-3) (Kriest et al., 2010)
real::bg_DIC_init=2.299 !(mmol m-3 -> mol m-3)
real::bg_ALK_init=2.420 !(mmol m-3 -> mol m-3)
real::bg_dt
integer::ioPO4,ioDOP,ioDIC,ioALK
integer::ioCO2,ioCO3,ioH,iopCO2
integer::iK1,iK2,iKw,iKp1,iKp2,iKp3,iKSi,iKb,iK0
integer::iaO2,iaCO2
integer::isPOP,isCaCO3

real,dimension(:),allocatable::tm_seasonal_scale
real,dimension(:),allocatable::tm_seasonal_rscale
integer,dimension(:),allocatable::tm_seasonal_n1
integer,dimension(:),allocatable::tm_seasonal_n2
integer::dt_count
integer::n_seasonal
integer::n_euphotic_boxes
integer::n_surface_boxes

logical::gen_restart_select=.false. ! default value that will reset if restart selected at cmd line
character(len=100)::gen_config_filename,gen_restart_filename

real::bg_C_to_P=116.0
real::bg_N_to_P=16.0

real,dimension(5,2)::Sc_coeffs ! 4 if using older values
real,dimension(6,2)::Bunsen_coeffs
real,dimension(7,2)::Sol_Orr

real::conv_sec_yr=60.0*60.0*24.0*365.25

integer::n_ATM_tracers
real::ATM_vol
real::ATM_mol
real::rho=1024.5 ! kg m-3
real::r_rho=1.0/1024.5 ! m3 kg-1

real::carbchem_tol=0.01 ! tolerance for H+ convergence

real::tm_save_intra_freq=1.0
integer::timeseries_count=1
integer::timeslice_count=1

! grid-dimensions: hard-coded!
integer::NZ = 15;
integer::NY = 64;
integer::NX = 128;

real::tm_dt

contains
! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!
subroutine load_namelist()

! local variables
integer::ios,nn,n,i
character(len=100)::arg

! read command line variables
do i=1,iargc()
	CALL getarg(i,arg)

	if(i.eq.1)then
		gen_config_filename=arg ! config
	else if(i.eq.2)then
		read(arg,*) gen_runtime_years ! runtime
	else if(i.eq.3)then
		gen_restart_filename=arg ! restart
		gen_restart_select=.true.
	end if
end do

! read in namelist file
print*,trim(gen_config_filename)
open(unit=20,file='../experiments/'//trim(gen_config_filename),status='old',action='read')
read(unit=20,nml=fml_namelist,iostat=ios)
close(unit=20)


end subroutine load_namelist

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

end module fml_lib
