module fml_lib
! ---------------------------------------------------------------------------------------!
!! Namelist definitions of parameters
!!
!! - definitions of key parameters in model
!! - default parameter values are set here
!! - parameter values can be overriden by the user
!! - loads user-defined parameters from experiment file
! ---------------------------------------------------------------------------------------!

implicit none
save

! ******************* namelist definitions ***********************!
! transport matrix parameters
CHARACTER(LEN=100)::tm_Aexp_filename,tm_Aimp_filename
!! sparse matrix netCDF filename
namelist /fml_namelist/ tm_Aexp_filename,tm_Aimp_filename
real::tm_native_dt=1200.0
!! native timestep of transport matrix
real::tm_dt_scale=1
!! multiplier to lengthen timestep
namelist /fml_namelist/tm_native_dt,tm_dt_scale
integer::tm_n_dt
!! number of timesteps
integer::tm_nbox
!! number of wet grid-boxes in matrix grid
namelist /fml_namelist/ tm_n_dt,tm_nbox
character(len=100)::tm_PO4restore_filename,tm_PO4uptake_filename
!! PO4 uptake forcing filenames
namelist /fml_namelist/ tm_PO4restore_filename,tm_PO4uptake_filename
character(len=100)::tm_bgc_data_filename='TMM_MIT_data.nc'
!! filename of transport matrix metdata
character(len=100)::tm_grid_data_filename='TMM_MIT_grid.nc'
!! filename of transport matrix grid data
namelist /fml_namelist/ tm_bgc_data_filename,tm_grid_data_filename
logical::tm_seasonal=.true.
!! flag for seasonal transport matrices
namelist /fml_namelist/ tm_seasonal
character(len=100)::tm_data_fileloc
!! location of transport matrix files
namelist /fml_namelist/ tm_data_fileloc
logical::tm_save_PO4_uptake=.false.
!! flag to save PO4 uptake for fixed uptake
namelist /fml_namelist/ tm_save_PO4_uptake

! biogeochemical parameters
real::bg_uptake_tau=30.0
!! michaelis-menten uptake rate (mol kg\(^{-1}\) yr\(^{-1}\))
namelist /fml_namelist/ bg_uptake_tau
real::bg_DOC_k=0.5
!! lifetime of dissolved organic matter (years)
namelist /fml_namelist/ bg_DOC_k
real::bg_dt_ratio= 1.0
!! biogeochemistry timestep ratio
namelist /fml_namelist/ bg_dt_ratio
real::bg_DOC_frac =0.66
!! fraction of biological uptake to dissolved organic carbon
namelist / fml_namelist / bg_DOC_frac
integer::bg_n_euphotic_lyrs = 2
!! number of vertical layers considered euphotic zone
namelist / fml_namelist / bg_n_euphotic_lyrs
character(len=100)::bg_uptake_function = 'restore'
!! biological uptake function
namelist / fml_namelist / bg_uptake_function
logical::bg_O_select = .false.
!! flag to select oxygen tracer
namelist / fml_namelist / bg_O_select
logical::bg_C_select = .false.
!! flag to select carbon tracers (DIC + ALK)
namelist / fml_namelist / bg_C_select
logical::bg_restore_atm_CO2=.false.
!! flag to restore atmospheric CO2
namelist / fml_namelist / bg_restore_atm_CO2
real::bg_restore_atm_CO2_target=278.0
!! atmospheric CO2 value to restore to (mixing ratio)
namelist / fml_namelist / bg_restore_atm_CO2_target
real::bg_gastransfer_a=6.97e-7
!! gas transfer velocity (m s\(^{-1}\))
namelist / fml_namelist / bg_gastransfer_a
real::bg_rain_ratio,bg_CaCO3_length_scale=2100.0
!! CaCO3 remineralisation length scale (m)
namelist / fml_namelist / bg_rain_ratio,bg_CaCO3_length_scale
real::bg_martin_remin_b=-0.858
!! Globally uniform Martin curve exponent (unitless)
namelist / fml_namelist / bg_martin_remin_b
logical::bg_martin_remin_spatial=.false.
!! flag to select spatially explicit Martin curves
namelist / fml_namelist / bg_martin_remin_spatial
character(len=100)::bg_martin_b_input_filename=''
!! filename of spatially explicit Martin curve exponents
namelist /fml_namelist / bg_martin_b_input_filename

! general model parameters
integer::gen_n_tracers=2
!! number of tracers
namelist / fml_namelist / gen_n_tracers
integer::gen_runtime_years=1
!! runtime (years)
namelist / fml_namelist / gen_runtime_years
character(len=100)::gen_save_timeseries_file
!! filename containing savepoints for timeseries output
character(len=100)::gen_save_timeslice_file
!! filename containing save points for netCDF output
namelist /fml_namelist/ gen_save_timeseries_file,gen_save_timeslice_file


! end of namelist definitions

! define sparse matrix type
type sparse
	real,allocatable,dimension(:,:)::val_n
	!! all transport matrix nonzeros
	real,allocatable,dimension(:)::val
	!! interpolated nonzeros
	integer(KIND=4),allocatable,dimension(:)::row
	!! row indices
	integer(KIND=4),allocatable,dimension(:)::col
	!! column indices
	integer::nnz
	!! number of nonzeros
	integer::nb
	!! matrix dimensions (assumed square)
	integer::n_time
	!! number of matrices
end type sparse

! define the sparse matrices
type(sparse)::Aexp
!! Explicit transport matrix
type(sparse)::Aimp
!! Implicit transport matrix
!type(sparse)::Aremin
!type(sparse)::I ! identity matrix
type(sparse)::Aconv ! convert indices
!! sparse matrix to convert indices

type(sparse)::Apow1 ! A**1 for set_TM_timestep()
!! temporary sparse matrix for multiplying Aimp

! ******************* allocatable ***********************!
real,dimension(:,:),allocatable::tm_seaice_frac
!! seaice fraction (fraction)
real,dimension(:,:),allocatable::tm_windspeed
!! windspeed (m s\(^{-1}\)
real,dimension(:,:),allocatable::tm_T
!! temperature (deg C)
real,dimension(:,:),allocatable::tm_S
!! salinity (PSU)
real,dimension(:,:),allocatable::tm_silica
!! dissolved silica mol kg\(^{-1}\)
real,dimension(:),allocatable::tm_area,tm_vol,tm_lon,tm_lat,tm_depth,tm_depth_btm
!! grid variables
integer,dimension(:),allocatable::tm_i,tm_j,tm_k,tm_wc
!! vector indies


real,dimension(:),allocatable::seaice_dt
!! interpolated seaice
real,dimension(:),allocatable::wind_dt
!! interpolated windspeed
real,dimension(:),allocatable::T_dt
!! interpolated temperature
real,dimension(:),allocatable::S_dt
!! interpolated salinity
real,dimension(:),allocatable::silica_dt
!! interpolated dissolved silica

real,dimension(:),allocatable::bg_martin_b
!! spatially explicit martin exponent
real,dimension(:,:),allocatable::bg_PO4_obs
!! PO4 observations
real,dimension(:,:),allocatable::bg_PO4_uptake
!! fixed PO4 uptake (mol kg\(^{-1}\) yr\(^{-1}\)

real,dimension(:,:),allocatable::tracers
!! state variable array
real,dimension(:,:),allocatable::tracers_1
!! state variable array
real,dimension(:,:),allocatable::C
!! carbonate chemistry array
real,dimension(:,:),allocatable::C_consts
!! carbonate chemistry constants
real,dimension(:,:),allocatable::J
!! source/sink
real,dimension(:,:),allocatable::Jatm
!! sorce/sink for atmosphere
real,dimension(:,:),allocatable::particles
!! particle array
real,dimension(:),allocatable::ATM
!! atmosphere array
real,dimension(:),allocatable::export_save
!! array to save export
real,dimension(:,:),allocatable::export_save_int
!! array to save integrated export
real,dimension(:,:),allocatable::diag
!! array of diagnosed quantities
real,dimension(:,:),allocatable::tracers_int
!! state variables integrated over time period
real,dimension(:),allocatable::ATM_int
!! atmosphere integrated over time period
real::t_int=0.0
!! ?
real,dimension(:,:),allocatable::diag_int
!! diagnosed quantities integrated over time period

integer,dimension(:),allocatable::iSur
!! indices of surface
real,dimension(:),allocatable::tm_timeseries,tm_timeslice
!! saving timepoints



! ******************* global variables ***********************!
real::gen_conv_d_yr=365.25
!! convert days to years
real::bg_DOC_rfrac
!! recipricol of bg_DOC_frac
REAL::bg_PO4_init=0.00217
!! initial PO4 (mmol m-3 -> mol m-3) (Kriest et al., 2010)
REAL::bg_DOC_init=0.0001*1E-03
!! inital DOP (mmol m-3 -> mol m-3) (Kriest et al., 2010)
real::bg_DIC_init=2.299
!! initial DIC (mmol m-3 -> mol m-3)
real::bg_ALK_init=2.420
!! initial Alkalinity (mmol m-3 -> mol m-3)
real::bg_dt
!! biogeochemistry timestep (year)
integer::ioPO4,ioDOP,ioDIC,ioALK
!!indices
integer::ioCO2,ioCO3,ioH,iopCO2
!! indices
integer::iK1,iK2,iKw,iKp1,iKp2,iKp3,iKSi,iKb,iK0
!! indices
integer::iaO2,iaCO2
!! indices
integer::isPOP,isCaCO3
!! indices

real,dimension(:),allocatable::tm_seasonal_scale
!! seasonal interpolation info
real,dimension(:),allocatable::tm_seasonal_rscale
!! seasonal interpolation info
integer,dimension(:),allocatable::tm_seasonal_n1
!! seasonal interpolation info
integer,dimension(:),allocatable::tm_seasonal_n2
!! seasonal interpolation info
integer::dt_count
!! counter for timesteps taken
integer::n_seasonal
!! number of seasonal matrices
integer::n_euphotic_boxes
!! number of vertical layers considered euphotic grid-boxes
integer::n_surface_boxes
!! number of surface grid-boxes

logical::gen_restart_select=.false.
!! flag to select restart
character(len=100)::gen_config_filename,gen_restart_filename

real::bg_C_to_P=116.0
!! Redfield C:P
real::bg_N_to_P=16.0
!! Redfield N:P

real,dimension(5,2)::Sc_coeffs ! 4 if using older values
!! Schmidt number array for gas exchange
real,dimension(6,2)::Bunsen_coeffs
!! Bunsen coefficients for gas exchange
real,dimension(7,2)::Sol_Orr
!! solubility coefficients for gas exchange

real::conv_sec_yr=60.0*60.0*24.0*365.25
!! convert seconds to years

integer::n_ATM_tracers
real::ATM_vol
real::ATM_mol
real::rho=1024.5 ! kg m-3
!! average density
real::r_rho=1.0/1024.5 ! m3 kg-1
!! recipricol of density

real::carbchem_tol=0.01
!! tolerance for H+ convergence in carbonate chemistry subroutine

real::tm_save_intra_freq=1.0
!! number of intra-year integration periods
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
! ---------------------------------------------------------------------------------------!
!! Load user-defined namelist parameters
! ---------------------------------------------------------------------------------------!

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
