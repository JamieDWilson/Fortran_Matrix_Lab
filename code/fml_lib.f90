module fml_lib

implicit none
save

! ******************* namelist definitions ***********************!
! tm parameters
CHARACTER(LEN=100)::tm_Aexp_filename,tm_Aimp_filename,tm_Aremin_filename
namelist /tm_namelist/ tm_Aexp_filename,tm_Aimp_filename, tm_Aremin_filename
real::tm_dt,tm_native_dt
namelist /tm_namelist/ tm_dt,tm_native_dt
integer::tm_n_dt
integer::tm_Aexp_nnz,tm_Aimp_nnz,tm_Aremin_nnz
integer::tm_nbox
namelist /tm_namelist/ tm_n_dt,tm_Aexp_nnz,tm_Aimp_nnz,tm_Aremin_nnz,tm_nbox
character(len=100)::tm_seaice_filename,tm_PO4restore_filename,tm_vol_filename,tm_PO4uptake_filename
namelist /tm_namelist/ tm_seaice_filename,tm_PO4restore_filename,tm_vol_filename,tm_PO4uptake_filename
character(len=100)::tm_bgc_data_filename
namelist /tm_namelist/ tm_bgc_data_filename
logical::tm_seasonal
namelist /tm_namelist/ tm_seasonal
character(len=100)::tm_data_fileloc
namelist /tm_namelist/ tm_data_fileloc

! biogeochemical parameters
real::bg_uptake_tau
REAL::bg_DOC_k
namelist /tm_namelist/ bg_uptake_tau,bg_DOC_k
integer::bg_dt_ratio
namelist /tm_namelist/ bg_dt_ratio
character(len=100)::bg_martin_b_input_filename
namelist /tm_namelist / bg_martin_b_input_filename
real::bg_DOC_frac
namelist / tm_namelist / bg_DOC_frac
integer::bg_n_euphotic_lyrs
namelist / tm_namelist / bg_n_euphotic_lyrs
logical::bg_PO4restore_select
namelist / tm_namelist / bg_PO4restore_select
logical::bg_O_select,bg_C_select
namelist / tm_namelist / bg_O_select,bg_C_select
logical::bg_restore_atm_CO2
namelist / tm_namelist / bg_restore_atm_CO2
real::bg_restore_atm_CO2_target
namelist / tm_namelist / bg_restore_atm_CO2_target

! general model parameters
integer::gen_n_tracers
namelist / tm_namelist / gen_n_tracers
integer::gen_runtime_years
namelist / tm_namelist / gen_runtime_years

! end of namelist definitions

! define sparse matrix type
type sparse
	real,allocatable,dimension(:,:)::val
	integer(KIND=4),allocatable,dimension(:)::row
	integer(KIND=4),allocatable,dimension(:)::col
	integer::nnz
end type sparse

! define the sparse matrices
type(sparse)::Aexp
type(sparse)::Aimp
type(sparse)::Aremin
!type(sparse)::I


! ******************* allocatable ***********************!
! real,dimension(:),allocatable::tm_Aexp,tm_Aimp,tm_Aremin
! integer,dimension(:),allocatable::tm_Aexp_col,tm_Aimp_col,tm_Aremin_col	
! integer,dimension(:),allocatable::tm_Aexp_row,tm_Aimp_row,tm_Aremin_row
!integer(kind=2),dimension(:),allocatable::tm_i,tm_j,tm_k
real,dimension(:,:),allocatable::tm_seaice_frac
real,dimension(:),allocatable::tm_vol
real,dimension(:,:),allocatable::tm_windspeed
real,dimension(:,:),allocatable::tm_T
real,dimension(:,:),allocatable::tm_S
real,dimension(:,:),allocatable::tm_silica
real,dimension(:),allocatable::tm_area

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
real,dimension(:,:),allocatable::particles
real,dimension(:),allocatable::ATM
real,dimension(:),allocatable::export
real,dimension(:,:),allocatable::tracers_PO4_int,tracers_DOP_int,EXPORT_int
real,dimension(:,:),allocatable::ATM_int

integer,dimension(:),allocatable::iSur


! ******************* global variables ***********************!
real::gen_conv_d_yr=365.25
real::bg_DOC_rfrac
REAL::bg_PO4_init=2.17/1e3 ! initial PO4 (mmol m-3 -> mol m-3) (Kriest et al., 2010)
REAL::bg_DOC_init=0.0001/1e3 ! inital DOP (mmol m-3 -> mol m-3) (Kriest et al., 2010)
real::bg_DIC_init=2299.0/1e3 !(mmol m-3 -> mol m-3) 
real::bg_ALK_init=2420.9/1e3 !(mmol m-3 -> mol m-3) 
real::bg_dt
integer::ioPO4,ioDOP,ioDIC,ioALK
integer::ioCO2,ioCO3,ioH
integer::iK1,iK2,iKw,iKp1,iKp2,iKp3,iKSi,iKb,iK0
integer::iaO2,iaCO2

real,dimension(:),allocatable::tm_seasonal_scale
real,dimension(:),allocatable::tm_seasonal_rscale
integer,dimension(:),allocatable::tm_seasonal_n1
integer,dimension(:),allocatable::tm_seasonal_n2
integer::dt_count
integer::n_seasonal
integer::n_euphotic_boxes
integer::n_surface_boxes

logical::gen_restart_select=.false.
character(len=100)::gen_config_filename,gen_restart_filename

real::bg_C_to_P=116.0
real::bg_N_to_P=16.0

real,dimension(5,2)::Sc_coeffs ! 4 if using older values
real,dimension(6,2)::Bunsen_coeffs
real,dimension(7,2)::Sol_Orr 

integer::n_ATM_tracers
real::ATM_vol
real::ATM_mol
real::bg_gastransfer_a=6.97e-7 ! see Orr et al., (2017), eqn. 13
real::rho=1024.5 ! kg m-3
real::r_rho=1.0/1024.5 ! m3 kg-1



contains
! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!
subroutine initialise_model()

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
read(unit=20,nml=tm_namelist,iostat=ios)
close(unit=20)

end subroutine initialise_model

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

end module fml_lib