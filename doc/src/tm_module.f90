module tm_module
! ---------------------------------------------------------------------------------------!
!! Transport matrix related subroutines
!!
!! - initialise model
!! - main timestepping subroutine
!! - various matrix operations
! ---------------------------------------------------------------------------------------!

use fml_lib
use io_module

implicit none
contains

! ---------------------------------------------------------------------------------------!

subroutine initialise_model()
! ---------------------------------------------------------------------------------------!
!! Initialises FML model
!!
!! 1) sets up array indices
!!
!! 2) allocates arrays
!!
!! 3) loads in matrix data
!!
!! 4) initialises parameters
!!
!! 5) initialises arrays
! ---------------------------------------------------------------------------------------!

! local variables
integer::n
character(8)::date
character(10)::time
character(5)::zone
integer,dimension(8)::value

! print header text to screen
call date_and_time(date,time,zone,value)
print*,
print*,'**************************'
print*,'*** Fortran (Transport) Matrix Lab ***'
print*,'**************************'
print*,''
print '(1x,I4,A1,I2,A1,I2)',value(1),'/',value(2),'/',value(3)
print '(1x,I2,A1,I2,A1,I2)',value(5),':',value(6),':',value(7)
print*,''
print*,'*************************'
print*,
print*,'Input Directory:'
print*,'../data'//'/'//trim(tm_data_fileloc)
PRINT*,
if(gen_restart_select)then
print*,'Restart Directory:'
print*,'../output/'//trim(gen_restart_filename)
PRINT*,
end if
print*,'Output Directory:'
print*,'../output/'//trim(gen_config_filename)
print*,
print*,'*************************'
print*,
print*,'Initialising model...'
print*,

! -- set-up output files-- !
call load_data_saving()
call initialise_output()

! -- set-up Transport Matrices -- !
call load_TM_metadata('../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aexp_filename),Aexp)
call load_TM_metadata('../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aimp_filename),Aimp)

allocate(Aexp%val_n(Aexp%nnz,Aexp%n_time))
allocate(Aexp%val(Aexp%nnz))
allocate(Aexp%row(Aexp%nb+1))
allocate(Aexp%col(Aexp%nnz))

allocate(Aimp%val_n(Aimp%nnz,Aimp%n_time))
allocate(Aimp%val(Aimp%nnz))
allocate(Aimp%row(Aimp%nb+1))
allocate(Aimp%col(Aimp%nnz))

call load_TM_netcdf('../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aexp_filename),Aexp)
call load_TM_netcdf('../data'//'/'//trim(tm_data_fileloc)//'/'//trim(tm_Aimp_filename),Aimp)

call set_TM_timestep()

tm_nbox=Aexp%nb
n_seasonal=Aexp%n_time

! -- load Transport Matrid grid -- !
allocate(tm_i(tm_nbox))
allocate(tm_j(tm_nbox))
allocate(tm_k(tm_nbox))
allocate(tm_lon(tm_nbox))
allocate(tm_lat(tm_nbox))
allocate(tm_depth(tm_nbox))
allocate(tm_depth_btm(tm_nbox))
allocate(tm_area(tm_nbox))
allocate(tm_vol(tm_nbox))
allocate(tm_wc(tm_nbox))
call load_TM_grid_data()

n_euphotic_boxes=0
n_surface_boxes=0
do n=1,tm_nbox
	if(tm_k(n).le.bg_n_euphotic_lyrs) n_euphotic_boxes=n_euphotic_boxes+1
	if(tm_k(n).eq.1) n_surface_boxes=n_surface_boxes+1
end do
call find_water_columns()

allocate(Aconv%val(tm_nbox))
allocate(Aconv%row(tm_nbox+1))
allocate(Aconv%col(tm_nbox))
call create_Aconv()

! -- load boundary condition data -- !
allocate(tm_T(n_euphotic_boxes,n_seasonal))
allocate(tm_S(n_euphotic_boxes,n_seasonal))
allocate(tm_silica(n_euphotic_boxes,n_seasonal))
allocate(bg_martin_b(tm_nbox))
allocate(tm_seaice_frac(n_euphotic_boxes,n_seasonal))
allocate(tm_windspeed(n_euphotic_boxes,n_seasonal))
call load_TM_bgc_data()

allocate(seaice_dt(n_euphotic_boxes))
allocate(wind_dt(n_euphotic_boxes))
allocate(T_dt(n_euphotic_boxes))
allocate(S_dt(n_euphotic_boxes))
allocate(silica_dt(n_euphotic_boxes))

! -- set-up model arrays -- !
call calc_seasonal_scaling()
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

isPOP=1
isCaCO3=2

n_ATM_tracers=2 ! n.b. 1 does not allow array operations
!if(bg_O_select) n_ATM_tracers=n_ATM_tracers+1
!IF(bg_C_select) n_ATM_tracers=n_ATM_tracers+1


allocate(tracers(tm_nbox,gen_n_tracers))
allocate(tracers_1(tm_nbox,gen_n_tracers))
allocate(C(tm_nbox,4))
allocate(C_consts(tm_nbox,9))
allocate(ATM(n_ATM_tracers))
allocate(J(tm_nbox,gen_n_tracers))
allocate(Jatm(n_surface_boxes,n_ATM_tracers))
allocate(particles(tm_nbox,gen_n_tracers))
allocate(tracers_int(tm_nbox,gen_n_tracers))
allocate(ATM_int(n_ATM_tracers))
allocate(diag(tm_nbox,6)) ! n.b. second dimension hard-coded currently
allocate(diag_int(tm_nbox,6)) ! n.b. second dimension hard-coded currently

! if saving output for another run
if(tm_save_PO4_uptake)then
	allocate(export_save(n_euphotic_boxes))
	allocate(export_save_int(n_euphotic_boxes,tm_n_dt))
endif

select case(trim(bg_uptake_function))
case('restore')
	allocate(bg_PO4_obs(n_euphotic_boxes,n_seasonal))
case('fixed')
	allocate(bg_PO4_uptake(n_euphotic_boxes,tm_n_dt))
end select

! -- initialise model arrays -- !
if(gen_restart_select)then
	call load_restart()
else
	tracers_1(:,ioPO4)=bg_PO4_init
	tracers_1(:,ioDOP)=bg_DOC_init
	tracers_1(:,ioDIC)=bg_DIC_init
	tracers_1(:,ioALK)=bg_ALK_init
	ATM(iaCO2)=278.0*1.0e-6
end if
tracers(:,:)=0.0
J(:,:)=0.0
tracers_int(:,:)=0.0
ATM_int(:)=0.0
diag_int(:,:)=0.0
Jatm(:,:)=0.0
dt_count=1 ! keep track of how many timesteps have passed in one year

! convert parameters to correct units

bg_DOC_rfrac=1.0-bg_DOC_frac ! reciprical of DOC fraction
!bg_DOC_k=1.0/bg_DOC_k ! year-1
bg_uptake_tau=(1.0/bg_uptake_tau)*gen_conv_d_yr ! days to years-1

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
!ATM(iaCO2)=278.0*1.0e-6

! print final header
print*,
print*,'*************************'
print*,
print*,'Running model...'
print*,
print*,'*************************'
print*,
print'(A7,A4,A7,A3,A8,A3,A5,A3,A8,A3,A9,A3)', &
'       ', &
'year', &
'       ', &
'PO4', &
'        ', &
'DOP', &
'     ', &
'DIC', &
'        ', &
'ALK', &
'         ', &
'CO2'

call print_to_screen(0,0.0)


end subroutine initialise_model

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine set_TM_timestep()
! ---------------------------------------------------------------------------------------!
!! Applies timestep to transport matrices
! ---------------------------------------------------------------------------------------!

! local variables
integer::n,nn
real::exponent,t

! Aexp = I+m(Aexp)
Aexp%val_n=Aexp%val_n*tm_dt_scale*tm_native_dt

do n=1,Aexp%nb
	do nn=Aexp%row(n),Aexp%row(n+1)-1
		if(Aexp%col(nn).eq.n) Aexp%val_n(nn,:)=1.0+Aexp%val_n(nn,:)
	end do
end do

! allocate temporary Aimp copy
! Aimp accumulates results
allocate(Apow1%val_n(Aimp%nnz,Aimp%n_time))
allocate(Apow1%val(Aimp%nnz))
allocate(Apow1%row(Aimp%nb+1))
allocate(Apow1%col(Aimp%nnz))

! copy Aimp
Apow1%val_n=Aimp%val_n
Apow1%val=Aimp%val
Apow1%row=Aimp%row
Apow1%col=Aimp%col
Apow1%nnz=Aimp%nnz
Apow1%nb=Aimp%nb
Apow1%n_time=Aimp%n_time

! calculate Aimp**tm_dt_scale
exponent=tm_dt_scale
do
	t=mod(exponent,2.0)
	exponent=floor(exponent/2.0)

	if(t.eq.1.0)then
		call amub(Aimp,Apow1)
	endif

	if(exponent.eq.0.0)then
		exit
	endif

	call amub(Apow1,Apow1)
enddo

! deallocate temporary Aimp copy
deallocate(Apow1%val_n)
deallocate(Apow1%val)
deallocate(Apow1%row)
deallocate(Apow1%col)

! timestep
tm_dt=tm_native_dt*tm_dt_scale*(1.0/(60.0*60.0*24.0*360.0)) ! yr
bg_dt=tm_dt*bg_dt_ratio ! BGC timestep
tm_n_dt=nint(1.0/tm_dt)

end subroutine set_TM_timestep

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine calc_seasonal_scaling
! ---------------------------------------------------------------------------------------!
!! Calculates indices for interpolating the seasonal transport matrices
! ---------------------------------------------------------------------------------------!

! local variables
integer::n,count,nn,n_dt_season

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

n_dt_season=int(real(tm_n_dt)/real(12))

count=1
do n=1,tm_n_dt,n_dt_season
	do nn=0,(n_dt_season)-1
		tm_seasonal_n1(n+nn)=count
		tm_seasonal_n2(n+nn)=count+1
	end do
	count=count+1
end do
where(tm_seasonal_n2==13) tm_seasonal_n2=1

do n=1,tm_n_dt,n_dt_season
	DO nn=0,(n_dt_season)-1
		tm_seasonal_rscale(nn+n)=real(nn)/(real(n_dt_season)-1)
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

subroutine timestep_fml()
! ---------------------------------------------------------------------------------------!
!! Integrate model forward in time
!!
!! - integrates model forward in time using fixed timesteps
!!
!! - calculates  \(\mathbf{A_{imp}} * ( \mathbf{A_{exp}}*\mathbf{c}+\mathbf{q}) \) where:
!!
!! \(\mathbf{A_{imp}}\) - implicit matrix
!!
!! \(\mathbf{A_{exp}}\) - explicit matrix
!!
!! \(\mathbf{c}\) - state variable vector
!!
!! \(\mathbf{q}\) - source/sink vector
! ---------------------------------------------------------------------------------------!

integer::n

! explicit step (Aexp*C)
do n=1,gen_n_tracers
	tracers(:,n)=amul(Aexp,tracers_1(:,n))
enddo

! source/sinks (+q)
tracers(:,:)=tracers(:,:)+J(:,:)*bg_dt
ATM(iaCO2)=ATM(iaCO2)+sum(Jatm(:,iaCO2))*bg_dt/ATM_mol

! implicit step (Aimp*...)
do n=1,gen_n_tracers
	tracers(:,n)=amul(Aimp,tracers(:,n))
enddo

tracers_1=tracers

end subroutine timestep_fml

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

! FUNCTION amul(A,Vector)
! ! output
! REAL,dimension(tm_nbox,gen_n_tracers)::amul
! ! dummy
! type(sparse),intent(in)::A
! REAL,INTENT(in),dimension(tm_nbox,gen_n_tracers)::Vector
! ! local
! integer::n,nn,i
! real::sum_val
! integer,dimension(2)::vector_size
! real,dimension(A%nnz)::val_tmp
!
! vector_size=shape(Vector)
!
! val_tmp=(tm_seasonal_scale(dt_count)*A%val(:,tm_seasonal_n1(dt_count)))&
! +&
! ((tm_seasonal_rscale(dt_count))*A%val(:,tm_seasonal_n2(dt_count)))
!
! DO i=1,vector_size(2)
! do n=1,tm_nbox
! 	sum_val=0.0
!
! 	do nn=A%row(n),A%row(n+1)-1
! 		sum_val=sum_val+val_tmp(nn)*Vector(A%col(nn),i)
! 	end do
! 	amul(n,i)=sum_val
! end do
! end do
!
! end FUNCTION

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

! FUNCTION amul_remin(A,Vector)
! ! output
! REAL,dimension(tm_nbox)::amul_remin
! ! dummy
! type(sparse),intent(in)::A
! REAL,INTENT(in),dimension(tm_nbox)::Vector
! ! local
! integer::n,nn,i
! real::sum_val
!
!
! do n=1,tm_nbox
! 	sum_val=0.0
!
! 	do nn=A%row(n),A%row(n+1)-1
! 		sum_val=sum_val+A%val(nn,1)*Vector(A%col(nn))
! 		!print*,nn,A%row(n),A%row(n+1)-1,A%val(nn,1),Vector(A%col(nn))
! 	end do
! 	amul_remin(n)=sum_val
! end do
! !stop
!
! end FUNCTION

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

FUNCTION amul(A,vector)
! ---------------------------------------------------------------------------------------!
!! sparse matrix - vector multiplication (CSR format)
!!
!! adapted from SPARSEKIT (https://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.html)
! ---------------------------------------------------------------------------------------!

! output
REAL,dimension(tm_nbox)::amul
! dummy
type(sparse),intent(in)::A
!! sparse matrix
REAL,INTENT(in),dimension(tm_nbox)::vector
!! vector

! local
integer::n,nn
real::sum_val

do n=1,size(vector)
	sum_val=0.0

	do nn=A%row(n),A%row(n+1)-1
		sum_val=sum_val+A%val(nn)*vector(A%col(nn))
	end do

	amul(n)=sum_val

end do

end FUNCTION

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

FUNCTION amul_transpose(A,vector)
! ---------------------------------------------------------------------------------------!
!! transpose sparse matrix - vector multiplication (CSR format)
!!
!! adapted from SPARSEKIT (https://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.html)
! ---------------------------------------------------------------------------------------!

! output
REAL,dimension(tm_nbox)::amul_transpose
! dummy
type(sparse),intent(in)::A
!! sparse matrix
REAL,INTENT(in),dimension(tm_nbox)::vector
!! vector

! local
integer::n,nn,i
real::sum_val
real,dimension(tm_nbox)::tmp

amul_transpose=0.0
do n=1,size(vector)

	do nn=A%row(n),A%row(n+1)-1
		amul_transpose(A%col(nn))=amul_transpose(A%col(nn))+Vector(n)*A%val(nn)
	end do

end do

end FUNCTION

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

function amub_nnz(A,B)
! ---------------------------------------------------------------------------------------!
!! sparse matrix multiplication (CSR format): find nonzeros
!!
!! adapted from SPARSEKIT (https://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.html)
! ---------------------------------------------------------------------------------------!

! output
integer::amub_nnz
! dummy
type(sparse),intent(in)::A
!! sparse matrix
type(sparse),intent(in)::B
!! sparse matrix

! local
integer::ncol
integer::ncolb
integer::nrow
integer::ii
integer,dimension(size(A%row)-1)::iw
integer::j
integer::jc
integer::jr
integer::k
integer::last
integer::ldg
integer,dimension(size(A%row)-1)::ndegr

! assuming square matrices of same size
ncol=size(A%row)-1
ncolb=size(A%row)-1
nrow=size(A%row)-1

iw(1:ncolb) = 0
ndegr(1:nrow) = 0

do ii = 1, nrow
!
!  For each row of A.
!
	ldg = 0
!
!  End-of-linked list.
!
	last = -1

	do j = A%row(ii), A%row(ii+1)-1
!
!  Row number to be added.
!
			jr = A%col(j)

			do k = B%row(jr), B%row(jr+1)-1
				 jc = B%col(k)
!
!  Add one element to the linked list.
!
				 if ( iw(jc) == 0 ) then
						ldg = ldg + 1
						iw(jc) = last
						last = jc
				 end if

			 end do

	end do

	ndegr(ii) = ldg
!
!  Reset IW to zero.
!
	do k = 1, ldg
		j = iw(last)
		iw(last) = 0
		last = j
	 end do

end do

! output nnz of A*B
amub_nnz = sum ( ndegr(1:nrow) )

end function

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine amub(A,B)
! ---------------------------------------------------------------------------------------!
!! sparse matrix multiplication A*B (CSR format)
!!
!! assumes square Matrix
!!
!! writes result to matrix A
!!
!! adapted from SPARSEKIT (https://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.html)
! ---------------------------------------------------------------------------------------!

! dummy
type(sparse),intent(inout)::A
!! sparse matrix
type(sparse),intent(in)::B
!! sparse matrix

! local
integer::n,nn,i,s
real::sum_val
real,dimension(A%nb)::tmp
integer::len
integer::ierr
integer,dimension(A%nb)::iw
integer::ka,ii,kb,jj,k
integer::jcol,jpos
real::scal
integer::nzmax,ncol,nrow
integer,dimension(size(A%row))::ic
integer,allocatable,dimension(:)::jc
real,allocatable,dimension(:,:)::c

ncol=A%nb ! size of col (assumes square matrix)
nrow=A%nb ! size of row (assumes square matrix)

! find nnz's of new matrix
nzmax=amub_nnz(A,B)

! allocate working arrays
allocate(jc(nzmax))
allocate(c(nzmax,A%n_time))

do s = 1, A%n_time

len = 0
ic(1) = 1
!
!  Initialize IW.
!
iw(1:ncol) = 0

do ii = 1, nrow
!
!  Row I.
!
	do ka = A%row(ii), A%row(ii+1)-1


		scal = A%val_n(ka,s)


		jj = B%col(ka)

		do kb = B%row(jj), B%row(jj+1)-1

				 jcol = B%col(kb)
				 jpos = iw(jcol)

				 if ( jpos == 0 ) then
						len = len + 1
						if ( nzmax < len ) then
							 ierr = ii
							 return
						end if
						jc(len) = jcol
						iw(jcol)= len
						c(len,s) = scal * B%val_n(kb,s)
				 else
							c(jpos,s) = c(jpos,s) + scal * B%val_n(kb,s)
				 end if

		end do

	end do

	do k = ic(ii), len
		iw(jc(k)) = 0
	end do

	ic(ii+1) = len + 1

end do

end do ! seasonal loop

! set A to C=A*B
if(nzmax.ne.B%nnz)then
	! deallocate(A%val_n)
	! deallocate(A%col)
	!
	! allocate(A%val_n(nzmax,A%n_time))
	! allocate(A%col(nzmax))
	print*,'fatal error - A*B results in different number of nonzeros'
	stop
end if

A%val_n=c
A%col=jc
A%row=ic
A%nnz=nzmax

!deallocate(c)
!deallocate(jc)

end subroutine amub

! ---------------------------------------------------------------------------------------!
! aplb
! - matrix + matrix (CSR Format)
! - adapted from SPARSEKIT (https://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.html)
! ---------------------------------------------------------------------------------------!

! FUNCTION aplb(A,B)
! ! output
! type(sparse),intent(inout)::amplb
! ! dummy
! type(sparse),intent(in)::A
! type(sparse),intent(in)::B
! ! local
! integer::n,nn,i
! real::sum_val
! real,dimension(tm_nbox)::tmp
!
! integer::len
! integer,dimension(size(A%row))::ic
! integer,dimension(tm_nbox)::iw
! integer::ka,ii,kb,jj,k
! real::scal
! integer::nzmax,ncol,nrow
!
! nzmax=size(A%val) ! number of nnz's
! ncol=tm_nbox ! size of col (assumes square matrix)
! nrow=tm_nbox ! size of row (assumes square matrix)
!
! ierr = 0
! len = 0
! ic(1) = 1
! iw(1:ncol) = 0
!
! do ii = 1, nrow
! !
! !  Row I.
! !
! 	 do ka = A%row(ii), A%row(ii+1)-1
!
! 			len = len + 1
! 			jcol = A%col(ka)
!
! 			if ( nzmax < len ) then
! 				ierr = ii
! 				return
! 			end if
!
! 			aplb%col(len) = jcol
! 			aplb%val(len) = A%val(ka)
! 			iw(jcol) = len
! 	 end do
!
! 	 do kb = B%row(ii), B%row(ii+1)-1
!
! 			jcol = B%col(kb)
! 			jpos = iw(jcol)
!
! 			if ( jpos == 0 ) then
!
! 				 len = len + 1
!
! 				 if ( nzmax < len ) then
! 					 ierr = ii
! 					 return
! 				 end if
!
! 				 aplb%col(len) = jcol
! 				 aplb%val(len) = B%val(kb)
! 				 iw(jcol)= len
! 			else
! 				 aplb%val(jpos) = aplb%val(jpos) + B%val(kb)
! 			end if
!
! 	 end do
!
! 	 do k = aplb%row(ii), len
! 		 iw(aplb%col(k)) = 0
! 	 end do
!
! 	 aplb%row(ii+1) = len+1
! end do
!
! return
! end
!
! end function aplb

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine print_to_screen(dum_t,dum_extra)
! ---------------------------------------------------------------------------------------!
!! Print output to screen during runtime
! ---------------------------------------------------------------------------------------!

integer::dum_t
!! timestep
real::dum_extra
!! not used

! local
real::vol_rtot

vol_rtot=1.0/sum(tm_vol)

print'(F11.2,F11.2,F11.2,F11.2,F11.2,F11.2)', &
real(dum_t/tm_n_dt), &
sum(tracers_1(:,ioPO4)*tm_vol)*vol_rtot*1.0e3, &
sum(tracers_1(:,ioDOP)*tm_vol)*vol_rtot*1.0e3, &
sum(tracers_1(:,ioDIC)*tm_vol)*vol_rtot*1.0e3, &
sum(tracers_1(:,ioALK)*tm_vol)*vol_rtot*1.0e3, &
ATM(iaCO2)*1.0e6



end subroutine print_to_screen

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine tm_vars_at_dt()
! ---------------------------------------------------------------------------------------!
!! Linearly interpolate seaice, windstress, T, S at current timestep
! ---------------------------------------------------------------------------------------!

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
! *** windspeed is m/s? so adjust this line of code *** !
!print*,'wind m s-1',wind_dt(2000),seaice_dt(2000)
wind_dt=(wind_dt*wind_dt)*bg_gastransfer_a*seaice_dt*conv_sec_yr

! seasonal matrix values
Aexp%val=(tm_seasonal_scale(dt_count)*Aexp%val_n(:,tm_seasonal_n1(dt_count))) &
+&
((tm_seasonal_rscale(dt_count))*Aexp%val_n(:,tm_seasonal_n2(dt_count)))

Aimp%val=(tm_seasonal_scale(dt_count)*Aimp%val_n(:,tm_seasonal_n1(dt_count))) &
+&
((tm_seasonal_rscale(dt_count))*Aimp%val_n(:,tm_seasonal_n2(dt_count)))




end subroutine tm_vars_at_dt

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!


subroutine integrate_output(loc_t,loc_save_count,loc_dt_count)
! ---------------------------------------------------------------------------------------!
!! Average output and write to files if appropriate
! ---------------------------------------------------------------------------------------!
integer,intent(in)::loc_t
!! timestep
integer,intent(inOUT)::loc_save_count
!! counter to keep track of number of save intervals
integer,intent(in)::loc_dt_count
!! number of timesteps through averaging interval

! local
real::scalar

!if(loc_t>=tm_timeseries(timeseries_count)*96-47 .and. loc_t<=tm_timeseries(timeseries_count)*96+48)then

	! integrate
scalar=bg_dt*tm_save_intra_freq

tracers_int(:,:)=tracers_int(:,:)+tracers(:,:)*scalar

ATM_int(:)=ATM_int(:)+ATM(:)*scalar

t_int=t_int+((loc_t-1.0)/(1.0/bg_dt))*scalar

diag_int=diag_int+diag(:,:)*scalar

if(loc_dt_count.eq.tm_n_dt)then ! at end of model year

	if(loc_t==tm_timeseries(timeseries_count)*96+48)then ! if within a timeseries save year
		call write_timeseries_output()
		timeseries_count=timeseries_count+1
	endif

	if(loc_t==tm_timeslice(timeslice_count)*96+48)then ! if within a timeslice save year
		call write_output_netcdf()
		timeslice_count=timeslice_count+1
	endif
	! write when reached end of time period
	! and reset integrating arrays
	!if(loc_t==tm_timeseries(timeseries_count)*96+48)then
		ATM_int(:)=0.0
		tracers_int(:,:)=0.0
		t_int=0.0
		diag_int=0.0
endif

! within the last model year, save export
if(tm_save_PO4_uptake)then

	! if(loc_t>=(gen_runtime_years*tm_n_dt)-95 .and. loc_t<=(gen_runtime_years*tm_n_dt))then ! if within final year
	! 	export_save_int(:,loc_save_count)=export_save_int(:,loc_save_count)+export_save(:)*bg_dt*(real(tm_n_dt)/real(n_seasonal))
	!
	! 	if(mod(real(loc_dt_count),real(tm_n_dt)/real(n_seasonal)).eq.0.0)then ! step through seasons
	! 		loc_save_count=loc_save_count+1
	! 	endif
	!
	! endif

	if(loc_t>=(gen_runtime_years*tm_n_dt)-95 .and. loc_t<=(gen_runtime_years*tm_n_dt))then ! if within final year
		export_save_int(:,loc_save_count)=export_save(:)
		loc_save_count=loc_save_count+1
	endif

endif


end subroutine integrate_output

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine find_water_columns()
! ---------------------------------------------------------------------------------------!
!! Creates index of water columns in model grid
! ---------------------------------------------------------------------------------------!

integer::n,nn
integer::i,j

do n=1,n_surface_boxes
	do nn=1,tm_nbox

		i=tm_i(n)
		j=tm_j(n)


		if(tm_i(nn).eq.i .and. tm_j(nn).eq.j) tm_wc(nn)=n
	enddo
enddo

end subroutine find_water_columns

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine create_Aconv()
! ---------------------------------------------------------------------------------------!
!! Creates a sparse matrix that re-orders vectors to water-column order for remin subroutines
! ---------------------------------------------------------------------------------------!

integer::n,nn,count
integer,dimension(tm_nbox)::tmp,tmp2


count=1
do n=1,n_surface_boxes
	do nn=1,tm_nbox
		if(tm_wc(nn)==n)then
			Aconv%val(count)=1.0
			Aconv%col(count)=nn
			Aconv%row(count)=count
			count=count+1
		endif
	enddo
enddo
Aconv%row(count)=count

end subroutine create_Aconv


end module tm_module
