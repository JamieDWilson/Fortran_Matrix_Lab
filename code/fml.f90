PROGRAM fml

use fml_lib
use tm_module
use bg_module
use io_module

implicit none

! local variables
integer::n,nn,t,p,count,save_count
real::start,finish,sum_val,start2,finish2
count=0
save_count=1

! load user parameters
call load_namelist()
! allocate arrays, assign parameter values
call initialise_model()


! run simulation
call cpu_time(start)

do t=1,gen_runtime_years*tm_n_dt

	J(:,:)=0.0
	particles(:,:)=0.0
	export(:)=0.0
	Jatm(:,:)=0.0

	call tm_vars_at_dt()

	call cpu_time(start2)

	if(mod(t,int(bg_dt_ratio))==0)THEN ! biogeochemistry source/sink

		if(bg_C_select)then
			call calc_C_consts()
			call calc_pCO2()
		end if

		call PO4_uptake()

		!call POP_remin()
		call POM_remin()

		call DOP_remin()

		call calc_gasexchange()

		call restore_atm_CO2()

	end if

	call integrate_model()

	call integrate_output(t,save_count,dt_count)

	call cpu_time(finish2)

	if(mod(t,tm_n_dt)==0.0)then
		call print_to_screen(t,finish2-start2)
	end if

	dt_count=dt_count+1

	! revert timestep counter to zero
	if(dt_count.gt.tm_n_dt)then
		dt_count=1
	end if

end do
call cpu_time(finish)
print*,'Time taken =',finish-start,'seconds'
print*,'*************************'
print*,

! write output to netcdf file
print*,'*************************'
call write_restart()
print*,'*************************'




END PROGRAM fml
