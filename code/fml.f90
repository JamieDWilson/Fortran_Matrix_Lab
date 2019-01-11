PROGRAM fml

use fml_lib
use tm_module
use bg_module

implicit none

! local variables
integer::n,nn,t,p,count,save_count
real::start,finish,sum_val,start2,finish2
count=0
save_count=1
! initialise parameters and arrays
call initialise_model()
! allocate arrays, assign parameter values
call setup_model()
! load transport matrix data
call load_TM_data()

print*,'*************************'
print*,'Running model...'
! run simulation
call cpu_time(start)

do t=1,gen_runtime_years*tm_n_dt

	J(:,:)=0.0
	particles(:,:)=0.0
	export(:)=0.0

	if(mod(t,bg_dt_ratio)==0)THEN ! circulation + biogeochemistry step
		
		call PO4_uptake()
		
		J(:,iPO4)=J(:,iPO4)+amul_remin(Aremin,(particles(:,iPO4))) ! POP remineralisation
		!print*,sum(amul_remin(Aremin,(particles(:,iPO4)))*tm_vol)
		!print*,sum(particles(:,iPO4)*tm_vol)
				
		call DOP_remin()
		
		call cpu_time(start2)
		tracers=amul(Aexp,tracers_1)
		call update_bgc()
		tracers=amul(Aimp,tracers)
		tracers_1=tracers
		call cpu_time(finish2)
		
	else ! circulation only step
		call cpu_time(start2)
		tracers=amul(Aexp,tracers_1)
		tracers=amul(Aimp,tracers)
		tracers_1=tracers

		call cpu_time(finish2)
		
 	end if
	
	! integrate output in last year of run
	if(t.gt.((gen_runtime_years*tm_n_dt)-tm_n_dt))then
		
		call integrate_output(save_count)

	end if
	
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
call write_output_netcdf()
print*,'*************************'




END PROGRAM fml