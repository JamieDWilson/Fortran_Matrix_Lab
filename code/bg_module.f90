module bg_module

use fml_lib

implicit none
contains

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine PO4_uptake()

integer::n
real::uptake
real,dimension(n_euphotic_boxes)::tmp_Fice,tmp_PO4

! linearly interpolate seaice arrays
tmp_Fice=(tm_seasonal_scale(dt_count)*tm_seaice_frac(:,tm_seasonal_n1(dt_count)))&
+&
((tm_seasonal_rscale(dt_count))*tm_seaice_frac(:,tm_seasonal_n2(dt_count)))

if(bg_PO4restore_select)then

	! linearly interpolate PO4 obs arrays
	tmp_PO4=(tm_seasonal_scale(dt_count)*bg_PO4_obs(:,tm_seasonal_n1(dt_count)))&
	+&
	((tm_seasonal_rscale(dt_count))*bg_PO4_obs(:,tm_seasonal_n2(dt_count)))

	do n=1,n_euphotic_boxes
		if(tracers_1(n,iPO4)>tmp_PO4(n))then
			uptake=tmp_Fice(n)*bg_uptake_tau*(tracers_1(n,iPO4)-tmp_PO4(n)) ! PO4 uptake
			export(n)=uptake ! export for saving output
			
			J(n,iPO4)=J(n,iPO4)-uptake ! PO4
			J(n,iDOP)=J(n,iDOP)+bg_DOC_frac*uptake ! DOP
			particles(n,iPO4)=particles(n,iPO4)+bg_DOC_rfrac*uptake ! POP
			
		end if
	end do
	
	
	else
	! linearly interpolate PO4 uptake arrays
	tmp_PO4=(tm_seasonal_scale(dt_count)*bg_PO4_uptake(:,tm_seasonal_n1(dt_count)))&
	+&
	(tm_seasonal_rscale(dt_count)*bg_PO4_uptake(:,tm_seasonal_n2(dt_count)))

	do n=1,n_euphotic_boxes
		if(tracers_1(n,iPO4)-(tmp_PO4(n)*bg_dt)>0.0)then ! if fixed export does not create negative tracer
			uptake=tmp_PO4(n) ! PO4 uptake
			export(n)=uptake! export for saving output
			
			J(n,iPO4)=J(n,iPO4)-uptake ! PO4
			J(n,iDOP)=J(n,iDOP)+bg_DOC_frac*uptake ! DOP
			particles(n,iPO4)=particles(n,iPO4)+bg_DOC_rfrac*uptake ! POP
			
		end if
	end do
	
end if
	


end subroutine PO4_uptake

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine DOP_remin()

integer::n
real::remin

do n=1,tm_nbox

	if(tracers_1(n,iDOP).gt.1e-8)then ! Kriest et al., (2010) . Also smaller concentrations slow down the matrix calculation significantly.
		remin=tracers_1(n,iDOP)*bg_DOC_k
		
		if(remin*bg_dt.lt.tracers_1(n,iDOP))then ! catch remineralisation making DOP go negative
			J(n,iPO4)=J(n,iPO4)+remin ! DOP remin -> PO4
			J(n,iDOP)=J(n,iDOP)-remin ! DOP remin <- DOP
		end if
	
	end if

end do

end subroutine DOP_remin

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine update_bgc()
integer::n 

		do n=1,tm_nbox
			tracers(n,1)=tracers(n,1)+J(n,1)*bg_dt
		end do
		
		do n=1,tm_nbox
			tracers(n,2)=tracers(n,2)+J(n,2)*bg_dt
		end do
			

end subroutine update_bgc

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine integrate_output(save_count)

integer,intent(inOUT)::save_count

tracers_PO4_int(:,save_count)=tracers_PO4_int(:,save_count)+tracers(:,iPO4)*tm_dt*n_seasonal
tracers_DOP_int(:,save_count)=tracers_DOP_int(:,save_count)+tracers(:,iDOP)*tm_dt*n_seasonal
EXPORT_int(:,save_count)=EXPORT_int(:,save_count)+export(:)*tm_dt*n_seasonal

if(mod(dt_count,tm_n_dt/n_seasonal).eq.0 .and. tm_seasonal)then
	save_count=save_count+1
end if

end subroutine integrate_output

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

end module bg_module