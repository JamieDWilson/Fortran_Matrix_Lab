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

subroutine calc_C_consts()

!integer::BT ! total boron
real::T,S
integer::n

! to do: add C_consts to fml_lib.f90
do n=1,n_surface_boxes

	T = tm_T(n,1)+273.15
	S = tm_S(n,1)
	
	! K1
	C_consts(n,iK1)=exp(2.83655-2307.1266/T-1.5529413*log(T) &
	-(0.207608410+4.0484/T)*sqrt(S) &
	+0.0846834*S-0.00654208*S**(3.0/2.0)+log(1.0-0.001005*S))
	
	! K2
	C_consts(n,iK2)=exp( &
	-9.226508-3351.6106/T-0.2005743*log(T) &
	-(0.106901773+23.9722/T)*sqrt(S) &
	+0.1130822*S-0.00846934*S**(3.0/2.0)+log(1.0-0.001005*S))
	
	! K0
	C_consts(n,iK0)=exp( &
	9345.17/T-60.2409+23.3585*log(T/100.0) &
	+S*(0.023517-0.00023656*T+0.0047036*(T/100.0)**2))
	
	! KB
	C_consts(n,iKb)=exp( &
	(-8966.90-2890.53*S**0.5-77.942*S+1.728*S**(3.0/2.0)-0.0996*S**2.0)/T &
	+148.0248+137.1942*S**0.5+1.62142*S &
	-(24.4344+25.085*S**0.5+0.2474*S)*log(T) &
	+0.053105*S**0.5*T)
	
	! Kw
	C_consts(n,iKw)=exp( &
	148.96502-13847.26/T-23.6521*log(T) &
	+(118.67/T-5.977+1.0495*log(T))*S**0.5-0.01615*S)
	
	! KSi - needs checking!
	!C_consts(n,iKSi)=exp( &
	!-8904.2/T+117.385-19.334*log(T) &
	!+(3.5913-458.79/T)*T**0.5 &
	!+(188.74/T-1.5998)*T &
	!+(0.07871-12.1652/T)*T**2 &
	!+log(1.0-0.001005*S))
	
	C_consts(n,iKsi)=1.0
	
	! KP1
	C_consts(n,iKp1)=exp( &
	-4576.752/T+115.525-18.453*log(T) &
	+(-106.736/T+0.69171)*S**0.5 &
	+(-0.65643/T-0.01844)*S)
	
	! KP2
	C_consts(n,iKp2)=exp( &
	-8814.715/T+172.0883-27.927*log(T) &
	+(-160.34/T+1.3566)*S**0.5 &
	+(0.37335/T-0.05778)*S)
	
	! KP3
	C_consts(n,iKp3)=exp( &
	-3070.75/T-18.141 &
	+(17.27039/T+2.81197)*S**0.5 &
	+(-44.99486/T-0.09984)*S)
	
	

enddo 
!print*,'T',tm_T(1,1)+273.15
!print*,'S',tm_S(1,1)
!print*,'K1',C_consts(1,iK1)
!print*,'K2',C_consts(1,iK2)
!print*,'K0',C_consts(1,iK0)
!print*,'Kw',C_consts(1,iKw)
!print*,'KB',C_consts(1,iKB)
!print*,'KSi',C_consts(1,iKSi)
!print*,'Kp1',C_consts(1,iKp1)
!print*,'Kp2',C_consts(1,iKp2)
!print*,'Kp3',C_consts(1,iKp3)

end subroutine calc_C_consts

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine calc_pCO2()
!............................................................
! Solve carbonate system for pC02
! M. Follows, T. Ito, S. Dutkiewicz (2006) in Ocean Modelling
!.............................................................

! local variables
real::pt,sit,ta,pCO2,dic,H,bt,k1,k2,k1p,k2p,k3p,kb,kw,ksi,k0 ! ff in original but not used, added k0
real::gamm,co2s,hg,cag,bohg,h3po4g,h2po4g,hpo4g,po4g,siooh3g,denom,dummy,fg
integer::n

! dic = dissolved inorganic carbon; pt = dissolved inorganic phosphorus
! sit = dissolved inorganic silica, bt = dissolved inorganic boron
! ta = total alkalinity; ca = carbonate alkalinity; H = [H+]
! pCO2 = partial pressure CO2; ff = fugacity of CO2
! k1, k2 = carbonate equilibrium coeffs; kw = dissociation of water
! klp, k2p, k3p = phosphate equilibrium coefficients
! ksi, kb = silicate and borate equilibrium coefficients
! Equilibrium relationships from DOE handbook (DOE, 1994):
! coefficients evaluated elsewhere and passed in.

do n=1,n_surface_boxes

! initialise variables
dic=tracers_1(n,iDIC)
ta=tracers_1(n,iALK)
pt=tracers_1(n,iPO4)
sit=0.0
bt=4.16e-4*(tm_S(n,1)/35.0) ! total boron conc. from ZW2001
k1=C_consts(n,iK1)
k2=C_consts(n,iK2)
kw=C_consts(n,iKw)
k1p=C_consts(n,iKp1)
k2p=C_consts(n,iKp2)
k3p=C_consts(n,iKp3)
ksi=C_consts(n,iKSi)
kb=C_consts(n,iKb)
k0=C_consts(n,iK0)

! First guess of [H+]: from last timestep *OR* fixed for cold start

if(C(n,iH).eq.-1.0)then
	hg=1.0e-9 ! cold start
	else
	hg=C(n,iH) ! previous timestep
endif

! estimate contributions to total alk from borate, silicate, phosphate
bohg = bt*kb/(hg + kb)
siooh3g = sit*ksi/(ksi + hg)
denom = hg*hg*hg + (k1p*hg*hg) + (k1p*k2p*hg) + (k1p*k2p*k3p)
h3po4g = (pt*hg*hg*hg)/denom
h2po4g = (pt*k1p*hg*hg)/denom
hpo4g = (pt*k1p*k2p*hg)/denom
po4g = (pt*k1p*k2p*k3p)/denom
! estimate carbonate alkalinity
fg = - bohg - (kw/hg) + hg - hpo4g - 2.0*po4g + h3po4g - siooh3g
cag = ta + fg
! improved estimate of hydrogen ion conc
gamm = dic/cag
dummy = (1.0-gamm)*(1.0-gamm)*k1*k1 - 4.0*k1*k2*(1.0 - 2.0*gamm)
H = 0.5*((gamm-1.0)*k1 + sqrt(dummy))
! evaluate [CO2*]
co2s = dic/(1.0 + (k1/H) + (k1*k2/(H*H)))
! evaluate surface pCO2
C(n,ipCO2) = co2s/k0
C(n,iH) = H

enddo

end subroutine calc_pCO2

end module bg_module