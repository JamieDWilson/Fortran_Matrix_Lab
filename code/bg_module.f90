module bg_module
! ---------------------------------------------------------------------------------------!
!! Subroutines related to biogeochemistry
! ---------------------------------------------------------------------------------------!

use fml_lib
use tm_module

implicit none
contains

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine PO4_uptake()
! ---------------------------------------------------------------------------------------!
!! Calculates PO4 uptake
!!
!! - calculates biological uptake of PO4
!! - updates other tracers via Redfield ratios
!! - calculates CaCO3 export via a rain-ratio
! ---------------------------------------------------------------------------------------!

integer::n
real::uptake
real::tmp_PO4
real::tmp_tracer

do n=1,n_euphotic_boxes

	uptake=0.0
	tmp_tracer=tracers_1(n,ioPO4)

	select case(trim(bg_uptake_function))
	case('restore')
		tmp_PO4=(tm_seasonal_scale(dt_count)*bg_PO4_obs(n,tm_seasonal_n1(dt_count)))&
		+ &
		((tm_seasonal_rscale(dt_count))*bg_PO4_obs(n,tm_seasonal_n2(dt_count)))
		if(tmp_tracer>tmp_PO4) uptake=seaice_dt(n)*bg_uptake_tau*(tmp_tracer-tmp_PO4) ! PO4 uptake
	case('fixed')
		uptake=0.0 ! set to zero initially, update if
		if(tmp_tracer-(bg_PO4_uptake(n,dt_count)*bg_dt)>0.0) uptake=bg_PO4_uptake(n,dt_count) ! PO4 uptake
	case('abiotic')
		uptake=0.0
	end select

	J(n,ioPO4)=J(n,ioPO4)-uptake ! PO4
	J(n,ioDOP)=J(n,ioDOP)+bg_DOC_frac*uptake ! DOP
	particles(n,isPOP)=particles(n,isPOP)+bg_DOC_rfrac*uptake ! POP


	if(bg_C_select)then
		! OM
		J(n,ioDIC)=J(n,ioDIC)-(uptake*bg_C_to_P)
		J(n,ioALK)=J(n,ioALK)+(uptake*bg_N_to_P)

		! CaCO3 precipitation
		particles(n,isCaCO3)=particles(n,isCaCO3)+(bg_DOC_rfrac*uptake*bg_C_to_P*bg_rain_ratio)
		J(n,ioDIC)=J(n,ioDIC)-(bg_DOC_rfrac*uptake*bg_C_to_P*bg_rain_ratio)
		J(n,ioALK)=J(n,ioALK)-(bg_DOC_rfrac*uptake*bg_C_to_P*bg_rain_ratio*2.0)
	end if

	diag(n,3)=uptake*tm_vol(n) ! mol
	diag(n,6)=(bg_DOC_rfrac*uptake*bg_C_to_P*bg_rain_ratio)*tm_vol(n)

	if(tm_save_PO4_uptake) export_save(n)=uptake! save export

end do


end subroutine PO4_uptake

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine DOP_remin()
! ---------------------------------------------------------------------------------------!
!! Dissolved organic phosphorus remineralisation
! ---------------------------------------------------------------------------------------!

integer::n
real::remin

do n=1,tm_nbox

	if(tracers_1(n,ioDOP).gt.1e-8)then ! Kriest et al., (2010) . Also smaller concentrations slow down the matrix calculation significantly.

		remin=tracers_1(n,ioDOP)*bg_DOC_k

		if(remin*bg_dt.lt.tracers_1(n,ioDOP))then ! catch remineralisation making DOP go negative
			J(n,ioPO4)=J(n,ioPO4)+remin ! DOP remin -> PO4
			J(n,ioDOP)=J(n,ioDOP)-remin ! DOP remin <- DOP

			if(bg_C_select)then
				J(n,ioDIC)=J(n,ioDIC)+(remin*bg_C_to_P)
				J(n,ioALK)=J(n,ioALK)-(remin*bg_N_to_P)
			end if

		end if

	end if

end do

end subroutine DOP_remin

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

! subroutine POP_remin()
!
! real,dimension(tm_nbox)::remin
!
! remin=amul(Aremin,(particles(:,isPOP))) ! POP remineralisation
!
! J(:,ioPO4)=J(:,ioPO4)+remin
!
! if(bg_C_select)then
! 	J(:,ioDIC)=J(:,ioDIC)+(remin*bg_C_to_P)
! 	J(:,ioALK)=J(:,ioALK)-(remin*bg_N_to_P)
! endif
!
! diag(:,4)=remin
!
! end subroutine POP_remin


! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine update_bgc()
! ---------------------------------------------------------------------------------------!
!! add biogeochemical sources/sinks to state arrays
! ---------------------------------------------------------------------------------------!
integer::n,n_tracer

do n_tracer=1,gen_n_tracers
	do n=1,tm_nbox
		tracers(n,n_tracer)=tracers(n,n_tracer)+J(n,n_tracer)*bg_dt
	enddo
enddo

ATM(iaCO2)=ATM(iaCO2)+sum(Jatm(:,iaCO2))*bg_dt/ATM_mol

end subroutine update_bgc

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine calc_C_consts()
! ---------------------------------------------------------------------------------------!
!! Calculate carbonate system constants
!! - uses DOE (1994)
! ---------------------------------------------------------------------------------------!

real::T,S,I
integer::n

! to do: add C_consts to fml_lib.f90
do n=1,n_surface_boxes

	T = T_dt(n)+273.15
	S = S_dt(n)

	! K1 (Roy et al., 1993)
	C_consts(n,iK1)=exp(2.83655-2307.1266/T-1.5529413*log(T) &
	-(0.207608410+4.0484/T)*sqrt(S) &
	+0.0846834*S-0.00654208*S**(3.0/2.0)+log(1.0-0.001005*S))

	! K2 (Roy et al., 1993)
	C_consts(n,iK2)=exp( &
	-9.226508-3351.6106/T-0.2005743*log(T) &
	-(0.106901773+23.9722/T)*sqrt(S) &
	+0.1130822*S-0.00846934*S**(3.0/2.0)+log(1.0-0.001005*S))

	! K0 (Weiss 1974)
	C_consts(n,iK0)=exp( &
	9345.17/T-60.2409+23.3585*log(T/100.0) &
	+S*(0.023517-0.00023656*T+0.0047036*(T/100.0)**2))

	! KB (Dickson 1990b(
	C_consts(n,iKb)=exp( &
	(-8966.90-2890.53*S**0.5-77.942*S+1.728*S**(3.0/2.0)-0.0996*S**2.0)/T &
	+148.0248+137.1942*S**0.5+1.62142*S &
	-(24.4344+25.085*S**0.5+0.2474*S)*log(T) &
	+0.053105*S**0.5*T)

	! Kw (Millero 1995)
	C_consts(n,iKw)=exp( &
	148.96502-13847.26/T-23.6521*log(T) &
	+(118.67/T-5.977+1.0495*log(T))*S**0.5-0.01615*S)

	!KSi (Millero 1995)
	I=(19.924*S)/(1000.0-1.005*S)
	!I=0.02*S

	C_consts(n,iKSi)=exp( &
	-8904.2/T+117.385-19.334*log(T) &
	+((-458.79/T+3.5913)*(I**(-0.5)) &
	+(188.74/T-1.5998))*I &
	+(-12.1652/T+0.07871)*(I**2) &
	+log(1.0-0.001005*S))

	! KP1 (Millero 1995)
	C_consts(n,iKp1)=exp( &
	-4576.752/T+115.525-18.453*log(T) &
	+(-106.736/T+0.69171)*S**0.5 &
	+(-0.65643/T-0.01844)*S)

	! KP2 (Millero 1995)
	C_consts(n,iKp2)=exp( &
	-8814.715/T+172.0883-27.927*log(T) &
	+(-160.34/T+1.3566)*S**0.5 &
	+(0.37335/T-0.05778)*S)

	! KP3 (millero 1995)
	C_consts(n,iKp3)=exp( &
	-3070.75/T-18.141 &
	+(17.27039/T+2.81197)*S**0.5 &
	+(-44.99486/T-0.09984)*S)



enddo
! compare with DOE (1994)
!print*,'T',T
!print*,'S',S
!print*,'K1',log(C_consts(1,iK1)),(C_consts(1,iK1))
!print*,'K2',log(C_consts(1,iK2)),(C_consts(1,iK2))
!print*,'K0',log(C_consts(1,iK0)),(C_consts(1,iK0))
!print*,'Kw',log(C_consts(1,iKw)),(C_consts(1,iKw))
!print*,'KB',log(C_consts(1,iKB)),(C_consts(1,iKB))
!print*,'KSi',log(C_consts(1,iKSi)),(C_consts(1,iKSi))
!print*,'Kp1',log(C_consts(1,iKp1)),(C_consts(1,iKp1))
!print*,'Kp2',log(C_consts(1,iKp2)),(C_consts(1,iKp2))
!print*,'Kp3',log(C_consts(1,iKp3)),(C_consts(1,iKp3))

end subroutine calc_C_consts

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine calc_pCO2()
! ---------------------------------------------------------------------------------------!
!! Solve carbonate system for pC02
!! - M. Follows, T. Ito, S. Dutkiewicz (2006) in Ocean Modelling
! ---------------------------------------------------------------------------------------!

! local variables
real::pt,sit,ta,pco2,dic,H,bt,k1,k2,k1p,k2p,k3p,kb,kw,ksi,k0 ! ff in original but not used, added k0
real::gamm,co2s,hg,cag,bohg,h3po4g,h2po4g,hpo4g,po4g,siooh3g,denom,dummy,fg
integer::n,n_loop
real::h_guess_1,h_guess_2

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
dic=tracers_1(n,ioDIC)*r_rho ! mol m-3 -> mol kg-1
ta=tracers_1(n,ioALK)*r_rho ! mol m-3 -> mol kg-1
pt=tracers_1(n,ioPO4)*r_rho ! mol m-3 -> mol kg-1
!sit=silica_dt(n)
sit=7.5/1.0e6 ! global mean surface (Orr et al., 2017) umol kg-1 -> mol kg-1
bt=0.0004106*(S_dt(n)/35.0) ! total boron conc. from ZW2001 (mol kg-1)
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

if(C(n,ioH).eq.-1.0)then
	hg=10e-8 ! cold start
	h_guess_1=1.0
	h_guess_2=hg
	else
	hg=C(n,ioH) ! previous timestep
	h_guess_1=1.0
	h_guess_2=hg
endif


n_loop=0
do while (abs(h_guess_1-h_guess_2)>carbchem_tol)

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

hg=H
h_guess_1=h_guess_2
h_guess_2=hg
n_loop=n_loop+1
if(n_loop.eq.10) exit ! get-out clause

end do

! evaluate [CO2*]
co2s = dic/(1.0 + (k1/H) + (k1*k2/(H*H)))
!evaluate surface pCO2
pco2=dic/k0*((1.0+k1/H+(k1*k2)/(H*H))**(-1.0))
C(n,ioCO2) = co2s*rho ! mol kg-1 -> mol m-3
C(n,ioH) = H
C(n,iopCO2) = pco2

enddo

end subroutine calc_pCO2

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine calc_gasexchange()
! ---------------------------------------------------------------------------------------!
!! Calculate gas exchange between atmosphere and ocean
!! - Orr et al., (2017) Geoscientific Model Development
! ---------------------------------------------------------------------------------------!

real::loc_T,loc_T2,loc_T3,loc_T4,loc_Tr100,loc_T100_2,loc_TK,loc_S,loc_T100
REAL::Sc,Bunsen,Sol,gasex
real::kw,CO2star,CO2starair
real,dimension(n_ATM_tracers)::atm_dt
integer::n

atm_dt=ATM! copy ATM so can integrate main tracer below

do n=1,n_surface_boxes

	loc_T=T_dt(n)
	loc_S=S_dt(n)
	if(loc_T<-2.0)loc_T=-2.0
	if(loc_T>40.0)loc_T=40.0
	loc_T2=loc_T*loc_T
	loc_T3=loc_T2*loc_T
	loc_T4=loc_T3*loc_T
	loc_TK=loc_T+273.15
	loc_Tr100=100.0/loc_TK
	loc_T100=loc_TK/100.0
	loc_T100_2=loc_T100*loc_T100

	if(bg_C_select)then

		Sc=Sc_coeffs(1,iaCO2) &
		+Sc_coeffs(2,iaCO2)*loc_T &
		+Sc_coeffs(3,iaCO2)*loc_T2 &
		+Sc_coeffs(4,iaCO2)*loc_T3 &
		+Sc_coeffs(5,iaCO2)*loc_T4

		Sol=exp( &
		Sol_Orr(1,iaCO2) &
		+Sol_Orr(2,iaCO2)*loc_Tr100 &
		+Sol_Orr(3,iaCO2)*log(loc_T100) &
		+Sol_Orr(4,iaCO2)*loc_T100_2 &
		+loc_S* &
		(Sol_Orr(5,iaCO2) &
		+Sol_Orr(6,iaCO2)*loc_T100 &
		+Sol_Orr(7,iaCO2)*loc_T100_2)) ! mol L-3 atm-1
		Sol=Sol*1000.0 ! mol L-1 atm-1 -> mol m-3 atm-1

		kw=wind_dt(n)*((Sc/660.0)**(-0.5)) ! in m yr-1 (conversion factors and parameter a computed in wind_dt)

		CO2starair=Sol*ATM(iaCO2) ! [CO2*]sat (mol m-3)

		CO2star=C(n,ioCO2) ! [CO2*] (mol m-3)

		gasex=kw*(CO2starair-CO2star)*(1.0/50.0) ! air -> sea, mol m-3 yr-1 (n.b. hard coded depth in m!!)

		J(n,ioDIC)=J(n,ioDIC)+gasex ! update ocean source/sink
		Jatm(n,iaCO2)=Jatm(n,iaCO2)-(gasex*tm_vol(n)) ! update atmosphere (mol yr-1)

		diag(n,2)=gasex*tm_vol(n) ! mol yr-1
		diag(n,1)=gasex*50.0 ! mol m-2 yr-1
		! if(n.eq.2000)then
		! print*,'n',n
		! print*,'co2starair',CO2starair
		! print*,'co2star',CO2star
		! print*,'wind',wind_dt(n)
		! print*,'Sc',Sc,((Sc/660.0)**(-0.5))
		! print*,'T',T_dt(n)
		! print*,'S',S_dt(n)
		! print*,'dic',tracers_1(n,ioDIC) ! mol m-3 -> mol kg-1
		! print*,'ta',tracers_1(n,ioALK) ! mol m-3 -> mol kg-1
		! print*,'phos',tracers_1(n,ioPO4)! mol m-3 -> mol kg-1
		! print*,'Si',(7.5/1.0e6)*rho
		! print*,'CO2atm',ATM(iaCO2)
		! print*,'gasex',kw*(CO2starair-0.17285),86.399*((Sc/660.0)**(-0.5))*(CO2starair-0.17285)
		! stop
		! endif


	endif

enddo

end subroutine calc_gasexchange

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine restore_atm_CO2()
! ---------------------------------------------------------------------------------------!
!! Applies a restoring forcing for CO2
! ---------------------------------------------------------------------------------------!

if(bg_restore_atm_CO2) ATM(iaCO2)=bg_restore_atm_CO2_target*1.0e-6

end subroutine restore_atm_CO2

! ---------------------------------------------------------------------------------------!

! ---------------------------------------------------------------------------------------!

subroutine water_column()
! ---------------------------------------------------------------------------------------!
!! Implicit remineralisation of particulate organic carbon and CaCO3
! ---------------------------------------------------------------------------------------!

integer::n,count,nn
integer,dimension(maxval(tm_wc))::loc_wc_start,loc_wc_end
integer,dimension(tm_nbox)::loc_wc
real,dimension(tm_nbox)::loc_particles,loc_vol,remin,loc_depth_btm,profile,loc_particles_copy,loc_b
real,dimension(tm_nbox)::loc_rvol
!real::loc_poc,loc_remin_tot,loc_poc_copy
real::layerratio,frac,pom_above
real::start,finish

! ************* Convert grid arrays to water column order ********************** !
loc_wc=amul(Aconv,real(tm_wc))
loc_vol=amul(Aconv,tm_vol)
loc_depth_btm=amul(Aconv,tm_depth_btm)
loc_b=amul(Aconv,bg_martin_b)

loc_rvol=1.0/loc_vol

! find water column starting/end points
loc_wc_start(1)=1
count=2
do n=2,tm_nbox
	if(loc_wc(n).gt.loc_wc(n-1))then
		loc_wc_start(count)=n
		loc_wc_end(count-1)=n-1
		count=count+1
	endif
enddo
loc_wc_end(maxval(tm_wc))=tm_nbox


! ************* POM Remin ********************** !
! pre-calculate curve
! loc_particles=amul(Aconv,particles(:,isPOP))
! remin=(loc_depth_btm/120.0)**(bg_martin_b)
! loc_particles=loc_particles*loc_vol ! mol m-3 -> mol

! new
loc_particles=amul(Aconv,particles(:,isPOP))
profile=(loc_depth_btm/120.0)**(loc_b)
loc_particles_copy=loc_particles
remin(:)=0.0
do n=1,maxval(loc_wc)
    if(loc_wc_end(n)-loc_wc_start(n)==0)then ! no water column below
        nn=loc_wc_end(n)
        remin(nn)=loc_particles(nn)
    elseif(loc_wc_end(n)-loc_wc_start(n)==1)then ! no water column below
        nn=loc_wc_end(n)
        layerratio=loc_vol(nn-1)*loc_rvol(nn)
				pom_above=loc_particles(nn-1)*layerratio
        remin(nn)=loc_particles(nn)+pom_above
    else
        do nn=loc_wc_start(n)+1,loc_wc_end(n) ! loop over water column
            layerratio=loc_vol(nn-1)*loc_rvol(nn)
						pom_above=loc_particles(nn-1)*layerratio
            if(nn<=loc_wc_start(n)+1)then ! surface
                loc_particles(nn)=loc_particles(nn)+pom_above! add particles from layer above
                loc_particles(nn-1)=0.0
            elseif(nn>loc_wc_start(n)+1 .and. nn<loc_wc_end(n))then ! interior
                frac=(profile(nn)/profile(nn-1))
                loc_particles(nn)=pom_above*frac
                remin(nn)=pom_above-loc_particles(nn)
            else
                loc_particles(nn)=0.0
                remin(nn)=pom_above-loc_particles(nn)
            endif
        enddo
    endif
    !print*,n,loc_wc_end(n)-loc_wc_start(n),sum(remin(loc_wc_start(n):loc_wc_end(n))* &
		!loc_vol(loc_wc_start(n):loc_wc_end(n))) &
		!/sum(loc_vol(loc_wc_start(n):loc_wc_end(n))) - &
    !sum(loc_particles_copy(loc_wc_start(n):loc_wc_end(n))*loc_vol(loc_wc_start(n):loc_wc_end(n))) &
		!/sum(loc_vol(loc_wc_start(n):loc_wc_end(n)))
enddo
!print*,n,sum(remin*loc_vol)/sum(loc_vol)-sum(loc_particles_copy*loc_vol)/sum(loc_vol)

remin=amul_transpose(Aconv,remin)
J(:,ioPO4)=J(:,ioPO4)+remin ! POP remineralisation
diag(:,4)=remin
if(bg_C_select)then
	! OM remineralisation
	J(:,ioDIC)=J(:,ioDIC)+remin*bg_C_to_P
	J(:,ioALK)=J(:,ioALK)-remin*bg_N_to_P
endif

! ************* CaCO3 Dissolution ********************** !
loc_particles=amul(Aconv,particles(:,isCaCO3))
profile=exp((120.0-loc_depth_btm)/bg_CaCO3_length_scale)
loc_particles_copy=loc_particles
remin(:)=0.0
do n=1,maxval(loc_wc)
    if(loc_wc_end(n)-loc_wc_start(n)==0)then ! no water column below
        nn=loc_wc_end(n)
        remin(nn)=loc_particles(nn)
    elseif(loc_wc_end(n)-loc_wc_start(n)==1)then ! no water column below
        nn=loc_wc_end(n)
        layerratio=loc_vol(nn-1)*loc_rvol(nn)
				pom_above=loc_particles(nn-1)*layerratio
        remin(nn)=loc_particles(nn)+pom_above
    else
        do nn=loc_wc_start(n)+1,loc_wc_end(n) ! loop over water column
            layerratio=loc_vol(nn-1)*loc_rvol(nn)
						pom_above=loc_particles(nn-1)*layerratio
            if(nn<=loc_wc_start(n)+1)then ! surface
                loc_particles(nn)=loc_particles(nn)+pom_above! add particles from layer above
                loc_particles(nn-1)=0.0
            elseif(nn>loc_wc_start(n)+1 .and. nn<loc_wc_end(n))then ! interior
                frac=(profile(nn)/profile(nn-1))
                loc_particles(nn)=pom_above*frac
                remin(nn)=pom_above-loc_particles(nn)
            else
                loc_particles(nn)=0.0
                remin(nn)=pom_above-loc_particles(nn)
            endif
        enddo
    endif
    !print*,n,loc_wc_end(n)-loc_wc_start(n),sum(remin(loc_wc_start(n):loc_wc_end(n))* &
		!loc_vol(loc_wc_start(n):loc_wc_end(n))) &
		!/sum(loc_vol(loc_wc_start(n):loc_wc_end(n))) - &
    !sum(loc_particles_copy(loc_wc_start(n):loc_wc_end(n))*loc_vol(loc_wc_start(n):loc_wc_end(n))) &
		!/sum(loc_vol(loc_wc_start(n):loc_wc_end(n)))
enddo

remin=amul_transpose(Aconv,remin)
J(:,ioDIC)=J(:,ioDIC)+remin ! CaCO3 dissolution
J(:,ioALK)=J(:,ioALK)+remin*2.0 ! CO32-
diag(:,5)=remin


! loc_particles=amul(Aconv,particles(:,isCaCO3))
! loc_particles=loc_particles*loc_vol ! mol m-3 -> mol
!
! remin=exp((120.0-loc_depth_btm)/bg_CaCO3_length_scale)
!
! do n=1,maxval(loc_wc)
! 	loc_remin_tot=0.0
! 	loc_poc=0.0
! 	loc_poc_copy=0.0
! 	count=1
! 	do nn=loc_wc_start(n),loc_wc_end(n)
! 		if(count.le.bg_n_euphotic_lyrs)then ! surface
! 			loc_poc=loc_poc+loc_particles(nn)! integrate POM in surface layers (mol)
! 			loc_particles(nn)=0.0 ! set flux to zero
! 			count=count+1
! 		elseif(count.gt.bg_n_euphotic_lyrs)then ! interior
! 			loc_particles(nn)=(remin(nn-1)-remin(nn))*loc_poc
! 			loc_remin_tot=loc_remin_tot+(remin(nn-1)-remin(nn))
! 			count=count+1
! 		endif
!
! 		! once reached the end of the wc, remin all remaining pom
! 		if(nn.eq.loc_wc_end(n)) loc_particles(nn)=loc_particles(nn)+(1.0-loc_remin_tot)*loc_poc
!
! 	enddo
!
! enddo
!
! loc_particles=loc_particles/loc_vol ! mol -> mol m-3
! ! reorder particles array
! particles(:,isCaCO3)=amul_transpose(Aconv,loc_particles)
!
! ! ************* Update other tracers ********************** !
!
! ! reorder particles array
! ! J(:,ioPO4)=J(:,ioPO4)+particles(:,isPOP) ! POP remineralisation
! ! diag(:,4)=particles(:,isPOP)
! ! diag(:,5)=particles(:,isCaCo3)
! J(:,ioPO4)=J(:,ioPO4)+particles(:,isPOP) ! POP remineralisation
! diag(:,4)=particles(:,isPOP)
! diag(:,5)=particles(:,isCaCo3)
!
! if(bg_C_select)then
! 	! OM remineralisation
! 	J(:,ioDIC)=J(:,ioDIC)+(particles(:,isPOP)*bg_C_to_P)
! 	J(:,ioALK)=J(:,ioALK)-(particles(:,isPOP)*bg_N_to_P)
!
! 	! CaCO3 dissolution
! 	J(:,ioDIC)=J(:,ioDIC)+(particles(:,isCaCO3))
! 	J(:,ioALK)=J(:,ioALK)+(particles(:,isCaCO3)*2.0)
! endif


end subroutine water_column


end module bg_module
