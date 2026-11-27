!> Thermodynamic relaxation for an ideal-gas MIXTURE / Noble-Abel-stiffened-gas pair at frozen composition
!> Same algorithm as relax_ig_nasg (p_relax, pT_relax, PThybrid), with the gas gamma, cv and formation energy q taken
!> per cell from the gas composition: at fixed composition an igmix is an ideal gas with gamma_m=cp_m/cv_m, cv_m and
!> q_m=sum(y*q). No mass is transferred, so the composition does not change during the relaxation.
!> Cell vector layout (amrmpcomp): Q(1:4)=[m_L,m_G,E_L,E_G], Q(5:7) unused, then the liquid species Q(8:6+liq%ns) and the
!> gas species Q(7+liq%ns:5+liq%ns+gas%ns): the first ns-1 partial densities of each phase, the last species is 1-sum.
module relax_igmix_nasg_frozen_class
   use precision,         only: WP
   use thermorelax_class, only: thermorelax,RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS,RELAX_VACUUM_GAS,RELAX_DEGENERATE
   use igmix_class,       only: igmix
   use nasg_class,        only: nasg
   implicit none
   private

   public :: relax_igmix_nasg_frozen
   public :: Prelax,PTrelax,PThybrid

   integer, parameter :: Prelax  =1   !< Mechanical relaxation only
   integer, parameter :: PTrelax =2   !< Mechanical + thermal relaxation
   integer, parameter :: PThybrid=3   !< Mechanical below Tratmax temperature contrast, mechanical+thermal above

   type, extends(thermorelax) :: relax_igmix_nasg_frozen
      class(igmix), pointer :: gas=>null()
      class(nasg),  pointer :: liq=>null()
      integer  :: model  =Prelax
      real(WP) :: RHOGmin=1.0e-2_WP   !< Skip mechanical relax when gas density falls below this
      real(WP) :: VFratmax=10.0_WP    !< Max per-call phase-volume change factor in p_relax (partial relax beyond)
      real(WP) :: Tratmax=10.0_WP     !< PThybrid: temperature contrast max(TG/TL,TL/TG) above which pT_relax is used
      real(WP) :: phist  =1.0_WP      !< Temporal weighting on equilibrium pressure
      real(WP) :: phi0   =0.0_WP      !< Temporal weighting on interface pressure
      integer  :: iG0=0               !< Index in Q of the first gas species partial density
   contains
      procedure :: initialize
      procedure :: apply
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: gas_y
   end type relax_igmix_nasg_frozen

contains

   !> Initialize: store EOS pointers and locate the gas species in the cell vector
   subroutine initialize(this,gas,liq)
      implicit none
      class(relax_igmix_nasg_frozen), intent(inout) :: this
      class(igmix), target, intent(in) :: gas
      class(nasg),  target, intent(in) :: liq
      this%gas=>gas; this%liq=>liq
      this%iG0=7+liq%ns
   end subroutine initialize

   !> Gas composition of the cell vector Q (last species 1-sum, clipped to [0,1])
   function gas_y(this,Q) result(y)
      implicit none
      class(relax_igmix_nasg_frozen), intent(in) :: this
      real(WP), dimension(:), intent(in) :: Q
      real(WP), dimension(this%gas%ns) :: y
      integer :: ns
      ns=this%gas%ns
      y=0.0_WP; y(ns)=1.0_WP
      if (ns.eq.1.or.Q(2).le.0.0_WP) return
      y(1:ns-1)=max(0.0_WP,min(Q(this%iG0:this%iG0+ns-2)/Q(2),1.0_WP))
      y(ns)=max(0.0_WP,1.0_WP-sum(y(1:ns-1)))
   end function gas_y

   !> Apply: gate on mixture cells, dispatch via model
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      use messager, only: die
      implicit none
      class(relax_igmix_nasg_frozen), intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP) :: TL,TG
      ! Mixture cells only
      if (VF.le.0.0_WP.or.VF.ge.1.0_WP) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Dispatch
      select case (this%model)
      case (Prelax);  call this%p_relax (dt,VF,Q,Pjump,ierr)
      case (PTrelax); call this%pT_relax(dt,VF,Q,Pjump,ierr)
      case (PThybrid)
         if (any(Q(1:4).le.0.0_WP)) then
            call this%p_relax(dt,VF,Q,Pjump,ierr)
         else
            TL=this%liq%get_T_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP])
            TG=this%gas%get_T_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=this%gas_y(Q))
            if (TL.gt.0.0_WP.and.TG.gt.0.0_WP.and.max(TG/TL,TL/TG).gt.this%Tratmax) then
               call this%pT_relax(dt,VF,Q,Pjump,ierr)
            else
               call this%p_relax(dt,VF,Q,Pjump,ierr)
            end if
         end if
      case default; call die('[relax_igmix_nasg_frozen apply] unknown model')
      end select
   end subroutine apply

   !> Mechanical relaxation (relax_ig_nasg p_relax with the cell's gas gamma and q)
   !> Conserves phasic masses Q(1:2), the species, and total internal energy Q(3)+Q(4)
   subroutine p_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_frozen), intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP) :: PG,PL,ZG,ZL,Pint,cJ,bL,rhoL,gG,qG
      real(WP) :: a,b,d,n1,n0,d1,d0,Peq,VFeq
      real(WP), dimension(this%gas%ns) :: y
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.this%RHOGmin) then; if (present(ierr)) ierr=RELAX_VACUUM_GAS; return; end if
      ! Gas mixture constants at the cell composition
      y=this%gas_y(Q)
      gG=sum(y*this%gas%cp)/sum(y*this%gas%cv)
      qG=sum(y*this%gas%q)
      bL=this%liq%b
      rhoL=Q(1)/(VF)
      ! Phasic pressures
      PL=this%liq%get_p_from_rho_e(rho=rhoL,e=Q(3)/Q(1),y=[1.0_WP])
      PG=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
      if (PL.le.-this%liq%pinf.or.1.0_WP-bL*rhoL.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (PG.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      ! Phasic acoustic impedances (rho*c)
      ZL=rhoL*this%liq%get_c_from_p_rho(p=PL,rho=rhoL,y=[1.0_WP])
      ZG=Q(2)/(1.0_WP-VF)*this%gas%get_c_from_p_rho(p=PG,rho=Q(2)/(1.0_WP-VF),y=y)
      cJ=ZL/(ZG+ZL)
      ! Interface pressure model
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Quadratic setup (pinf_g = 0; liquid co-volume b enters n1,n0)
      n1=VF*this%phist+Q(1)*bL/(this%liq%gamma-1.0_WP)
      n0=VF*(this%phi0*Pint-this%phist*cJ*Pjump)+Q(3)-Q(1)*this%liq%q+Q(1)*bL*this%liq%gamma/(this%liq%gamma-1.0_WP)*this%liq%pinf
      d1=this%phist+1.0_WP/(this%liq%gamma-1.0_WP)
      d0=this%phi0*Pint-this%phist*cJ*Pjump+this%liq%gamma/(this%liq%gamma-1.0_WP)*this%liq%pinf
      a=d1*(1.0_WP/(gG-1.0_WP)+this%phist*VF)+n1*(-1.0_WP/(gG-1.0_WP)-this%phist)
      b=d1*(-Pjump/(gG-1.0_WP)+Q(2)*qG-Q(4)+VF*(this%phi0*Pint-this%phist*cJ*Pjump))+n1*(Pjump/(gG-1.0_WP)-this%phi0*Pint+this%phist*cJ*Pjump)+d0*(1.0_WP/(gG-1.0_WP)+this%phist*VF)+n0*(-1.0_WP/(gG-1.0_WP)-this%phist)
      d=d0*(-Pjump/(gG-1.0_WP)+Q(2)*qG-Q(4)+VF*(this%phi0*Pint-this%phist*cJ*Pjump))+n0*(Pjump/(gG-1.0_WP)-this%phi0*Pint+this%phist*cJ*Pjump)
      ! Equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Soundness check (pinf_g = 0)
      if (Peq.le.-this%liq%pinf) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (Peq-Pjump.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      ! Equilibrium VF (strict bounds)
      VFeq=(n1*Peq+n0)/(d1*Peq+d0)
      if (VFeq.lt.0.0_WP.or.VFeq.gt.1.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Clamp per-call phase-volume change to a factor VFratmax (bounds consistent for VFratmax>=1)
      VFeq=max(1.0_WP-this%VFratmax*(1.0_WP-VF),VF/this%VFratmax, &
      &        min(VFeq,1.0_WP-(1.0_WP-VF)/this%VFratmax,this%VFratmax*VF))
      ! Update Q with p*dV work at the relaxed interface pressure
      Q(3)=Q(3)-(this%phist*Peq+this%phi0*Pint-this%phist*cJ*Pjump)*(VFeq-VF)
      Q(4)=Q(4)+(this%phist*Peq+this%phi0*Pint-this%phist*cJ*Pjump)*(VFeq-VF)
      VF=VFeq
      if (present(ierr)) ierr=RELAX_OK
   end subroutine p_relax

   !> Thermal relaxation (relax_ig_nasg pT_relax with the cell's gas gamma, cv and q)
   !> Conserves phasic masses Q(1:2), the species, and total internal energy Q(3)+Q(4)
   subroutine pT_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_frozen), intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP) :: a,b,d,Peq,VFeq,Eth
      real(WP) :: cv1,cv2,g1,g2,q2,pinf,R1,R2,bL,ombm
      real(WP), dimension(this%gas%ns) :: y
      ! Mechanical relaxation first (its own ierr is not propagated; pT owns the final verdict)
      call this%p_relax(dt,VF,Q,Pjump)
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Shorthands, the gas ones at the cell composition
      y=this%gas_y(Q)
      cv1=this%liq%cv; cv2=sum(y*this%gas%cv)
      g1=this%liq%gamma; g2=sum(y*this%gas%cp)/cv2; q2=sum(y*this%gas%q); pinf=this%liq%pinf
      R1=cv1*(g1-1.0_WP); R2=cv2*(g2-1.0_WP)
      bL=this%liq%b
      ombm=1.0_WP-bL*Q(1)
      if (ombm.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      ! Thermal internal energy (formation energies removed); invariant under thermal relax
      Eth=Q(3)+Q(4)-Q(1)*this%liq%q-Q(2)*q2
      ! Quadratic for liquid equilibrium pressure Peq (gas pressure Peq-Pjump), TL=TG, pinf_g=0
      a=ombm*(Q(1)*cv1+Q(2)*cv2)
      b=ombm*(Q(1)*cv1*g1*pinf+Q(2)*cv2*pinf-Pjump*(Q(1)*cv1+Q(2)*cv2))-Eth*(Q(1)*R1+Q(2)*R2)
      d=-ombm*Pjump*(Q(1)*cv1*g1*pinf+Q(2)*cv2*pinf)+Pjump*Eth*Q(1)*R1-Eth*Q(2)*R2*pinf
      if (b**2-4.0_WP*a*d.lt.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      if (Peq.le.-this%liq%pinf) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (Peq-Pjump.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      VFeq=bL*Q(1)+ombm*Q(1)*R1*(Peq-Pjump)/(Q(1)*R1*(Peq-Pjump)+Q(2)*R2*(Peq+pinf))
      if (VFeq.lt.0.0_WP.or.VFeq.gt.1.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      Q(3)=(VFeq-bL*Q(1))*(Peq+g1*pinf)/(g1-1.0_WP)+Q(1)*this%liq%q
      Q(4)=(1.0_WP-VFeq)*(Peq-Pjump  )/(g2-1.0_WP)+Q(2)*q2
      VF=VFeq
      if (present(ierr)) ierr=RELAX_OK
   end subroutine pT_relax

end module relax_igmix_nasg_frozen_class
