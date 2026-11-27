!> Global one-step irreversible reaction  F + nu O -> (1+nu) P  (mass basis), Westbrook-Dryer rate law:
!>    omega = Apre exp(-Ta/T) [F]^aF [O]^bO   (mol of F per m^3 per s; concentrations [k]=rho*Y_k/M_k in mol/m^3)
!>    dY_F/dt = -omega M_F/rho,   dY_O/dt = nu dY_F/dt,   dY_P/dt = -(1+nu) dY_F/dt;  other species (inerts) unchanged.
!> Advanced over dt at fixed (rho,T) with midpoint sub-steps that consume at most a fraction 'frac' of the limiting
!> reactant each, clipped at its exhaustion, so mass fractions stay in [0,1] and sum to 1. The heat release is carried by
!> the mixture EOS (igmix formation energies q).
module global_onestep_class
   use precision,      only: WP
   use reaction_class, only: reaction
   implicit none
   private
   public :: global_onestep

   type, extends(reaction) :: global_onestep
      real(WP) :: Apre=0.0_WP    !< Pre-exponential factor in SI concentration units, (mol/m^3)^(1-a-b)/s
      real(WP) :: Ta=0.0_WP      !< Activation temperature Ea/Ru (K)
      real(WP) :: aF=1.0_WP,bO=1.0_WP  !< Reaction orders in [F] and [O]
      real(WP) :: Tign=0.0_WP    !< No reaction below this temperature
      real(WP) :: nu=1.0_WP      !< Mass of O consumed per mass of F
      real(WP) :: MF=1.0_WP,MO=1.0_WP   !< Molar masses of F and O (kg/mol)
      integer  :: iF=1,iO=2,iP=3 !< Species indices of fuel, oxidizer and (lumped) product
      real(WP) :: frac=0.1_WP    !< Largest fraction of the limiting reactant consumed per sub-step
      integer  :: nsub_max=2000  !< Sub-step cap (the last sub-step then consumes what it can)
   contains
      procedure :: initialize
      procedure :: advance
      procedure :: get_rate
      procedure :: omega
   end type global_onestep

contains

   subroutine initialize(this,ns,iF,iO,iP,nu,MF,MO,Apre,Ta,aF,bO,Tign)
      use messager, only: die
      implicit none
      class(global_onestep), intent(inout) :: this
      integer,  intent(in) :: ns,iF,iO,iP
      real(WP), intent(in) :: nu,MF,MO,Apre,Ta,aF,bO
      real(WP), intent(in), optional :: Tign
      this%ns=ns; this%iF=iF; this%iO=iO; this%iP=iP
      this%nu=nu; this%MF=MF; this%MO=MO; this%Apre=Apre; this%Ta=Ta; this%aF=aF; this%bO=bO
      this%Tign=0.0_WP; if (present(Tign)) this%Tign=Tign
      if (any([iF,iO,iP].lt.1).or.any([iF,iO,iP].gt.ns).or.iF.eq.iO.or.iF.eq.iP.or.iO.eq.iP) &
      &  call die('[global_onestep initialize] bad species indices')
   end subroutine initialize

   !> Fuel consumption rate (mol/m^3/s) at (rho,T,Y_F,Y_O)
   pure real(WP) function omega(this,rho,T,YF,YO) result(w)
      implicit none
      class(global_onestep), intent(in) :: this
      real(WP), intent(in) :: rho,T,YF,YO
      w=0.0_WP
      if (T.lt.this%Tign.or.T.le.0.0_WP.or.YF.le.0.0_WP.or.YO.le.0.0_WP) return
      w=this%Apre*exp(-this%Ta/T)*(rho*YF/this%MF)**this%aF*(rho*YO/this%MO)**this%bO
   end function omega

   subroutine advance(this,rho,T,y,dt)
      implicit none
      class(global_onestep), intent(in) :: this
      real(WP), intent(in) :: rho,T,dt
      real(WP), dimension(:), intent(inout) :: y
      real(WP) :: t_,h,r,lim,dYF,YF,YO
      integer :: n
      t_=0.0_WP; n=0
      do while (t_.lt.dt)
         YF=y(this%iF); YO=y(this%iO)
         lim=min(YF,YO/this%nu)
         if (lim.le.0.0_WP) exit
         r=this%omega(rho,T,YF,YO)*this%MF/rho          ! -dY_F/dt
         if (r.le.0.0_WP) exit
         n=n+1
         h=dt-t_
         if (n.lt.this%nsub_max) h=min(h,this%frac*lim/r)
         ! Midpoint rate over the sub-step
         dYF=min(0.5_WP*h*r,lim)
         r=this%omega(rho,T,YF-dYF,YO-this%nu*dYF)*this%MF/rho
         dYF=min(h*r,lim)
         y(this%iF)=y(this%iF)-dYF
         y(this%iO)=y(this%iO)-this%nu*dYF
         y(this%iP)=y(this%iP)+(1.0_WP+this%nu)*dYF
         t_=t_+h
      end do
      y(this%iF)=max(y(this%iF),0.0_WP); y(this%iO)=max(y(this%iO),0.0_WP)
   end subroutine advance

   !> Inverse time to consume the limiting reactant (floored at 1e-2 so that a nearly burnt cell does not limit dt)
   real(WP) function get_rate(this,rho,T,y) result(rr)
      implicit none
      class(global_onestep), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      rr=this%omega(rho,T,y(this%iF),y(this%iO))*this%MF/rho/max(min(y(this%iF),y(this%iO)/this%nu),1.0e-2_WP)
   end function get_rate

end module global_onestep_class
