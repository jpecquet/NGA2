!> One-step irreversible Arrhenius reaction R -> P, first order in the reactant mass fraction:
!>    dY_R/dt = -k(T) Y_R,   k(T) = A exp(-Ta/T) for T >= Tign, 0 below (cold-boundary cut-off)
!> Advanced exactly over dt at fixed T (exponential integrator), so Y_R stays in [0,1] for any dt.
!> The heat release Q is carried by the mixture EOS as q_R-q_P (igmix formation energies).
module onestep_arrhenius_class
   use precision,     only: WP
   use reaction_class, only: reaction
   implicit none
   private
   public :: onestep_arrhenius

   type, extends(reaction) :: onestep_arrhenius
      real(WP) :: A=0.0_WP      !< Pre-exponential factor (1/time)
      real(WP) :: Ta=0.0_WP     !< Activation temperature Ea/R
      real(WP) :: Tign=0.0_WP   !< Ignition temperature: no reaction below
      integer  :: iR=1,iP=2     !< Reactant and product species indices
   contains
      procedure :: initialize
      procedure :: advance
      procedure :: get_rate
      procedure :: get_k
   end type onestep_arrhenius

contains

   !> Initialize with the mixture species count, the rate constants, and the reactant/product indices
   subroutine initialize(this,ns,A,Ta,Tign,iR,iP)
      use messager, only: die
      implicit none
      class(onestep_arrhenius), intent(inout) :: this
      integer,  intent(in) :: ns
      real(WP), intent(in) :: A,Ta
      real(WP), intent(in), optional :: Tign
      integer,  intent(in), optional :: iR,iP
      this%ns=ns; this%A=A; this%Ta=Ta
      this%Tign=0.0_WP; if (present(Tign)) this%Tign=Tign
      this%iR=1; if (present(iR)) this%iR=iR
      this%iP=2; if (present(iP)) this%iP=iP
      if (this%iR.lt.1.or.this%iR.gt.ns.or.this%iP.lt.1.or.this%iP.gt.ns.or.this%iR.eq.this%iP) &
      &  call die('[onestep_arrhenius initialize] bad reactant/product indices')
   end subroutine initialize

   !> Rate constant k(T)
   pure real(WP) function get_k(this,T) result(k)
      implicit none
      class(onestep_arrhenius), intent(in) :: this
      real(WP), intent(in) :: T
      k=0.0_WP
      if (T.ge.this%Tign.and.T.gt.0.0_WP) k=this%A*exp(-this%Ta/T)
   end function get_k

   !> Y_R <- Y_R exp(-k dt), the converted mass goes to the product
   subroutine advance(this,rho,T,y,dt)
      implicit none
      class(onestep_arrhenius), intent(in) :: this
      real(WP), intent(in) :: rho,T,dt
      real(WP), dimension(:), intent(inout) :: y
      real(WP) :: dY
      dY=max(y(this%iR),0.0_WP)*(1.0_WP-exp(-this%get_k(T)*dt))
      y(this%iR)=y(this%iR)-dY
      y(this%iP)=y(this%iP)+dY
   end subroutine advance

   !> Inverse chemical time: k(T) where there is reactant
   real(WP) function get_rate(this,rho,T,y) result(r)
      implicit none
      class(onestep_arrhenius), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      r=0.0_WP
      if (y(this%iR).gt.1.0e-6_WP) r=this%get_k(T)
   end function get_rate

end module onestep_arrhenius_class
