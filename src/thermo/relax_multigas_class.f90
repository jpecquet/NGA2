!> Phase-change relaxation of a liquid and ONE of its vapors in a multi-species ideal-gas mixture.
!> Wraps a vapor/carrier relaxation model of the relax_igmix_sg family (e.g. relax_igmix_nasg_max_ent_finite), which solves
!> for a gas of exactly two species. During one relaxation only the vapor mass changes: every other gas species keeps its
!> mass, so together they form a carrier of frozen composition. Per cell, the wrapper builds that carrier as the second
!> species of a two-species igmix (mass-weighted cv, cp, q, qp; R=cp-cv, gamma=cp/cv; exact for an igmix), calls the inner
!> model on the vapor/carrier cell vector, and maps the new liquid/gas masses, energies, VF and vapor mass back.
!> Cell vector layout (amrmpcomp): Q(1:4)=[m_L,m_G,E_L,E_G], Q(5:7) unused, liquid species Q(8:6+liq%ns) (none: liq%ns=1),
!> then the first gas%ns-1 gas partial densities. The vapor must be one of them (not the implicit last species).
module relax_multigas_class
   use precision,            only: WP
   use thermorelax_class,    only: thermorelax,RELAX_OK,RELAX_DEGENERATE
   use relax_igmix_sg_class, only: relax_igmix_sg
   use stiffened_gas_class,  only: stiffened_gas
   use igmix_class,          only: igmix
   implicit none
   private
   public :: relax_multigas

   type, extends(thermorelax) :: relax_multigas
      class(igmix),          pointer :: gas=>null()     !< Full gas mixture
      class(stiffened_gas),  pointer :: liq=>null()     !< Liquid (single species)
      class(relax_igmix_sg), pointer :: inner=>null()   !< Vapor/carrier model, initialized on gas2 by this wrapper
      type(igmix),           pointer :: gas2=>null()    !< Two-species scratch gas [vapor,carrier], carrier reset per cell
      integer :: iv=0    !< Vapor index in the full gas
      integer :: iG0=0   !< Index in Q of the first gas partial density
   contains
      procedure :: initialize
      procedure :: apply
   end type relax_multigas

contains

   !> Store the materials, build the vapor/carrier scratch gas, and initialize the inner model on it (indV=1, indA=2)
   subroutine initialize(this,liq,gas,iv,inner)
      use messager, only: die
      implicit none
      class(relax_multigas), intent(inout) :: this
      class(stiffened_gas),  target, intent(in) :: liq
      class(igmix),          target, intent(in) :: gas
      integer,               intent(in) :: iv
      class(relax_igmix_sg), target, intent(inout) :: inner
      integer :: c
      if (liq%ns.ne.1) call die('[relax_multigas initialize] the liquid must have a single species')
      if (iv.lt.1.or.iv.gt.gas%ns-1) call die('[relax_multigas initialize] the vapor must be a transported gas species (1..ns-1)')
      this%liq=>liq; this%gas=>gas; this%iv=iv; this%iG0=7+liq%ns
      ! Scratch gas: vapor, and the carrier initialized as the last species (reset per cell)
      c=gas%ns
      allocate(this%gas2)
      call this%gas2%initialize(gamma=[gas%gamma(iv),gas%gamma(c)],cv=[gas%cv(iv),gas%cv(c)],q=[gas%q(iv),gas%q(c)], &
      &                         qp=[gas%qp(iv),gas%qp(c)],species_names=['vapor  ','carrier'],name='vapor+carrier')
      call inner%initialize(liq=liq,gas=this%gas2,indV=1,indA=2)
      this%inner=>inner
   end subroutine initialize

   !> Relax the cell: pseudo-binary vapor/carrier problem solved by the inner model, then mapped back
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_multigas),  intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP), dimension(this%gas%ns) :: y,w
      real(WP), dimension(this%iG0) :: Qin
      real(WP) :: yc,cvc,cpc,qc,qpc
      integer :: ns,k,ivQ
      ns=this%gas%ns; ivQ=this%iG0+this%iv-1
      ! Full gas composition (last species 1-sum)
      y=0.0_WP; y(ns)=1.0_WP
      if (Q(2).gt.0.0_WP) then
         y(1:ns-1)=max(0.0_WP,min(Q(this%iG0:this%iG0+ns-2)/Q(2),1.0_WP))
         y(ns)=max(0.0_WP,1.0_WP-sum(y(1:ns-1)))
      end if
      ! Carrier = all species but the vapor, with their current proportions (the last species alone if there is none)
      w=y; w(this%iv)=0.0_WP; yc=sum(w)
      if (yc.gt.1.0e-14_WP) then
         w=w/yc
      else
         w=0.0_WP; w(ns)=1.0_WP
      end if
      cvc=sum(w*this%gas%cv); cpc=sum(w*this%gas%cp); qc=sum(w*this%gas%q); qpc=sum(w*this%gas%qp)
      this%gas2%cv(2)=cvc; this%gas2%cp(2)=cpc; this%gas2%R(2)=cpc-cvc; this%gas2%gamma(2)=cpc/cvc
      this%gas2%q(2)=qc; this%gas2%qp(2)=qpc
      this%inner%Ma=1.0_WP/this%gas2%R(2)
      ! Vapor/carrier cell vector for the inner model: Q(1:7) and the vapor partial density
      Qin(1:this%iG0-1)=Q(1:this%iG0-1)
      Qin(this%iG0)=Q(ivQ)
      call this%inner%apply(dt=dt,VF=VF,Q=Qin,Pjump=Pjump,ierr=ierr)
      ! Map back: masses and energies of the phases, and the vapor; the other gas species keep their masses
      Q(1:4)=Qin(1:4)
      Q(ivQ)=Qin(this%iG0)
   end subroutine apply

end module relax_multigas_class
