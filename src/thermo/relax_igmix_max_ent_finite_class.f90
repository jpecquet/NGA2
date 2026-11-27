!> Finite-rate maximum-entropy relaxation (relaxation_ent_max.tex, Finite-rate relaxation): rate law dx~/dt=R*F^,
!> R=diag(lambda,vartheta,nu). The transfers faster than the finite one are instantaneous, the slower ones inactive;
!> a finite stage moves dx~=omega*K^-1*F^0, omega=1-exp(-r*K*dt) (eq. finite_solution), with K its pivot at the step start
!> (Pelanti 2022 eqs. 54-56, 65-67, 75-77).
module relax_igmix_max_ent_finite_class
   use precision,                 only: WP
   use relax_igmix_max_ent_inst_class, only: relax_igmix_max_ent_inst,phasic_state,max_ent_apply
   use relax_igmix_sg_class,      only: Prelax,PTrelax,PTgrelax
   use thermorelax_class,         only: RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS
   implicit none
   private

   public :: relax_igmix_max_ent_finite
   public :: restore_faster_eq,energy_split

   !> Closure after the finite step: restore the faster equilibria by their stages (Activation), or Pelanti's
   !> mixture-energy relation at the new V_l and masses (Pelanti 2022 eqs. 59, 70, 82)
   integer, parameter :: restore_faster_eq=1,energy_split=2

   type, extends(relax_igmix_max_ent_inst) :: relax_igmix_max_ent_finite
      real(WP) :: MRR =huge(1.0_WP)   !< Mechanical rate lambda [1/(Pa s)]; huge = instantaneous
      real(WP) :: TRR =huge(1.0_WP)   !< Thermal rate vartheta [W/(m^3 K)]; huge = instantaneous
      real(WP) :: CRR =huge(1.0_WP)   !< Chemical rate nu [kg s/m^5]; huge = instantaneous
      real(WP) :: tauC=-1.0_WP        !< Chemical relaxation time; >0 sets nu=1/(K^C**tauC), K^C at the step start, so omega=1-exp(-dt/tauC)
      !> Default restore_faster_eq: the faster stages end exactly at their equilibria. energy_split is Pelanti's one-pass closure
      !> (F^_M=0 only), whose second-order energy error stays in the phase of smaller heat capacity: T_l/=T_g grows with the transfer
      integer  :: closure=restore_faster_eq
      !> p_I: inherited default pI_end, p* for the instantaneous mechanical stage (NR_stage_pI_end); the finite mechanical step
      !> holds the impedance p_I of its start state, between the pressures (tex: Choice of p_I)
   contains
      procedure :: apply
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: pTg_relax
      procedure :: finite_step
      procedure :: mixture_energy_closure
   end type relax_igmix_max_ent_finite

contains

   !> Rate-separation checks, then the instantaneous class's apply on this, so that it dispatches to the finite workers below
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      use messager, only: die
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      logical :: nuc0
      if (this%model.ne.Prelax.and.this%MRR.lt.huge(1.0_WP)) call die('[relax_igmix_max_ent_finite apply] pT and pTg require MRR=huge')
      if (this%model.eq.PTgrelax.and.this%TRR.lt.huge(1.0_WP)) call die('[relax_igmix_max_ent_finite apply] pTg requires TRR=huge')
      ! No nucleation at zero chemical rate: the seed is itself a phase change
      nuc0=this%do_nucleate
      if (this%CRR.le.0.0_WP.and.this%tauC.le.0.0_WP) this%do_nucleate=.false.
      call max_ent_apply(this,dt,VF,Q,Pjump,ierr)
      this%do_nucleate=nuc0
   end subroutine apply

   !> Mechanical relaxation at rate lambda
   subroutine p_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      integer :: ier
      if (this%MRR.ge.huge(1.0_WP)) then
         call this%relax_M(VF,Q,Pjump,ier)
      else
         call this%finite_step(1,dt,VF,Q,Pjump,ier)
      end if
      if (present(ierr)) ierr=ier
   end subroutine p_relax

   !> Instantaneous mechanical stage, then thermal relaxation at rate vartheta
   subroutine pT_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      integer :: ier
      if (this%TRR.ge.huge(1.0_WP)) then
         call this%relax_MT(VF,Q,Pjump,ier)
      else
         call this%relax_M(VF,Q,Pjump,ier)
         if (ier.eq.RELAX_OK) call this%finite_step(2,dt,VF,Q,Pjump,ier)
      end if
      if (present(ierr)) ierr=ier
   end subroutine pT_relax

   !> Instantaneous mechanical and thermal stages, then chemical relaxation at rate nu
   subroutine pTg_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      integer :: ier
      call this%require_unclamped()
      if (this%CRR.ge.huge(1.0_WP).and.this%tauC.le.0.0_WP) then
         ! Instantaneous chemical stage as in the inst class, through the stages of this type
         call this%relax_MT(VF,Q,Pjump,ier)
         if (ier.eq.RELAX_OK) call this%relax_C(VF,Q,Pjump,ier)
      else
         call this%relax_MT(VF,Q,Pjump,ier)
         if (ier.eq.RELAX_OK) call this%finite_step(3,dt,VF,Q,Pjump,ier)
      end if
      if (present(ierr)) ierr=ier
   end subroutine pTg_relax

   !> Finite-rate step of stage ns from a state where the faster forces vanish (eq. finite_solution), then the closure
   subroutine finite_step(this,ns,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      integer,                           intent(in)    :: ns
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,                           intent(out)   :: ierr
      type(phasic_state) :: st
      real(WP), dimension(3,3) :: D,J
      real(WP), dimension(3) :: dx,K
      real(WP) :: rate,omega
      call this%get_state(VF,Q,Pjump,st,ierr)
      if (ierr.ne.RELAX_OK) return
      if (ns.eq.3) then
         ! No transfer at zero chemical rate
         if (this%CRR.le.0.0_WP.and.this%tauC.le.0.0_WP) return
         ! A dry gas is seeded along d_C (tex: Dilute vapor), then F_M=F_T=0 restored
         if (st%dry_gas) then
            D=this%get_dir(st,Pjump)
            call this%transfer(D,[0.0_WP,0.0_WP,this%xv_dry_seed*st%ma*this%gas%R(this%indA)/this%gas%R(this%indV)],Pjump,VF,Q)
            call this%relax_MT(VF,Q,Pjump,ierr)
            if (ierr.ne.RELAX_OK) return
            call this%get_state(VF,Q,Pjump,st,ierr)
            if (ierr.ne.RELAX_OK) return
         end if
      end if
      D=this%get_dir(st,Pjump,ns)
      J=this%get_jac(st,D,Pjump)
      ! Pivot K of the stage at the step start (eq. stiff_def)
      call this%eliminate(ns,J,st%F,1.0_WP,dx,K,ierr)
      if (ierr.ne.RELAX_OK) return
      select case (ns)
      case (1)
         rate=this%MRR; this%KM=K(1)
      case (2)
         rate=this%TRR; this%KT=K(2)
      case default
         ! nu=CRR, or nu=1/(K^C**tauC) so that nu*K^C**dt=dt/tauC
         rate=this%CRR; this%KC=K(3)
         if (this%tauC.gt.0.0_WP) rate=1.0_WP/(K(3)*this%tauC)
         this%CRR_loc=rate; this%tauC_loc=1.0_WP/(rate*K(3))
      end select
      ! omega=1-exp(-r*K*dt), 1 once r*K*dt>=50
      omega=0.0_WP
      if (rate.gt.0.0_WP) then
         omega=1.0_WP
         if (log(rate)+log(K(ns)*dt).lt.log(50.0_WP)) omega=1.0_WP-exp(-rate*K(ns)*dt)
      end if
      call this%eliminate(ns,J,st%F,omega,dx,K,ierr)
      if (ierr.ne.RELAX_OK) return
      call this%transfer(D,dx,Pjump,VF,Q)
      ! Closure of the faster equilibria
      if (ns.eq.1) return
      select case (this%closure)
      case (energy_split)
         call this%mixture_energy_closure(VF,Q,Pjump,ierr)
      case default
         if (ns.eq.2) then
            call this%relax_M(VF,Q,Pjump,ierr)
         else
            call this%relax_MT(VF,Q,Pjump,ierr)
         end if
      end select
   end subroutine finite_step

   !> Pelanti 2022 eqs. 59, 70, 82: V_l and the masses kept, energies split so that p_l-p_g=p_j with E_l+E_g fixed
   !> (tex: Update of pressure and temperature, option 3: heat dE~ along d_T until F^_M=0)
   subroutine mixture_energy_closure(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,                           intent(out)   :: ierr
      type(phasic_state) :: st
      real(WP) :: dE
      integer  :: it
      ! Newton in dE~ for a general equation of state on the root that defines option 3; one step for a stiffened-gas/NASG liquid and an ideal gas
      do it=1,this%itmax
         call this%get_state(VF,Q,Pjump,st,ierr)
         if (ierr.ne.RELAX_OK) return
         if (this%converged(st,1,Pjump)) return
         dE=st%F(1)/(st%Zl(1,2)+st%Zg(1,2))
         Q(3)=Q(3)-dE
         Q(4)=Q(4)+dE
      end do
      ierr=RELAX_FAILED
   end subroutine mixture_energy_closure

end module relax_igmix_max_ent_finite_class
