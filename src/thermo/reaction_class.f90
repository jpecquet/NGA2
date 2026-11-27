!> Abstract gas-phase reaction model: advances the mass fractions of a mixture at fixed density and internal energy.
!> With formation energies in the mixture EOS (igmix q), the heat release follows from the change of composition at fixed
!> internal energy, so a model only changes the mass fractions.
module reaction_class
   use precision, only: WP
   implicit none
   private
   public :: reaction

   type, abstract :: reaction
      integer :: ns=0   !< Number of species the model acts on (the gas mixture's ns)
   contains
      procedure(advance_iface), deferred :: advance
      procedure(rate_iface),    deferred :: get_rate
   end type reaction

   abstract interface
      !> Advance the full mass-fraction vector y over dt at density rho and temperature T (y stays in [0,1] and sums to 1)
      subroutine advance_iface(this,rho,T,y,dt)
         import :: reaction,WP
         class(reaction), intent(in) :: this
         real(WP), intent(in) :: rho,T,dt
         real(WP), dimension(:), intent(inout) :: y
      end subroutine advance_iface
      !> Inverse chemical time scale at (rho,T,y), for the time-step limit (0 if nothing reacts)
      real(WP) function rate_iface(this,rho,T,y)
         import :: reaction,WP
         class(reaction), intent(in) :: this
         real(WP), intent(in) :: rho,T
         real(WP), dimension(:), intent(in) :: y
      end function rate_iface
   end interface

end module reaction_class
