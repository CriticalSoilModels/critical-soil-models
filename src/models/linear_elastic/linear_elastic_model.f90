! Linear elastic model — simplest concrete implementation of csm_model_t.
! No yield surface, no state variables. Always elastic.
! Used to validate the abstract base + integrator stack end-to-end.

module mod_linear_elastic_model
   use stdlib_kinds,      only: dp
   use mod_csm_model,     only: csm_model_t
   use mod_elastic_utils, only: calc_stiffness_GK, calc_K_from_G_nu
   implicit none

   type, extends(csm_model_t) :: linear_elastic_model_t
      real(dp) :: G   !! Shear modulus [kPa]
      real(dp) :: nu  !! Poisson's ratio [-]
   contains
      procedure :: yield_fn          => le_yield_fn
      procedure :: flow_rule         => le_flow_rule
      procedure :: plastic_potential => le_plastic_potential
      procedure :: update_hardening  => le_update_hardening
      procedure :: elastic_stiffness => le_elastic_stiffness
      procedure :: snapshot          => le_snapshot
      procedure :: restore           => le_restore
   end type linear_elastic_model_t

contains

   function le_yield_fn(self, sig) result(F)
      !! Always returns -1: linear elastic material never yields
      class(linear_elastic_model_t), intent(in) :: self
      real(dp), intent(in) :: sig(6)
      real(dp) :: F
      F = -1.0_dp
   end function le_yield_fn

   function le_flow_rule(self, sig) result(dF_by_dsig)
      !! Never called — yield function is always negative
      class(linear_elastic_model_t), intent(in) :: self
      real(dp), intent(in) :: sig(6)
      real(dp) :: dF_by_dsig(6)
      dF_by_dsig = 0.0_dp
   end function le_flow_rule

   function le_plastic_potential(self, sig) result(dg_plas_by_dsig)
      !! Never called — yield function is always negative
      class(linear_elastic_model_t), intent(in) :: self
      real(dp), intent(in) :: sig(6)
      real(dp) :: dg_plas_by_dsig(6)
      dg_plas_by_dsig = 0.0_dp
   end function le_plastic_potential

   subroutine le_update_hardening(self, deps_p)
      !! No-op — no internal state to update
      class(linear_elastic_model_t), intent(inout) :: self
      real(dp), intent(in) :: deps_p(6)
   end subroutine le_update_hardening

   function le_elastic_stiffness(self) result(stiff_e)
      !! Isotropic linear elastic stiffness matrix in Voigt notation
      class(linear_elastic_model_t), intent(in) :: self
      real(dp) :: stiff_e(6,6)

      stiff_e = calc_stiffness_GK(self%G, calc_K_from_G_nu(self%G, self%nu))

   end function le_elastic_stiffness

   subroutine le_snapshot(self, saved)
      !! No state to save — allocates an empty array
      class(linear_elastic_model_t), intent(in)    :: self
      real(dp), allocatable,         intent(out)   :: saved(:)
      allocate(saved(0))
   end subroutine le_snapshot

   subroutine le_restore(self, saved)
      !! No state to restore — no-op
      class(linear_elastic_model_t), intent(inout) :: self
      real(dp),                      intent(in)    :: saved(:)
   end subroutine le_restore

end module mod_linear_elastic_model
