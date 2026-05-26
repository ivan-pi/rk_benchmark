! Self-contained Robertson ODE plugin compiled as a CMake MODULE library.
! Exported via bind(c) so it can be loaded at runtime with dlopen/dlsym.
! No dependency on rk_kinds or rk_types: uses iso_c_binding directly.
module robertson_plugin_mod
  use iso_c_binding, only: dp => c_double, c_ptr, c_f_pointer
  implicit none
  private

  ! Mirror of robertson_ctx in robertson_models.f90 (bind(c) guarantees layout).
  type, bind(c) :: robertson_ctx_t
    real(dp) :: k1, k2, k3
  end type

contains

  subroutine rob_cptr_plugin(neqn, t, y, ydot, ctx) bind(c, name="rob_cptr")
    integer,  intent(in)  :: neqn
    real(dp), intent(in)  :: t, y(neqn)
    real(dp), intent(out) :: ydot(neqn)
    type(c_ptr), value    :: ctx
    type(robertson_ctx_t), pointer :: d
    call c_f_pointer(ctx, d)
    ydot(1) = -d%k1*y(1) + d%k2*y(2)*y(3)
    ydot(2) =  d%k1*y(1) - d%k2*y(2)*y(3) - d%k3*y(2)**2
    ydot(3) =  d%k3*y(2)**2
  end subroutine

end module robertson_plugin_mod
