module rk_c_solver
  use rk_kinds, only: c_int, c_ptr, dp
  use rk_types, only: func_cptr, rk_stats
  use iso_c_binding, only: c_double, c_funloc, c_funptr
  implicit none
  private

  public :: rk23_cptr_external

  interface
    subroutine rk23_cptr_external_impl(neqn, fun, t, y, tend, h, atol, rtol, work, ctx, idid, stats) &
        bind(c, name="rk23_cptr_external")
      import :: c_double, c_funptr, c_int, c_ptr, rk_stats
      integer(c_int), intent(in)    :: neqn
      type(c_funptr), value         :: fun
      real(c_double), intent(inout) :: t, y(neqn), h
      real(c_double), intent(in)    :: tend, atol(neqn), rtol
      real(c_double), intent(inout) :: work(neqn, 5)
      type(c_ptr), value            :: ctx
      integer(c_int), intent(out)   :: idid
      type(rk_stats), intent(out)   :: stats
    end subroutine rk23_cptr_external_impl
  end interface

contains

  subroutine rk23_cptr_external(neqn, fun, t, y, tend, h, atol, rtol, work, ctx, idid, stats)
    integer, intent(in) :: neqn
    procedure(func_cptr) :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol(neqn), rtol
    real(dp), intent(inout) :: work(neqn, 5)
    type(c_ptr), value      :: ctx
    integer,  intent(out)   :: idid
    type(rk_stats), intent(out) :: stats

    integer(c_int) :: neqn_c, idid_c

    neqn_c = neqn
    call rk23_cptr_external_impl(neqn_c, c_funloc(fun), t, y, tend, h, atol, rtol, work, ctx, idid_c, stats)
    idid = idid_c
  end subroutine rk23_cptr_external

end module rk_c_solver
