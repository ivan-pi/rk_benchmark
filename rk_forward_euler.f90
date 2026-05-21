module rk_forward_euler
  use rk_kinds
  use rk_types
  use robertson_models, only: rob_direct
  implicit none

contains

  subroutine run_euler_direct(neqn, n_steps, h, t, y)
    integer,  intent(in)    :: neqn, n_steps
    real(dp), intent(in)    :: h
    real(dp), intent(inout) :: t, y(neqn)
    real(dp) :: dy(neqn)
    integer  :: rep

    do rep = 1, n_steps
      call rob_direct(neqn, t, y, dy)
      y = y + h * dy
      t = t + h
    end do
  end subroutine run_euler_direct

  subroutine run_euler_f77(neqn, n_steps, h, t, y, fun)
    integer,  intent(in)    :: neqn, n_steps
    real(dp), intent(in)    :: h
    real(dp), intent(inout) :: t, y(neqn)
    external :: fun
    real(dp) :: dy(neqn)
    integer  :: rep

    do rep = 1, n_steps
      call fun(neqn, t, y, dy)
      y = y + h * dy
      t = t + h
    end do
  end subroutine run_euler_f77

  subroutine run_euler_par(neqn, n_steps, h, t, y, fun, rpar, ipar)
    integer,  intent(in)    :: neqn, n_steps
    real(dp), intent(in)    :: h
    real(dp), intent(inout) :: t, y(neqn)
    procedure(func_par)     :: fun
    real(dp), intent(inout) :: rpar(*)
    integer,  intent(inout) :: ipar(*)
    real(dp) :: dy(neqn)
    integer  :: rep

    do rep = 1, n_steps
      call fun(neqn, t, y, dy, rpar, ipar)
      y = y + h * dy
      t = t + h
    end do
  end subroutine run_euler_par

  subroutine run_euler_cptr(neqn, n_steps, h, t, y, fun, ctx)
    integer,  intent(in)    :: neqn, n_steps
    real(dp), intent(in)    :: h
    real(dp), intent(inout) :: t, y(neqn)
    procedure(func_cptr)    :: fun
    type(c_ptr), value      :: ctx
    real(dp) :: dy(neqn)
    integer  :: rep

    do rep = 1, n_steps
      call fun(neqn, t, y, dy, ctx)
      y = y + h * dy
      t = t + h
    end do
  end subroutine run_euler_cptr

  subroutine run_euler_tb(neqn, n_steps, h, t, y, functor)
    integer,  intent(in)    :: neqn, n_steps
    real(dp), intent(in)    :: h
    real(dp), intent(inout) :: t, y(neqn)
    class(ode_functor), intent(inout) :: functor
    real(dp) :: dy(neqn)
    integer  :: rep

    do rep = 1, n_steps
      call functor%eval(neqn, t, y, dy)
      y = y + h * dy
      t = t + h
    end do
  end subroutine run_euler_tb

  subroutine run_euler_class_star(neqn, n_steps, h, t, y, fun, ctx)
    integer,  intent(in)    :: neqn, n_steps
    real(dp), intent(in)    :: h
    real(dp), intent(inout) :: t, y(neqn)
    procedure(func_class_star) :: fun
    class(*), intent(inout) :: ctx
    real(dp) :: dy(neqn)
    integer  :: rep

    do rep = 1, n_steps
      call fun(neqn, t, y, dy, ctx)
      y = y + h * dy
      t = t + h
    end do
  end subroutine run_euler_class_star

end module rk_forward_euler
