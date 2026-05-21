module robertson_models
  use rk_kinds
  use rk_types
  use iso_fortran_env, only: error_unit
  implicit none

  ! Robertson coefficients
  real(dp), parameter :: k1 = 0.04_dp, k2 = 1.0e4_dp, k3 = 3.0e7_dp

  type, bind(c) :: robertson_ctx
    real(c_double) :: k1, k2, k3
  end type

  type :: robertson_params
    real(c_double) :: k1 = 0.04_dp, k2 = 1.0e4_dp, k3 = 3.0e7_dp
  end type

  type, extends(ode_functor) :: robertson_functor
    real(dp) :: k1 = 0.04_dp, k2 = 1.0e4_dp, k3 = 3.0e7_dp
  contains
    procedure, pass(this) :: eval => rob_eval
  end type

contains

  subroutine rob_par(neqn, t, y, ydot, rpar, ipar)
    integer, intent(in) :: neqn
    real(dp), intent(in) :: t, y(neqn)
    real(dp), intent(out) :: ydot(neqn)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    ydot(1) = -rpar(1)*y(1) + rpar(2)*y(2)*y(3)
    ydot(2) =  rpar(1)*y(1) - rpar(2)*y(2)*y(3) - rpar(3)*y(2)**2
    ydot(3) =  rpar(3)*y(2)**2
  end subroutine

  subroutine rob_cptr(neqn, t, y, ydot, ctx) bind(c)
    integer(c_int), value :: neqn
    real(c_double), value :: t
    real(c_double), intent(in) :: y(neqn)
    real(c_double), intent(out) :: ydot(neqn)
    type(c_ptr), value :: ctx
    type(robertson_ctx), pointer :: d
    call c_f_pointer(ctx, d)
    ydot(1) = -d%k1*y(1) + d%k2*y(2)*y(3)
    ydot(2) =  d%k1*y(1) - d%k2*y(2)*y(3) - d%k3*y(2)**2
    ydot(3) =  d%k3*y(2)**2
  end subroutine

  subroutine rob_eval(this, neqn, t, y, ydot)
    class(robertson_functor), intent(inout) :: this
    integer, intent(in) :: neqn
    real(dp), intent(in) :: t, y(neqn)
    real(dp), intent(out) :: ydot(neqn)
    ydot(1) = -this%k1*y(1) + this%k2*y(2)*y(3)
    ydot(2) =  this%k1*y(1) - this%k2*y(2)*y(3) - this%k3*y(2)**2
    ydot(3) =  this%k3*y(2)**2
  end subroutine

  subroutine rob_direct(neqn, t, y, ydot)
    integer, intent(in) :: neqn
    real(dp), intent(in) :: t, y(neqn)
    real(dp), intent(out) :: ydot(neqn)
    ydot(1) = -k1*y(1) + k2*y(2)*y(3)
    ydot(2) =  k1*y(1) - k2*y(2)*y(3) - k3*y(2)**2
    ydot(3) =  k3*y(2)**2
  end subroutine

  subroutine rob_class_star(neqn, t, y, ydot, ctx)
    integer,  intent(in)  :: neqn
    real(dp), intent(in)  :: t, y(neqn)
    real(dp), intent(out) :: ydot(neqn)
    class(*), intent(inout) :: ctx

    ! Dynamic type resolution
    select type (d => ctx)
    class is (robertson_params)
      ydot(1) = -d%k1*y(1) + d%k2*y(2)*y(3)
      ydot(2) =  d%k1*y(1) - d%k2*y(2)*y(3) - d%k3*y(2)**2
      ydot(3) =  d%k3*y(2)**2
    class default
      write(error_unit,'(A)') "robertson_models: wrong type of ctx argument in rob_class_star"
      error stop 99
    end select
  end subroutine rob_class_star

end module robertson_models
