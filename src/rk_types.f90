module rk_types
  use rk_kinds
  implicit none

  integer, parameter :: CTRL_I = 1
  integer, parameter :: CTRL_PI = 2

  ! Statistics collected during a single integration
  type, bind(c) :: rk_stats
    integer(c_int) :: accepted = 0   ! accepted steps
    integer(c_int) :: rejected = 0   ! rejected steps
    integer(c_int) :: nfev     = 0   ! function evaluations
  end type rk_stats

  type :: step_controller
    integer  :: kind = CTRL_I
    real(dp) :: k1 = 0.7_dp
    real(dp) :: k2 = 0.4_dp
    real(dp) :: safety = 0.9_dp
    real(dp) :: fac_min = 0.2_dp
    real(dp) :: fac_max = 5.0_dp
    real(dp) :: err_prev = 1.0_dp
    logical  :: have_prev = .false.
  contains
    procedure :: next_h
    procedure :: accept
    procedure :: reject
  end type step_controller

  abstract interface
    subroutine func_simple(neqn, t, y, ydot)
      import :: dp
      integer,  intent(in)  :: neqn
      real(dp), intent(in)  :: t, y(neqn)
      real(dp), intent(out) :: ydot(neqn)
    end subroutine func_simple

    subroutine func_par(neqn, t, y, ydot, rpar, ipar)
      import :: dp
      integer,  intent(in)  :: neqn
      real(dp), intent(in)  :: t, y(neqn)
      real(dp), intent(out) :: ydot(neqn)
      real(dp), intent(inout) :: rpar(*)
      integer,  intent(inout) :: ipar(*)
    end subroutine func_par

    subroutine func_cptr(neqn, t, y, ydot, ctx) bind(c)
      import :: c_double, c_int, c_ptr
      integer(c_int), value :: neqn
      real(c_double), value :: t
      real(c_double), intent(in)  :: y(neqn)
      real(c_double), intent(out) :: ydot(neqn)
      type(c_ptr), value :: ctx
    end subroutine func_cptr

    subroutine func_class_star(neqn, t, y, ydot, ctx)
      import :: dp
      integer,  intent(in)  :: neqn
      real(dp), intent(in)  :: t, y(neqn)
      real(dp), intent(out) :: ydot(neqn)
      class(*), intent(inout) :: ctx
    end subroutine func_class_star
  end interface

  type, abstract :: ode_functor
  contains
    procedure(func_tb), deferred, pass(this) :: eval
  end type ode_functor

  abstract interface
    subroutine func_tb(this, neqn, t, y, ydot)
      import :: dp, ode_functor
      class(ode_functor), intent(inout) :: this
      integer,  intent(in)  :: neqn
      real(dp), intent(in)  :: t, y(neqn)
      real(dp), intent(out) :: ydot(neqn)
    end subroutine func_tb
  end interface

contains

  pure function next_h(this, err, h, p, step_rejected) result(h_new)
    class(step_controller), intent(in) :: this
    real(dp), intent(in) :: err, h
    integer,  intent(in) :: p
    logical,  intent(in) :: step_rejected
    real(dp) :: h_new
    real(dp) :: fac, q, err_now

    q = real(p + 1, dp)
    err_now = max(err, 1.0e-10_dp)

    select case (this%kind)
    case (CTRL_PI)
      if (this%have_prev) then
        fac = this%safety * err_now**(-this%k1 / q) * max(this%err_prev, 1.0e-10_dp)**(-this%k2 / q)
      else
        fac = this%safety * err_now**(-1.0_dp / q)
      end if
    case default
      fac = this%safety * err_now**(-1.0_dp / q)
    end select

    if (step_rejected) fac = min(1.0_dp, fac)
    h_new = h * max(this%fac_min, min(this%fac_max, fac))
  end function next_h

  pure subroutine accept(this, err)
    class(step_controller), intent(inout) :: this
    real(dp), intent(in) :: err

    this%err_prev = max(err, 1.0e-10_dp)
    this%have_prev = .true.
  end subroutine accept

  pure subroutine reject(this)
    class(step_controller), intent(inout) :: this

    this%have_prev = .false.
  end subroutine reject

end module rk_types
