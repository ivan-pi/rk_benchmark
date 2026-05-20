module rk_types
  use rk_kinds
  implicit none

  ! Statistics collected during a single integration
  type, bind(c) :: rk_stats
    integer(c_int) :: accepted = 0   ! accepted steps
    integer(c_int) :: rejected = 0   ! rejected steps
    integer(c_int) :: nfev     = 0   ! function evaluations
  end type rk_stats

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
      real(c_double), intent(in)  :: y(*)
      real(c_double), intent(out) :: ydot(*)
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

end module rk_types
