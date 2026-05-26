! ==============================================================================
! 1-D Reaction-Diffusion Models for the RK Benchmark Suite
!
! Problem: du/dt = D * d^2u/dx^2 - lambda * u   on x in (0,1), t > 0
!          u(0,t) = u(1,t) = 0  (Dirichlet BCs)
!          u(x,0) = sin(pi*x)   (initial condition)
!
! Spatial discretisation: N interior points x_i = i/(N+1), i = 1..N
! Central-difference Laplacian with zero ghost values at i=0, i=N+1.
!
! The initial condition is the lowest FD eigenfunction.  Its exact
! time evolution gives a pure exponential decay:
!
!   u_i(t) = sin(pi*i/(N+1)) * exp(-(mu1 + lambda)*t)
!
! where  mu1 = D*(N+1)^2 * (2 - 2*cos(pi/(N+1)))  is the lowest
! eigenvalue of the discrete Laplacian operator.
!
! Because the exact solution satisfies the discrete ODE exactly (for
! the single-mode initial condition), error against rd_exact reflects
! only the time-integration error, not spatial truncation error.
! ==============================================================================
module reaction_diffusion_models
  use rk_kinds
  use rk_types
  use iso_fortran_env, only: error_unit
  implicit none

  ! Fixed model parameters used by the module-level callbacks
  ! (rd_direct, for F77-style and RCI strategies).
  real(dp), parameter :: rd_D      = 0.1_dp
  real(dp), parameter :: rd_lambda = 1.0_dp

  ! C-interoperable context for the C-pointer callback strategy
  type, bind(c) :: rd_ctx
    real(c_double) :: D
    real(c_double) :: lambda
  end type

  ! Context type for the class(*) polymorphic callback strategy
  type :: rd_params
    real(dp) :: D      = rd_D
    real(dp) :: lambda = rd_lambda
  end type

  ! OOP functor type for the type-bound-procedure strategy
  type, extends(ode_functor) :: rd_functor
    real(dp) :: D      = rd_D
    real(dp) :: lambda = rd_lambda
  contains
    procedure, pass(this) :: eval => rd_eval
  end type

contains

  ! ============================================================================
  ! Core tridiagonal diffusion-reaction kernel (autonomous, pure).
  ! Computes ydot_i = D*(y_{i-1} - 2*y_i + y_{i+1})/dx^2 - lambda*y_i
  ! with y_0 = y_{N+1} = 0 (Dirichlet BCs), dx = 1/(neqn+1).
  ! ============================================================================
  pure subroutine rd_rhs_core(neqn, D, lambda, y, ydot)
    integer,  intent(in)  :: neqn
    real(dp), intent(in)  :: D, lambda
    real(dp), intent(in)  :: y(neqn)
    real(dp), intent(out) :: ydot(neqn)

    real(dp) :: inv_dx2
    integer  :: i

    ! inv_dx2 = D / dx^2 = D * (neqn+1)^2
    inv_dx2 = D * real(neqn + 1, dp)**2

    if (neqn == 1) then
      ydot(1) = -2.0_dp * inv_dx2 * y(1) - lambda * y(1)
      return
    end if

    ! Left boundary (i=1): left ghost y_0 = 0
    ydot(1) = inv_dx2 * (-2.0_dp*y(1) + y(2)) - lambda*y(1)

    ! Interior points
    do i = 2, neqn - 1
      ydot(i) = inv_dx2 * (y(i-1) - 2.0_dp*y(i) + y(i+1)) - lambda*y(i)
    end do

    ! Right boundary (i=neqn): right ghost y_{neqn+1} = 0
    ydot(neqn) = inv_dx2 * (y(neqn-1) - 2.0_dp*y(neqn)) - lambda*y(neqn)
  end subroutine rd_rhs_core


  ! ============================================================================
  ! 1. Module-parameter callback (used for F77-style and RCI strategies).
  !    Reads D and lambda from the module-level parameters rd_D, rd_lambda.
  ! ============================================================================
  subroutine rd_direct(neqn, t, y, ydot)
    integer,  intent(in)  :: neqn
    real(dp), intent(in)  :: t, y(neqn)
    real(dp), intent(out) :: ydot(neqn)
    call rd_rhs_core(neqn, rd_D, rd_lambda, y, ydot)
  end subroutine rd_direct


  ! ============================================================================
  ! 2. RPAR/IPAR callback: rpar(1) = D, rpar(2) = lambda
  ! ============================================================================
  subroutine rd_par(neqn, t, y, ydot, rpar, ipar)
    integer,  intent(in)    :: neqn
    real(dp), intent(in)    :: t, y(neqn)
    real(dp), intent(out)   :: ydot(neqn)
    real(dp), intent(inout) :: rpar(*)
    integer,  intent(inout) :: ipar(*)
    call rd_rhs_core(neqn, rpar(1), rpar(2), y, ydot)
  end subroutine rd_par


  ! ============================================================================
  ! 3. C-pointer callback: context carries D and lambda
  ! ============================================================================
  subroutine rd_cptr(neqn, t, y, ydot, ctx) bind(c)
    integer(c_int), value       :: neqn
    real(c_double), value       :: t
    real(c_double), intent(in)  :: y(neqn)
    real(c_double), intent(out) :: ydot(neqn)
    type(c_ptr), value          :: ctx
    type(rd_ctx), pointer       :: d
    call c_f_pointer(ctx, d)
    call rd_rhs_core(int(neqn), d%D, d%lambda, y, ydot)
  end subroutine rd_cptr


  ! ============================================================================
  ! 4. OOP functor type-bound procedure
  ! ============================================================================
  subroutine rd_eval(this, neqn, t, y, ydot)
    class(rd_functor), intent(inout) :: this
    integer,  intent(in)  :: neqn
    real(dp), intent(in)  :: t, y(neqn)
    real(dp), intent(out) :: ydot(neqn)
    call rd_rhs_core(neqn, this%D, this%lambda, y, ydot)
  end subroutine rd_eval


  ! ============================================================================
  ! 5. Unlimited-polymorphic class(*) callback
  ! ============================================================================
  subroutine rd_class_star(neqn, t, y, ydot, ctx)
    integer,  intent(in)    :: neqn
    real(dp), intent(in)    :: t, y(neqn)
    real(dp), intent(out)   :: ydot(neqn)
    class(*), intent(inout) :: ctx
    select type (d => ctx)
    class is (rd_params)
      call rd_rhs_core(neqn, d%D, d%lambda, y, ydot)
    class default
      write(error_unit,'(A)') &
        "reaction_diffusion_models: unexpected ctx type in rd_class_star"
      error stop 99
    end select
  end subroutine rd_class_star


  ! ============================================================================
  ! Exact time-evolution of the discrete ODE system for the lowest eigenmode.
  !
  ! y_i(t) = sin(pi*i/(neqn+1)) * exp(-(mu1 + lambda)*t)
  !
  ! where mu1 = D*(neqn+1)^2 * (2 - 2*cos(pi/(neqn+1))) is the smallest
  ! eigenvalue (in magnitude) of the discrete diffusion operator.
  ! This satisfies the ODE system exactly (no spatial truncation error).
  ! ============================================================================
  pure subroutine rd_exact(neqn, t, D, lambda, y_ex)
    integer,  intent(in)  :: neqn
    real(dp), intent(in)  :: t, D, lambda
    real(dp), intent(out) :: y_ex(neqn)

    real(dp) :: pi, mu1, decay
    integer  :: i

    pi    = acos(-1.0_dp)
    mu1   = D * real(neqn+1, dp)**2 * (2.0_dp - 2.0_dp*cos(pi/real(neqn+1, dp)))
    decay = exp(-(mu1 + lambda)*t)
    do i = 1, neqn
      y_ex(i) = sin(pi * real(i, dp) / real(neqn+1, dp)) * decay
    end do
  end subroutine rd_exact

end module reaction_diffusion_models
