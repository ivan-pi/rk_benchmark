! ==============================================================================
! Runge-Kutta Benchmark: Strategy Pattern / Callback Overhead in Fortran
! FSAL Property: Direct mapping without temporary buffers
! ==============================================================================
module rk_solvers
  use rk_kinds, only: dp
  use rk_types
  use iso_c_binding, only: c_double, c_f_procpointer, c_funptr, c_int, c_ptr
  implicit none
  private

  public :: rk23_simple, rk23_par, rk23_cptr, rk23_cptr_external_impl, rk23_tb
  public :: rk23_rci, rk23_class_star
  public :: rk_stats

  interface
    subroutine rk23_cptr_external_impl(neqn, fun, t, y, tend, h, atol, rtol, work, ctx, idid, stats) &
        bind(c, name="rk23_cptr_external")
      import :: c_double, c_funptr, c_int, c_ptr, rk_stats
      integer(c_int), value        :: neqn
      type(c_funptr), value        :: fun
      real(c_double), intent(inout) :: t, h
      real(c_double), intent(inout) :: y(*)
      real(c_double), value        :: tend, rtol
      real(c_double), intent(in)   :: atol(*)
      real(c_double), intent(inout) :: work(*)
      type(c_ptr), value           :: ctx
      integer(c_int), intent(out)  :: idid
      type(rk_stats), intent(out)  :: stats
    end subroutine rk23_cptr_external_impl
  end interface


  interface
    function cbrt(x) bind(c,name="cbrt")
      use, intrinsic :: iso_c_binding, only: c_double
      real(c_double), value :: x
      real(c_double) :: cbrt
    end function
  end interface

contains

  pure function weighted_norm(n, y_old, y_next, err_vec, a_tol, r_tol) result(err_val)
    integer, intent(in)  :: n
    real(dp), intent(in) :: y_old(n), y_next(n), err_vec(n), a_tol(n), r_tol
    real(dp) :: err_val, scale_i, sum_sq
    integer :: i

    sum_sq = 0.0_dp
    do i = 1, n
      scale_i = a_tol(i) + max(abs(y_old(i)), abs(y_next(i))) * r_tol
      sum_sq = sum_sq + (err_vec(i) / scale_i)**2
    end do
    err_val = sqrt(sum_sq / n)
  end function weighted_norm

  ! ============================================================================
  ! 1. F77-Style External Callback Integrator (implicit interface)
  !    The callback API uses only Fortran 77 features: no intent declarations,
  !    explicit-size arrays, and an external dummy argument with implicit
  !    interface.  Per Fortran 2008, an internal procedure may still be passed
  !    as the actual argument.
  ! ============================================================================
  subroutine rk23_simple(neqn, fun, t, y, tend, h, atol, rtol, work, idid, stats)
    integer, intent(in) :: neqn
    external :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol(neqn), rtol
    real(dp), intent(inout) :: work(neqn, 5)
    integer,  intent(out)   :: idid
    type(rk_stats), intent(out) :: stats

    real(dp) :: err, fac, y_new(neqn)
    logical  :: step_rejected

    idid = 0
    stats = rk_stats()
    ! Evaluate initial k1 directly into column 1
    call fun(neqn, t, y, work(:,1))
    stats%nfev = stats%nfev + 1

    integrate: do while (t < tend)
      if (t + h > tend) h = tend - t
      step_rejected = .false.

      attempt: do
        if (h <= spacing(t)) then
          idid = -1
          exit integrate
        end if

        call rk23_step(neqn, fun, t, h, y, work(:,1), work(:,2), work(:,3), work(:,4), work(:,5), atol, rtol, y_new, err)
        stats%nfev = stats%nfev + 3
        fac = 0.9_dp * cbrt(1.0_dp / max(err, 1.0e-10_dp))

        if (err < 1.0_dp) then
          if (step_rejected) fac = min(1.0_dp, fac)
          t = t + h
          y = y_new
          work(:,1) = work(:,4) ! FSAL: k4 of accepted step becomes k1 of next
          h = h * max(0.2_dp, min(5.0_dp, fac))
          stats%accepted = stats%accepted + 1
          exit attempt
        end if

        ! Rejected: t and y are unchanged, so work(:,1) is still the valid k1
        h = h * max(0.2_dp, min(5.0_dp, fac))
        step_rejected = .true.
        stats%rejected = stats%rejected + 1
      end do attempt
    end do integrate

  contains
    subroutine rk23_step(n, fn, t_cur, dt, y_cur, k1, k2, k3, k4, tmp, a_tol, r_tol, y_next, err_val)
      integer,  intent(in)  :: n
      external :: fn
      real(dp), intent(in)  :: t_cur, dt, y_cur(n), k1(n), a_tol(n), r_tol
      real(dp), intent(out) :: k2(n), k3(n), k4(n), tmp(n), y_next(n), err_val

      tmp = y_cur + dt * 0.5_dp * k1
      call fn(n, t_cur + 0.5_dp*dt, tmp, k2)

      tmp = y_cur + dt * 0.75_dp * k2
      call fn(n, t_cur + 0.75_dp*dt, tmp, k3)

      y_next = y_cur + dt * ((2.0_dp/9.0_dp)*k1 + (1.0_dp/3.0_dp)*k2 + (4.0_dp/9.0_dp)*k3)
      call fn(n, t_cur + dt, y_next, k4)

      tmp = dt * (-5.0_dp/72.0_dp * k1 + 1.0_dp/12.0_dp * k2 + &
                   1.0_dp/9.0_dp * k3 - 1.0_dp/8.0_dp * k4)
      err_val = weighted_norm(n, y_cur, y_next, tmp, a_tol, r_tol)
    end subroutine rk23_step

  end subroutine rk23_simple


  ! ============================================================================
  ! 2. SLATEC-like Callback (rpar, ipar) Integrator
  ! ============================================================================
  subroutine rk23_par(neqn, fun, t, y, tend, h, atol, rtol, work, rpar, ipar, idid, stats)
    integer, intent(in) :: neqn
    procedure(func_par) :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol(neqn), rtol
    real(dp), intent(inout) :: work(neqn, 5)
    real(dp), intent(inout) :: rpar(*)
    integer,  intent(inout) :: ipar(*)
    integer,  intent(out)   :: idid
    type(rk_stats), intent(out) :: stats

    real(dp) :: err, fac, y_new(neqn)
    logical  :: step_rejected

    idid = 0
    stats = rk_stats()
    call fun(neqn, t, y, work(:,1), rpar, ipar)
    stats%nfev = stats%nfev + 1

    integrate: do while (t < tend)
      if (t + h > tend) h = tend - t
      step_rejected = .false.

      attempt: do
        if (h <= spacing(t)) then
          idid = -1
          exit integrate
        end if

        call rk23_step(neqn, fun, t, h, y, work(:,1), work(:,2), work(:,3), &
                       work(:,4), work(:,5), atol, rtol, rpar, ipar, y_new, err)
        stats%nfev = stats%nfev + 3
        fac = 0.9_dp * cbrt(1.0_dp / max(err, 1.0e-10_dp))

        if (err < 1.0_dp) then
          if (step_rejected) fac = min(1.0_dp, fac)
          t = t + h
          y = y_new
          work(:,1) = work(:,4)
          h = h * max(0.2_dp, min(5.0_dp, fac))
          stats%accepted = stats%accepted + 1
          exit attempt
        end if

        h = h * max(0.2_dp, min(5.0_dp, fac))
        step_rejected = .true.
        stats%rejected = stats%rejected + 1
      end do attempt
    end do integrate

  contains
    subroutine rk23_step(n, fn, t_cur, dt, y_cur, k1, k2, k3, k4, tmp, a_tol, r_tol, rp, ip, y_next, err_val)
      integer,  intent(in)  :: n
      procedure(func_par) :: fn
      real(dp), intent(in)  :: t_cur, dt, y_cur(n), k1(n), a_tol(n), r_tol
      real(dp), intent(inout) :: rp(*)
      integer,  intent(inout) :: ip(*)
      real(dp), intent(out) :: k2(n), k3(n), k4(n), tmp(n), y_next(n), err_val

      tmp = y_cur + dt * 0.5_dp * k1
      call fn(n, t_cur + 0.5_dp*dt, tmp, k2, rp, ip)

      tmp = y_cur + dt * 0.75_dp * k2
      call fn(n, t_cur + 0.75_dp*dt, tmp, k3, rp, ip)

      y_next = y_cur + dt * ((2.0_dp/9.0_dp)*k1 + (1.0_dp/3.0_dp)*k2 + (4.0_dp/9.0_dp)*k3)
      call fn(n, t_cur + dt, y_next, k4, rp, ip)

      tmp = dt * (-5.0_dp/72.0_dp * k1 + 1.0_dp/12.0_dp * k2 + &
                   1.0_dp/9.0_dp * k3 - 1.0_dp/8.0_dp * k4)
      err_val = weighted_norm(n, y_cur, y_next, tmp, a_tol, r_tol)
    end subroutine rk23_step

  end subroutine rk23_par


  ! ============================================================================
  ! 3. C-Pointer Callback Integrator
  ! ============================================================================
  subroutine rk23_cptr(neqn, fun, t, y, tend, h, atol, rtol, work, ctx, idid, stats)
    integer(c_int), value :: neqn
    type(c_funptr), value :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol(neqn), rtol
    real(dp), intent(inout) :: work(neqn, 5)
    type(c_ptr), value      :: ctx
    integer(c_int), intent(out) :: idid
    type(rk_stats), intent(out) :: stats
 
    integer :: neqn_f
    procedure(func_cptr), pointer :: fun_proc
    real(dp) :: err, fac, y_new(neqn)
    logical  :: step_rejected
 
    neqn_f = int(neqn, kind(neqn_f))
    idid = 0
    stats = rk_stats()
    call c_f_procpointer(fun, fun_proc)
    call fun_proc(neqn, t, y, work(:,1), ctx)
    stats%nfev = stats%nfev + 1

    integrate: do while (t < tend)
      if (t + h > tend) h = tend - t
      step_rejected = .false.

      attempt: do
        if (h <= spacing(t)) then
          idid = -1_c_int
          exit integrate
        end if
 
        call rk23_step(neqn_f, fun_proc, t, h, y, work(:,1), work(:,2), work(:,3), &
                       work(:,4), work(:,5), atol, rtol, ctx, y_new, err)
        stats%nfev = stats%nfev + 3
        fac = 0.9_dp * cbrt(1.0_dp / max(err, 1.0e-10_dp))

        if (err < 1.0_dp) then
          if (step_rejected) fac = min(1.0_dp, fac)
          t = t + h
          y = y_new
          work(:,1) = work(:,4)
          h = h * max(0.2_dp, min(5.0_dp, fac))
          stats%accepted = stats%accepted + 1
          exit attempt
        end if

        h = h * max(0.2_dp, min(5.0_dp, fac))
        step_rejected = .true.
        stats%rejected = stats%rejected + 1
      end do attempt
    end do integrate

  contains
    subroutine rk23_step(n, fn, t_cur, dt, y_cur, k1, k2, k3, k4, tmp, a_tol, r_tol, cctx, y_next, err_val)
      integer,  intent(in)  :: n
      procedure(func_cptr) :: fn
      real(dp), intent(in)  :: t_cur, dt, y_cur(n), k1(n), a_tol(n), r_tol
      type(c_ptr), value    :: cctx
      real(dp), intent(out) :: k2(n), k3(n), k4(n), tmp(n), y_next(n), err_val

      integer(c_int) :: n_c

      n_c = n
      tmp = y_cur + dt * 0.5_dp * k1
      call fn(n_c, t_cur + 0.5_dp*dt, tmp, k2, cctx)

      tmp = y_cur + dt * 0.75_dp * k2
      call fn(n_c, t_cur + 0.75_dp*dt, tmp, k3, cctx)

      y_next = y_cur + dt * ((2.0_dp/9.0_dp)*k1 + (1.0_dp/3.0_dp)*k2 + (4.0_dp/9.0_dp)*k3)
      call fn(n_c, t_cur + dt, y_next, k4, cctx)

      tmp = dt * (-5.0_dp/72.0_dp * k1 + 1.0_dp/12.0_dp * k2 + &
                   1.0_dp/9.0_dp * k3 - 1.0_dp/8.0_dp * k4)
      err_val = weighted_norm(n, y_cur, y_next, tmp, a_tol, r_tol)
    end subroutine rk23_step

  end subroutine rk23_cptr


  ! ============================================================================
  ! 4. Functor OOP Class Integrator
  ! ============================================================================
  subroutine rk23_tb(neqn, fun, t, y, tend, h, atol, rtol, work, idid, stats)
    integer, intent(in) :: neqn
    class(ode_functor), intent(inout) :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol(neqn), rtol
    real(dp), intent(inout) :: work(neqn, 5)
    integer,  intent(out)   :: idid
    type(rk_stats), intent(out) :: stats

    real(dp) :: err, fac, y_new(neqn)
    logical  :: step_rejected

    idid = 0
    stats = rk_stats()
    call fun%eval(neqn, t, y, work(:,1))
    stats%nfev = stats%nfev + 1

    integrate: do while (t < tend)
      if (t + h > tend) h = tend - t
      step_rejected = .false.

      attempt: do
        if (h <= spacing(t)) then
          idid = -1
          exit integrate
        end if

        call rk23_step(neqn, fun, t, h, y, work(:,1), work(:,2), work(:,3), work(:,4), work(:,5), atol, rtol, y_new, err)
        stats%nfev = stats%nfev + 3
        fac = 0.9_dp * cbrt(1.0_dp / max(err, 1.0e-10_dp))

        if (err < 1.0_dp) then
          if (step_rejected) fac = min(1.0_dp, fac)
          t = t + h
          y = y_new
          work(:,1) = work(:,4)
          h = h * max(0.2_dp, min(5.0_dp, fac))
          stats%accepted = stats%accepted + 1
          exit attempt
        end if

        h = h * max(0.2_dp, min(5.0_dp, fac))
        step_rejected = .true.
        stats%rejected = stats%rejected + 1
      end do attempt
    end do integrate

  contains
    subroutine rk23_step(n, fn, t_cur, dt, y_cur, k1, k2, k3, k4, tmp, a_tol, r_tol, y_next, err_val)
      integer,  intent(in)  :: n
      class(ode_functor), intent(inout) :: fn
      real(dp), intent(in)  :: t_cur, dt, y_cur(n), k1(n), a_tol(n), r_tol
      real(dp), intent(out) :: k2(n), k3(n), k4(n), tmp(n), y_next(n), err_val

      tmp = y_cur + dt * 0.5_dp * k1
      call fn%eval(n, t_cur + 0.5_dp*dt, tmp, k2)

      tmp = y_cur + dt * 0.75_dp * k2
      call fn%eval(n, t_cur + 0.75_dp*dt, tmp, k3)

      y_next = y_cur + dt * ((2.0_dp/9.0_dp)*k1 + (1.0_dp/3.0_dp)*k2 + (4.0_dp/9.0_dp)*k3)
      call fn%eval(n, t_cur + dt, y_next, k4)

      tmp = dt * (-5.0_dp/72.0_dp * k1 + 1.0_dp/12.0_dp * k2 + &
                   1.0_dp/9.0_dp * k3 - 1.0_dp/8.0_dp * k4)
      err_val = weighted_norm(n, y_cur, y_next, tmp, a_tol, r_tol)
    end subroutine rk23_step

  end subroutine rk23_tb


  ! ============================================================================
  ! 5. Fully self-contained Reverse Communication Interface (RCI)
  ! Caller passes k1 into work(:,1) before entering the loop.
  ! ============================================================================
  subroutine rk23_rci(stage, neqn, t, y, tend, h, atol, rtol, work, &
                      t_eval, y_eval, idid, stats, step_rejected)
    integer,  intent(inout) :: stage
    integer,  intent(in)    :: neqn
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol(neqn), rtol
    real(dp), intent(inout) :: work(neqn, 5)
    real(dp), intent(out)   :: t_eval, y_eval(neqn)
    integer,  intent(out)   :: idid
    type(rk_stats), intent(inout) :: stats
    logical,  intent(inout) :: step_rejected

    real(dp) :: err, fac, err_vec(neqn)

    drive: do
      select case(stage)
      case(1) ! Step begins: work(:,1) is already valid. Prepare inputs for k2.
        if (t >= tend) then; stage = 6; return; end if
        if (t + h > tend) h = tend - t
        if (h <= spacing(t)) then; idid = -1; stage = 6; return; end if

        t_eval = t + 0.5_dp * h; y_eval = y + h * 0.5_dp * work(:,1)
        stage = 2; return ! Yield to evaluate f2 into work(:,2)

      case(2) ! After f2, prepare inputs for k3
        t_eval = t + 0.75_dp * h; y_eval = y + h * 0.75_dp * work(:,2)
        stage = 3; return ! Yield to evaluate f3 into work(:,3)

      case(3) ! After f3, calculate y_new into work(:,5) and ask for k4
        work(:,5) = y + h * ((2.0_dp/9.0_dp)*work(:,1) + (1.0_dp/3.0_dp)*work(:,2) + (4.0_dp/9.0_dp)*work(:,3))
        t_eval = t + h; y_eval = work(:,5)
        stage = 4; return ! Yield to evaluate f4 into work(:,4)

      case(4) ! After f4: calculate error and finalize step
        err_vec = h * (-5.0_dp/72.0_dp * work(:,1) + 1.0_dp/12.0_dp * work(:,2) + &
                        1.0_dp/9.0_dp * work(:,3) - 1.0_dp/8.0_dp * work(:,4))

        err = weighted_norm(neqn, y, work(:,5), err_vec, atol, rtol)
        fac = 0.9_dp * cbrt(1.0_dp / max(err, 1.0e-10_dp))

        stats%nfev = stats%nfev + 3
        if (err < 1.0_dp) then
          if (step_rejected) fac = min(1.0_dp, fac)
          t = t + h
          y = work(:,5)
          work(:,1) = work(:,4) ! FSAL: k4 becomes k1
          step_rejected = .false. ! reset for the next step
          stats%accepted = stats%accepted + 1
        else
          step_rejected = .true.
          stats%rejected = stats%rejected + 1
        end if

        h = h * max(0.2_dp, min(5.0_dp, fac))
        stage = 1 ! Cycle internally back to start of step (or retry)

      case default
        idid = -2; stage = 6; return
      end select
    end do drive
  end subroutine rk23_rci


  ! ============================================================================
  ! 6. Unlimited Polymorphic Class(*) Callback Integrator
  ! ============================================================================
  subroutine rk23_class_star(neqn, fun, t, y, tend, h, atol, rtol, work, ctx, idid, stats)
    integer, intent(in) :: neqn
    procedure(func_class_star) :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol(neqn), rtol
    real(dp), intent(inout) :: work(neqn, 5)
    class(*), intent(inout) :: ctx
    integer,  intent(out)   :: idid
    type(rk_stats), intent(out) :: stats

    real(dp) :: err, fac, y_new(neqn)
    logical  :: step_rejected

    idid = 0
    stats = rk_stats()
    call fun(neqn, t, y, work(:,1), ctx)
    stats%nfev = stats%nfev + 1

    integrate: do while (t < tend)
      if (t + h > tend) h = tend - t
      step_rejected = .false.

      attempt: do
        if (h <= spacing(t)) then
          idid = -1
          exit integrate
        end if

        call rk23_step(neqn, fun, t, h, y, work(:,1), work(:,2), work(:,3), work(:,4), work(:,5), atol, rtol, ctx, y_new, err)
        stats%nfev = stats%nfev + 3
        fac = 0.9_dp * cbrt(1.0_dp / max(err, 1.0e-10_dp))

        if (err < 1.0_dp) then
          if (step_rejected) fac = min(1.0_dp, fac)
          t = t + h
          y = y_new
          work(:,1) = work(:,4)
          h = h * max(0.2_dp, min(5.0_dp, fac))
          stats%accepted = stats%accepted + 1
          exit attempt
        end if

        h = h * max(0.2_dp, min(5.0_dp, fac))
        step_rejected = .true.
        stats%rejected = stats%rejected + 1
      end do attempt
    end do integrate

  contains
    subroutine rk23_step(n, fn, t_cur, dt, y_cur, k1, k2, k3, k4, tmp, a_tol, r_tol, cctx, y_next, err_val)
      integer,  intent(in)  :: n
      procedure(func_class_star) :: fn
      real(dp), intent(in)  :: t_cur, dt, y_cur(n), k1(n), a_tol(n), r_tol
      class(*), intent(inout) :: cctx
      real(dp), intent(out) :: k2(n), k3(n), k4(n), tmp(n), y_next(n), err_val

      tmp = y_cur + dt * 0.5_dp * k1
      call fn(n, t_cur + 0.5_dp*dt, tmp, k2, cctx)

      tmp = y_cur + dt * 0.75_dp * k2
      call fn(n, t_cur + 0.75_dp*dt, tmp, k3, cctx)

      y_next = y_cur + dt * ((2.0_dp/9.0_dp)*k1 + (1.0_dp/3.0_dp)*k2 + (4.0_dp/9.0_dp)*k3)
      call fn(n, t_cur + dt, y_next, k4, cctx)

      tmp = dt * (-5.0_dp/72.0_dp * k1 + 1.0_dp/12.0_dp * k2 + &
                   1.0_dp/9.0_dp * k3 - 1.0_dp/8.0_dp * k4)
      err_val = weighted_norm(n, y_cur, y_next, tmp, a_tol, r_tol)
    end subroutine rk23_step

  end subroutine rk23_class_star

end module rk_solvers
