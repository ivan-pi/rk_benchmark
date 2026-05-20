! ==============================================================================
! Runge-Kutta Benchmark: Strategy Pattern / Callback Overhead in Fortran
! FSAL Property: Direct mapping without temporary buffers
! ==============================================================================
module rk_solvers
  use rk_kinds, only: dp
  use rk_types
  implicit none
  private

  public :: rk23_simple, rk23_par, rk23_cptr, rk23_tb
  public :: rk23_rci, rk23_class_star

contains

  pure function weighted_norm(n, y_old, y_next, err_vec, a_tol, r_tol) result(err_val)
    integer, intent(in)  :: n
    real(dp), intent(in) :: y_old(n), y_next(n), err_vec(n), a_tol, r_tol
    real(dp) :: err_val, scale_vec(n)

    scale_vec = a_tol + max(abs(y_old), abs(y_next)) * r_tol
    err_val = sqrt( sum((err_vec / scale_vec)**2) / n )
  end function weighted_norm

  ! ============================================================================
  ! 1. Internal / Simple Callback Integrator
  ! ============================================================================
  subroutine rk23_simple(neqn, fun, t, y, tend, h, atol, rtol, work, idid)
    integer, intent(in) :: neqn
    procedure(func_simple) :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol, rtol
    real(dp), intent(inout) :: work(neqn, 5)
    integer,  intent(out)   :: idid

    real(dp) :: err, fac, y_new(neqn)

    idid = 0
    ! Evaluate initial k1 directly into column 1
    call fun(neqn, t, y, work(:,1))

    integrate: do while (t < tend)
      if (t + h > tend) h = tend - t

      attempt: do
        if (h <= spacing(t)) then
          idid = -1
          exit integrate
        end if

        call rk23_step(h, y_new, err)
        fac = 0.9_dp * (1.0_dp / max(err, 1.0e-10_dp))**(1.0_dp/3.0_dp)

        if (err <= 1.0_dp) then
          t = t + h
          y = y_new
          work(:,1) = work(:,4) ! FSAL: k4 of accepted step becomes k1 of next
          h = h * max(0.2_dp, min(5.0_dp, fac))
          exit attempt
        end if

        ! Rejected: t and y are unchanged, so work(:,1) is still the valid k1
        h = h * max(0.2_dp, min(5.0_dp, fac))
      end do attempt
    end do integrate

  contains
    subroutine rk23_step(dt, y_next, err_val)
      real(dp), intent(in)  :: dt
      real(dp), intent(out) :: y_next(neqn), err_val
      real(dp) :: err_vec(neqn)

      work(:,5) = y + dt * 0.5_dp * work(:,1)
      call fun(neqn, t + 0.5_dp*dt, work(:,5), work(:,2))

      work(:,5) = y + dt * 0.75_dp * work(:,2)
      call fun(neqn, t + 0.75_dp*dt, work(:,5), work(:,3))

      y_next = y + dt * ((2.0_dp/9.0_dp)*work(:,1) + (1.0_dp/3.0_dp)*work(:,2) + (4.0_dp/9.0_dp)*work(:,3))
      call fun(neqn, t + dt, y_next, work(:,4))

      err_vec = dt * (-5.0_dp/72.0_dp * work(:,1) + 1.0_dp/12.0_dp * work(:,2) + &
                       1.0_dp/9.0_dp * work(:,3) - 1.0_dp/8.0_dp * work(:,4))
      err_val = weighted_norm(neqn, y, y_next, err_vec, atol, rtol)
    end subroutine rk23_step

  end subroutine rk23_simple


  ! ============================================================================
  ! 2. SLATEC-like Callback (rpar, ipar) Integrator
  ! ============================================================================
  subroutine rk23_par(neqn, fun, t, y, tend, h, atol, rtol, work, rpar, ipar, idid)
    integer, intent(in) :: neqn
    procedure(func_par) :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol, rtol
    real(dp), intent(inout) :: work(neqn, 5)
    real(dp), intent(inout) :: rpar(:)
    integer,  intent(inout) :: ipar(:)
    integer,  intent(out)   :: idid

    real(dp) :: err, fac, y_new(neqn)

    idid = 0
    call fun(neqn, t, y, work(:,1), rpar, ipar)

    integrate: do while (t < tend)
      if (t + h > tend) h = tend - t

      attempt: do
        if (h <= spacing(t)) then
          idid = -1
          exit integrate
        end if

        call rk23_step(h, y_new, err)
        fac = 0.9_dp * (1.0_dp / max(err, 1.0e-10_dp))**(1.0_dp/3.0_dp)

        if (err <= 1.0_dp) then
          t = t + h
          y = y_new
          work(:,1) = work(:,4)
          h = h * max(0.2_dp, min(5.0_dp, fac))
          exit attempt
        end if

        h = h * max(0.2_dp, min(5.0_dp, fac))
      end do attempt
    end do integrate

  contains
    subroutine rk23_step(dt, y_next, err_val)
      real(dp), intent(in)  :: dt
      real(dp), intent(out) :: y_next(neqn), err_val
      real(dp) :: err_vec(neqn)

      work(:,5) = y + dt * 0.5_dp * work(:,1)
      call fun(neqn, t + 0.5_dp*dt, work(:,5), work(:,2), rpar, ipar)

      work(:,5) = y + dt * 0.75_dp * work(:,2)
      call fun(neqn, t + 0.75_dp*dt, work(:,5), work(:,3), rpar, ipar)

      y_next = y + dt * ((2.0_dp/9.0_dp)*work(:,1) + (1.0_dp/3.0_dp)*work(:,2) + (4.0_dp/9.0_dp)*work(:,3))
      call fun(neqn, t + dt, y_next, work(:,4), rpar, ipar)

      err_vec = dt * (-5.0_dp/72.0_dp * work(:,1) + 1.0_dp/12.0_dp * work(:,2) + &
                       1.0_dp/9.0_dp * work(:,3) - 1.0_dp/8.0_dp * work(:,4))
      err_val = weighted_norm(neqn, y, y_next, err_vec, atol, rtol)
    end subroutine rk23_step

  end subroutine rk23_par


  ! ============================================================================
  ! 3. C-Pointer Callback Integrator
  ! ============================================================================
  subroutine rk23_cptr(neqn, fun, t, y, tend, h, atol, rtol, work, ctx, idid)
    integer, intent(in) :: neqn
    procedure(func_cptr) :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol, rtol
    real(dp), intent(inout) :: work(neqn, 5)
    type(c_ptr), value      :: ctx
    integer,  intent(out)   :: idid

    real(dp) :: err, fac, y_new(neqn)

    idid = 0
    call fun(neqn, t, y, work(:,1), ctx)

    integrate: do while (t < tend)
      if (t + h > tend) h = tend - t

      attempt: do
        if (h <= spacing(t)) then
          idid = -1
          exit integrate
        end if

        call rk23_step(h, y_new, err)
        fac = 0.9_dp * (1.0_dp / max(err, 1.0e-10_dp))**(1.0_dp/3.0_dp)

        if (err <= 1.0_dp) then
          t = t + h
          y = y_new
          work(:,1) = work(:,4)
          h = h * max(0.2_dp, min(5.0_dp, fac))
          exit attempt
        end if

        h = h * max(0.2_dp, min(5.0_dp, fac))
      end do attempt
    end do integrate

  contains
    subroutine rk23_step(dt, y_next, err_val)
      real(dp), intent(in)  :: dt
      real(dp), intent(out) :: y_next(neqn), err_val
      real(dp) :: err_vec(neqn)

      work(:,5) = y + dt * 0.5_dp * work(:,1)
      call fun(neqn, t + 0.5_dp*dt, work(:,5), work(:,2), ctx)

      work(:,5) = y + dt * 0.75_dp * work(:,2)
      call fun(neqn, t + 0.75_dp*dt, work(:,5), work(:,3), ctx)

      y_next = y + dt * ((2.0_dp/9.0_dp)*work(:,1) + (1.0_dp/3.0_dp)*work(:,2) + (4.0_dp/9.0_dp)*work(:,3))
      call fun(neqn, t + dt, y_next, work(:,4), ctx)

      err_vec = dt * (-5.0_dp/72.0_dp * work(:,1) + 1.0_dp/12.0_dp * work(:,2) + &
                       1.0_dp/9.0_dp * work(:,3) - 1.0_dp/8.0_dp * work(:,4))
      err_val = weighted_norm(neqn, y, y_next, err_vec, atol, rtol)
    end subroutine rk23_step

  end subroutine rk23_cptr


  ! ============================================================================
  ! 4. Functor OOP Class Integrator
  ! ============================================================================
  subroutine rk23_tb(neqn, fun, t, y, tend, h, atol, rtol, work, idid)
    integer, intent(in) :: neqn
    class(ode_functor), intent(inout) :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol, rtol
    real(dp), intent(inout) :: work(neqn, 5)
    integer,  intent(out)   :: idid

    real(dp) :: err, fac, y_new(neqn)

    idid = 0
    call fun%eval(neqn, t, y, work(:,1))

    integrate: do while (t < tend)
      if (t + h > tend) h = tend - t

      attempt: do
        if (h <= spacing(t)) then
          idid = -1
          exit integrate
        end if

        call rk23_step(h, y_new, err)
        fac = 0.9_dp * (1.0_dp / max(err, 1.0e-10_dp))**(1.0_dp/3.0_dp)

        if (err <= 1.0_dp) then
          t = t + h
          y = y_new
          work(:,1) = work(:,4)
          h = h * max(0.2_dp, min(5.0_dp, fac))
          exit attempt
        end if

        h = h * max(0.2_dp, min(5.0_dp, fac))
      end do attempt
    end do integrate

  contains
    subroutine rk23_step(dt, y_next, err_val)
      real(dp), intent(in)  :: dt
      real(dp), intent(out) :: y_next(neqn), err_val
      real(dp) :: err_vec(neqn)

      work(:,5) = y + dt * 0.5_dp * work(:,1)
      call fun%eval(neqn, t + 0.5_dp*dt, work(:,5), work(:,2))

      work(:,5) = y + dt * 0.75_dp * work(:,2)
      call fun%eval(neqn, t + 0.75_dp*dt, work(:,5), work(:,3))

      y_next = y + dt * ((2.0_dp/9.0_dp)*work(:,1) + (1.0_dp/3.0_dp)*work(:,2) + (4.0_dp/9.0_dp)*work(:,3))
      call fun%eval(neqn, t + dt, y_next, work(:,4))

      err_vec = dt * (-5.0_dp/72.0_dp * work(:,1) + 1.0_dp/12.0_dp * work(:,2) + &
                       1.0_dp/9.0_dp * work(:,3) - 1.0_dp/8.0_dp * work(:,4))
      err_val = weighted_norm(neqn, y, y_next, err_vec, atol, rtol)
    end subroutine rk23_step

  end subroutine rk23_tb


  ! ============================================================================
  ! 5. Fully self-contained Reverse Communication Interface (RCI)
  ! Caller passes k1 into work(:,1) before entering the loop.
  ! ============================================================================
  subroutine rk23_rci(stage, neqn, t, y, tend, h, atol, rtol, work, &
                      t_eval, y_eval, idid)
    integer,  intent(inout) :: stage
    integer,  intent(in)    :: neqn
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol, rtol
    real(dp), intent(inout) :: work(neqn, 5)
    real(dp), intent(out)   :: t_eval, y_eval(neqn)
    integer,  intent(out)   :: idid

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
        fac = 0.9_dp * (1.0_dp / max(err, 1.0e-10_dp))**(1.0_dp/3.0_dp)

        if (err <= 1.0_dp) then
          t = t + h
          y = work(:,5)
          work(:,1) = work(:,4) ! FSAL: k4 becomes k1
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
  subroutine rk23_class_star(neqn, fun, t, y, tend, h, atol, rtol, work, ctx, idid)
    integer, intent(in) :: neqn
    procedure(func_class_star) :: fun
    real(dp), intent(inout) :: t, y(neqn), h
    real(dp), intent(in)    :: tend, atol, rtol
    real(dp), intent(inout) :: work(neqn, 5)
    class(*), intent(inout) :: ctx
    integer,  intent(out)   :: idid

    real(dp) :: err, fac, y_new(neqn)

    idid = 0
    call fun(neqn, t, y, work(:,1), ctx)

    integrate: do while (t < tend)
      if (t + h > tend) h = tend - t

      attempt: do
        if (h <= spacing(t)) then
          idid = -1
          exit integrate
        end if

        call rk23_step(h, y_new, err)
        fac = 0.9_dp * (1.0_dp / max(err, 1.0e-10_dp))**(1.0_dp/3.0_dp)

        if (err <= 1.0_dp) then
          t = t + h
          y = y_new
          work(:,1) = work(:,4)
          h = h * max(0.2_dp, min(5.0_dp, fac))
          exit attempt
        end if

        h = h * max(0.2_dp, min(5.0_dp, fac))
      end do attempt
    end do integrate

  contains
    subroutine rk23_step(dt, y_next, err_val)
      real(dp), intent(in)  :: dt
      real(dp), intent(out) :: y_next(neqn), err_val
      real(dp) :: err_vec(neqn)

      work(:,5) = y + dt * 0.5_dp * work(:,1)
      call fun(neqn, t + 0.5_dp*dt, work(:,5), work(:,2), ctx)

      work(:,5) = y + dt * 0.75_dp * work(:,2)
      call fun(neqn, t + 0.75_dp*dt, work(:,5), work(:,3), ctx)

      y_next = y + dt * ((2.0_dp/9.0_dp)*work(:,1) + (1.0_dp/3.0_dp)*work(:,2) + (4.0_dp/9.0_dp)*work(:,3))
      call fun(neqn, t + dt, y_next, work(:,4), ctx)

      err_vec = dt * (-5.0_dp/72.0_dp * work(:,1) + 1.0_dp/12.0_dp * work(:,2) + &
                       1.0_dp/9.0_dp * work(:,3) - 1.0_dp/8.0_dp * work(:,4))
      err_val = weighted_norm(neqn, y, y_next, err_vec, atol, rtol)
    end subroutine rk23_step

  end subroutine rk23_class_star

end module rk_solvers