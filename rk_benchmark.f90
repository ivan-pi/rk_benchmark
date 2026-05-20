program rk_benchmark
  use rk_kinds
  use rk_types
  use rk_solvers
  use robertson_models
  use iso_fortran_env, only: error_unit
  implicit none

  integer,  parameter :: neqn = 3, N_runs = 100
  real(dp), parameter :: y_init(neqn) = [1.0_dp, 0.0_dp, 0.0_dp]
  real(dp), parameter :: t_start = 0.0_dp, t_end = 100.0_dp
  real(dp), parameter :: atol(neqn) = 1.0e-8_dp, rtol = 1.0e-8_dp
  character(len=30), parameter :: plot_labels(6) = [ &
    "F77 Ext. (implicit iface)      ", &
    "Callback with RPAR/IPAR        ", &
    "Callback C-Style (ctx)         ", &
    "Functor Method (OOP)           ", &
    "Reverse Communication          ", &
    "Class(*) Select Type           "  ]
  character(len=*), parameter :: plot_data_file = "build/mean_time_per_step.tsv"

  real(dp) :: y(neqn), t, h
  real(dp), target :: work(neqn, 5)
  integer :: idid, i
  type(rk_stats) :: stats

  ! Callback specifics
  real(dp) :: rpar(3)
  integer  :: ipar(1)
  type(robertson_ctx), target :: c_data
  type(robertson_params) :: params
  type(c_ptr) :: p_data
  type(robertson_functor) :: sys
  integer(8) :: t1, t2, count_rate
  real(dp)   :: elapsed, mean_elapsed
  real(dp)   :: elapsed_all(6)   ! store mean elapsed time per integration
  real(dp)   :: us_per_step_all(6)
  integer    :: plot_unit, io_stat

  work = 0.0_dp

  call system_clock(count_rate=count_rate)

  write(*,'(A)') "RK23 Final Refactored Benchmark (Clean FSAL Property)"
  write(*,'(A,I0)') "Integrations per test: ", N_runs
  write(*,'(A)') repeat("-", 80)
  write(*,'(A4,A30,A10,A8,A8,A8,A12,A12)') &
    "", "Interface", "Mean(s)", "Steps", "Rej", "NFev", "us/step", "us/NFev"
  write(*,'(A)') repeat("-", 80)

  ! ----------------------------------------------------------------------------
  ! 1. F77-Style External Callback (implicit interface)
  !    rk23_simple uses "external fun" – no explicit interface block.
  !    We still pass an internal procedure per Fortran 2008 semantics.
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = 1.0e-3_dp
    call rk23_simple(neqn, rhs_internal, t, y, t_end, h, atol, rtol, work, idid, stats)
    if (idid == -1) write(error_unit,*) "1. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(1) = mean_elapsed
  call print_row(1, plot_labels(1), mean_elapsed, stats, y, us_per_step_all(1))


  ! ----------------------------------------------------------------------------
  ! 2. SLATEC-like Callback (RPAR/IPAR)
  ! ----------------------------------------------------------------------------
  rpar = [0.04_dp, 1.0e4_dp, 3.0e7_dp]
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = 1.0e-3_dp
    call rk23_par(neqn, rob_par, t, y, t_end, h, atol, rtol, work, rpar, ipar, idid, stats)
    if (idid == -1) write(error_unit,*) "2. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(2) = mean_elapsed
  call print_row(2, plot_labels(2), mean_elapsed, stats, y, us_per_step_all(2))


  ! ----------------------------------------------------------------------------
  ! 3. C-Style Pointer
  ! ----------------------------------------------------------------------------
  c_data%k1 = 0.04_dp; c_data%k2 = 1.0e4_dp; c_data%k3 = 3.0e7_dp
  p_data = c_loc(c_data)
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = 1.0e-3_dp
    call rk23_cptr(neqn, rob_cptr, t, y, t_end, h, atol, rtol, work, p_data, idid, stats)
    if (idid == -1) write(error_unit,*) "3. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(3) = mean_elapsed
  call print_row(3, plot_labels(3), mean_elapsed, stats, y, us_per_step_all(3))


  ! ----------------------------------------------------------------------------
  ! 4. Functor OOP
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = 1.0e-3_dp
    call rk23_tb(neqn, sys, t, y, t_end, h, atol, rtol, work, idid, stats)
    if (idid == -1) write(error_unit,*) "4. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(4) = mean_elapsed
  call print_row(4, plot_labels(4), mean_elapsed, stats, y, us_per_step_all(4))


  ! ----------------------------------------------------------------------------
  ! 5. Full RCI State-Machine Loop
  ! ----------------------------------------------------------------------------
  block
    integer  :: stage
    real(dp) :: t_eval, y_eval(neqn)

    call system_clock(t1)
    do i = 1, N_runs
      t = t_start; y = y_init; h = 1.0e-3_dp
      work = 0.0_dp
      idid = 0
      stats = rk_stats()

      ! Evaluate the initial k1 directly into column 1
      call rob_direct(neqn, t, y, work(:, 1))
      stats%nfev = stats%nfev + 1
      stage = 1

      integrate: do
        call rk23_rci(stage, neqn, t, y, t_end, h, atol, rtol, work, &
                      t_eval, y_eval, idid, stats)

        select case(stage)
        case(2:4)
          call rob_direct(neqn, t_eval, y_eval, work(:, stage))
        case(6)
          if (idid == -1) write(error_unit, '(A)') "5. Step-size underflow!"
          exit integrate
        case default
          write(error_unit,'(A,I0)') "rk23_rci unexpected stage = ", stage
          error stop 1
        end select
      end do integrate
    end do
    call system_clock(t2)
    elapsed = real(t2-t1, dp)/real(count_rate, dp)
    mean_elapsed = elapsed / real(N_runs, dp)
    elapsed_all(5) = mean_elapsed
    call print_row(5, plot_labels(5), mean_elapsed, stats, y, us_per_step_all(5))
  end block

  ! ----------------------------------------------------------------------------
  ! 6. Class(*) / Unlimited Polymorphic Context
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = 1.0e-3_dp
    call rk23_class_star(neqn, rob_class_star, t, y, t_end, h, atol, rtol, work, params, idid, stats)
    if (idid == -1) write(error_unit,*) "6. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(6) = mean_elapsed
  call print_row(6, plot_labels(6), mean_elapsed, stats, y, us_per_step_all(6))

  write(*,'(A)') repeat("-", 80)
  write(*,'(A)') "Notes: Mean(s) is the mean time for one integration over all runs."
  write(*,'(A)') "       Steps and NFev are from the last run; Final Y is printed as a cross-check."
  write(*,'(A)') "       RK23 uses 3 evals per step attempt."

  open(newunit=plot_unit, file=plot_data_file, status='replace', action='write', iostat=io_stat)
  if (io_stat /= 0) then
    write(error_unit,'(A,A)') "Warning: failed to open plot data file: ", plot_data_file
  else
    write(plot_unit,'(A)') "# Interface	us_per_step"
    do i = 1, 6
      write(plot_unit,'(A,A,F12.4)') trim(plot_labels(i)), achar(9), us_per_step_all(i)
    end do
    close(plot_unit)
    write(*,'(A,A)') "Wrote machine-readable plot data: ", trim(plot_data_file)
  end if

  ! ----------------------------------------------------------------------------
  ! Callback overhead analysis
  !   penalty(i) = elapsed_all(i) / elapsed_all(1), i = 2..6
  !   Values > 1: method i is slower than the F77 baseline (overhead penalty)
  !   Values < 1: method i is faster than the F77 baseline
  !   Geometric mean of these ratios is the overall callback overhead score.
  ! ----------------------------------------------------------------------------
  block
    real(dp) :: penalty(5), log_sum, geo_mean
    integer  :: j
    character(len=30), parameter :: cb_labels(5) = [ &
      "Callback with RPAR/IPAR       ", &
      "Callback C-Style (ctx)        ", &
      "Functor Method (OOP)          ", &
      "Reverse Communication         ", &
      "Class(*) Select Type          "  ]

    write(*,'(A)') ""
    write(*,'(A)') repeat("-", 80)
    write(*,'(A)') "Callback overhead vs. F77 External (Test 1, implicit interface):"
    write(*,'(A4,A30,A12)') "", "Interface", "Penalty"
    write(*,'(A)') repeat("-", 80)

    log_sum = 0.0_dp
    do j = 1, 5
      penalty(j) = elapsed_all(j + 1) / elapsed_all(1)
      log_sum = log_sum + log(penalty(j))
      write(*,'(I2,A2,A30,F12.4)') j+1, ". ", cb_labels(j), penalty(j)
    end do

    geo_mean = exp(log_sum / 5.0_dp)
    write(*,'(A)') repeat("-", 80)
    write(*,'(A,F8.4)') &
      "Geometric mean penalty (callback overhead score): ", geo_mean
    write(*,'(A)') &
      "  score > 1.0 : callbacks slower  (overhead relative to F77)"
    write(*,'(A)') &
      "  score < 1.0 : callbacks faster than F77 implicit interface"
    write(*,'(A)') ""
    write(*,'(A)') "Note: scores may vary between runs due to runtime load, cache effects, etc."
    write(*,'(A)') "      Results must be interpreted with care."
  end block

contains

  subroutine rhs_internal(n, t_ev, y_ev, f)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: t_ev, y_ev(n)
    real(dp), intent(out) :: f(n)
    f(1) = -0.04_dp * y_ev(1) + 1.0e4_dp * y_ev(2) * y_ev(3)
    f(2) =  0.04_dp * y_ev(1) - 1.0e4_dp * y_ev(2) * y_ev(3) - 3.0e7_dp * y_ev(2)**2
    f(3) =  3.0e7_dp * y_ev(2)**2
  end subroutine rhs_internal

  subroutine print_row(id, label, mean_elapsed, s, y_fin, us_per_step_out)
    integer,          intent(in) :: id
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: mean_elapsed
    type(rk_stats),   intent(in) :: s
    real(dp),         intent(in) :: y_fin(neqn)
    real(dp),         intent(out), optional :: us_per_step_out

    real(dp) :: us_per_step, us_per_nfev

    if (s%accepted > 0) then
      us_per_step = mean_elapsed * 1.0e6_dp / real(s%accepted, dp)
    else
      us_per_step = 0.0_dp
    end if
    if (s%nfev > 0) then
      us_per_nfev = mean_elapsed * 1.0e6_dp / real(s%nfev, dp)
    else
      us_per_nfev = 0.0_dp
    end if

    write(*,'(I2,A2,A30,F10.4,I8,I8,I8,F12.4,F12.4)') &
      id, ". ", label, mean_elapsed, s%accepted, s%rejected, s%nfev, &
      us_per_step, us_per_nfev
    write(*,'(A36,A,3ES12.4)') "", "Final Y:", y_fin
    if (present(us_per_step_out)) us_per_step_out = us_per_step
  end subroutine print_row

end program rk_benchmark
