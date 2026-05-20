program rk_benchmark
  use rk_kinds
  use rk_types
  use rk_solvers, only: rk23_simple, rk23_par, rk23_tb, rk23_rci, rk23_class_star, rk23_cptr, rk_stats
  use robertson_models
  use iso_fortran_env, only: error_unit
  implicit none

  integer,  parameter :: neqn = 3
  integer             :: N_runs
  real(dp), parameter :: y_init(neqn) = [1.0_dp, 0.0_dp, 0.0_dp]
  real(dp), parameter :: t_start = 0.0_dp, t_end = 100.0_dp
  real(dp), parameter :: atol(neqn) = 1.0e-8_dp, rtol = 1.0e-8_dp
  real(dp), parameter :: h_init = 1.0e-4_dp
  character(len=30), parameter :: plot_labels(6) = [ &
    "F77 Ext. (implicit iface)      ", &
    "Callback with RPAR/IPAR        ", &
    "Callback C-Style (ctx)         ", &
    "Functor Method (OOP)           ", &
    "Reverse Communication          ", &
    "Class(*) Select Type           "  ]
  character(len=*), parameter :: plot_data_file = "build/mean_time_per_step.dat"

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
  integer    :: plot_unit, io_stat
  logical    :: plot_data_enabled
  character(len=20) :: cli_arg
  integer    :: cli_stat
  real(dp)   :: y_ref(neqn)        ! reference final state from strategy 1 (warm-up)
  logical    :: val_ok              ! per-strategy consistency check result
  integer    :: n_fail              ! count of failed consistency checks
  character(len=8) :: vstatus       ! validation status string for output
  real(dp), parameter :: val_atol = 1.0e-6_dp, val_rtol = 1.0e-6_dp

  work = 0.0_dp

  N_runs = 100
  if (command_argument_count() >= 1) then
    call get_command_argument(1, cli_arg)
    read(cli_arg, *, iostat=cli_stat) N_runs
    if (cli_stat /= 0 .or. N_runs < 1) then
      write(error_unit,'(A)') "Usage: rk_benchmark [N_runs]"
      error stop 1
    end if
  end if

  call system_clock(count_rate=count_rate)

  plot_data_enabled = .false.
  open(newunit=plot_unit, file=plot_data_file, status='replace', action='write', iostat=io_stat)
  if (io_stat /= 0) then
    write(error_unit,'(A,A)') "Warning: failed to open plot data file: ", plot_data_file
  else
    plot_data_enabled = .true.
    write(plot_unit,'(A)') "# id us_per_step label"
  end if

  write(*,'(A)') "RK23 Final Refactored Benchmark (Clean FSAL Property)"
  write(*,'(A,I0)') "Integrations per test: ", N_runs

  ! ============================================================================
  ! Validation / Warm-Up Pass
  !   Each strategy is run once to verify consistency against strategy 1.
  !   Tolerance: |y(i) - y_ref(i)| <= val_atol + val_rtol * |y_ref(i)|
  ! ============================================================================
  write(*,'(A)') ""
  write(*,'(A)') "Validation / Warm-Up Pass"
  write(*,'(A)') repeat("-", 80)
  write(*,'(A4,A30,A8,3A12)') "", "Interface", "Status", "Y(1)", "Y(2)", "Y(3)"
  write(*,'(A)') repeat("-", 80)
  n_fail = 0

  ! 1 – reference run
  work = 0.0_dp
  t = t_start; y = y_init; h = h_init
  call rk23_simple(neqn, rhs_internal, t, y, t_end, h, atol, rtol, work, idid, stats)
  if (idid == -1) write(error_unit,*) "warm-up 1: step-size underflow!"
  y_ref = y
  write(*,'(I2,A2,A30,A8,3ES12.4)') 1, ". ", plot_labels(1), "REF     ", y_ref

  ! 2
  rpar = [0.04_dp, 1.0e4_dp, 3.0e7_dp]
  work = 0.0_dp
  t = t_start; y = y_init; h = h_init
  call rk23_par(neqn, rob_par, t, y, t_end, h, atol, rtol, work, rpar, ipar, idid, stats)
  if (idid == -1) write(error_unit,*) "warm-up 2: step-size underflow!"
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 2, ". ", plot_labels(2), vstatus, y

  ! 3
  c_data%k1 = 0.04_dp; c_data%k2 = 1.0e4_dp; c_data%k3 = 3.0e7_dp
  p_data = c_loc(c_data)
  work = 0.0_dp
  t = t_start; y = y_init; h = h_init
  call run_cptr_solver(t, y, h, idid, stats)
  if (idid == -1) write(error_unit,*) "warm-up 3: step-size underflow!"
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 3, ". ", plot_labels(3), vstatus, y

  ! 4
  work = 0.0_dp
  t = t_start; y = y_init; h = h_init
  call rk23_tb(neqn, sys, t, y, t_end, h, atol, rtol, work, idid, stats)
  if (idid == -1) write(error_unit,*) "warm-up 4: step-size underflow!"
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 4, ". ", plot_labels(4), vstatus, y

  ! 5 – RCI
  t = t_start; y = y_init; h = h_init
  call run_rci(t, y, h, idid, stats, "warm-up 5")
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 5, ". ", plot_labels(5), vstatus, y

  ! 6
  work = 0.0_dp
  t = t_start; y = y_init; h = h_init
  call rk23_class_star(neqn, rob_class_star, t, y, t_end, h, atol, rtol, work, params, idid, stats)
  if (idid == -1) write(error_unit,*) "warm-up 6: step-size underflow!"
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 6, ". ", plot_labels(6), vstatus, y

  write(*,'(A)') repeat("-", 80)
  if (n_fail == 0) then
    write(*,'(A)') "All strategies consistent with strategy 1 (within tolerance)."
  else
    write(*,'(I0,A)') n_fail, " strategy/strategies failed the consistency check!"
  end if

  ! ============================================================================
  ! Benchmark Pass
  ! ============================================================================
  write(*,'(A)') ""
  write(*,'(A)') "Benchmark Pass"
  write(*,'(A)') repeat("-", 80)
  write(*,'(A4,A30,A10,A8,A8,A8,A12)') &
    "", "Interface", "Mean(s)", "Steps", "Rej", "NFev", "us/step"
  write(*,'(A)') repeat("-", 80)

  ! ----------------------------------------------------------------------------
  ! 1. F77-Style External Callback (implicit interface)
  !    rk23_simple uses "external fun" – no explicit interface block.
  !    We still pass an internal procedure per Fortran 2008 semantics.
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = h_init
    call rk23_simple(neqn, rhs_internal, t, y, t_end, h, atol, rtol, work, idid, stats)
    if (idid == -1) write(error_unit,*) "1. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(1) = mean_elapsed
  call print_row(1, plot_labels(1), mean_elapsed, stats, plot_unit, plot_data_enabled)


  ! ----------------------------------------------------------------------------
  ! 2. SLATEC-like Callback (RPAR/IPAR)
  ! ----------------------------------------------------------------------------
  rpar = [0.04_dp, 1.0e4_dp, 3.0e7_dp]
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = h_init
    call rk23_par(neqn, rob_par, t, y, t_end, h, atol, rtol, work, rpar, ipar, idid, stats)
    if (idid == -1) write(error_unit,*) "2. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(2) = mean_elapsed
  call print_row(2, plot_labels(2), mean_elapsed, stats, plot_unit, plot_data_enabled)


  ! ----------------------------------------------------------------------------
  ! 3. C-Style Pointer
  ! ----------------------------------------------------------------------------
  c_data%k1 = 0.04_dp; c_data%k2 = 1.0e4_dp; c_data%k3 = 3.0e7_dp
  p_data = c_loc(c_data)
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = h_init
    call run_cptr_solver(t, y, h, idid, stats)
    if (idid == -1) write(error_unit,*) "3. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(3) = mean_elapsed
  call print_row(3, plot_labels(3), mean_elapsed, stats, plot_unit, plot_data_enabled)


  ! ----------------------------------------------------------------------------
  ! 4. Functor OOP
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = h_init
    call rk23_tb(neqn, sys, t, y, t_end, h, atol, rtol, work, idid, stats)
    if (idid == -1) write(error_unit,*) "4. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(4) = mean_elapsed
  call print_row(4, plot_labels(4), mean_elapsed, stats, plot_unit, plot_data_enabled)


  ! ----------------------------------------------------------------------------
  ! 5. Full RCI State-Machine Loop
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = h_init
    call run_rci(t, y, h, idid, stats, "5")
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(5) = mean_elapsed
  call print_row(5, plot_labels(5), mean_elapsed, stats, plot_unit, plot_data_enabled)

  ! ----------------------------------------------------------------------------
  ! 6. Class(*) / Unlimited Polymorphic Context
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = h_init
    call rk23_class_star(neqn, rob_class_star, t, y, t_end, h, atol, rtol, work, params, idid, stats)
    if (idid == -1) write(error_unit,*) "6. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp)/real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(6) = mean_elapsed
  call print_row(6, plot_labels(6), mean_elapsed, stats, plot_unit, plot_data_enabled)

  write(*,'(A)') repeat("-", 80)
  write(*,'(A)') "Notes: Mean(s) is the mean time for one integration over all runs."
  write(*,'(A)') "       Steps and NFev are from the last run."
  write(*,'(A)') "       RK23 uses 3 evals per step attempt."
  if (plot_data_enabled) then
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

  pure logical function is_close(y, y_ref, atol_v, rtol_v)
    real(dp), intent(in) :: y(:), y_ref(:), atol_v, rtol_v
    is_close = all(abs(y - y_ref) <= atol_v + rtol_v * abs(y_ref))
  end function is_close

  ! RHS callback passed to rk23_simple (F77-style, implicit interface)
  subroutine rhs_internal(n, t_ev, y_ev, f)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: t_ev, y_ev(n)
    real(dp), intent(out) :: f(n)
    f(1) = -0.04_dp * y_ev(1) + 1.0e4_dp * y_ev(2) * y_ev(3)
    f(2) =  0.04_dp * y_ev(1) - 1.0e4_dp * y_ev(2) * y_ev(3) - 3.0e7_dp * y_ev(2)**2
    f(3) =  3.0e7_dp * y_ev(2)**2
  end subroutine rhs_internal

  subroutine print_row(id, label, mean_elapsed, s, plot_unit_out, plot_data_enabled_out)
    integer,          intent(in) :: id
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: mean_elapsed
    type(rk_stats),   intent(in) :: s
    integer,          intent(in), optional :: plot_unit_out
    logical,          intent(in), optional :: plot_data_enabled_out

    real(dp) :: us_per_step

    if (s%accepted > 0) then
      us_per_step = mean_elapsed * 1.0e6_dp / real(s%accepted, dp)
    else
      us_per_step = 0.0_dp
    end if

    write(*,'(I2,A2,A30,F10.4,I8,I8,I8,F12.4)') &
      id, ". ", label, mean_elapsed, s%accepted, s%rejected, s%nfev, &
      us_per_step
    if (present(plot_unit_out) .and. present(plot_data_enabled_out)) then
      if (plot_data_enabled_out) then
        write(plot_unit_out,'(I2,1X,F12.4,1X,A)') id, us_per_step, '"'//trim(label)//'"'
      end if
    end if
  end subroutine print_row

  ! Runs a single RK23 integration using the reverse communication interface.
  ! Encapsulates the state-machine loop so it can be called from both the
  ! validation pass and the benchmark pass without duplication.
  subroutine run_cptr_solver(t, y, h, idid, stats)
#ifdef USE_EXTERNAL_C_RK23
    use rk_solvers, only: rk23_cptr_external
#endif
    real(dp), intent(inout) :: t, y(neqn), h
    integer,  intent(out)   :: idid
    type(rk_stats), intent(out) :: stats

#ifdef USE_EXTERNAL_C_RK23
    call rk23_cptr_external(neqn, rob_cptr, t, y, t_end, h, atol, rtol, work, p_data, idid, stats)
#else
    call rk23_cptr(neqn, rob_cptr, t, y, t_end, h, atol, rtol, work, p_data, idid, stats)
#endif
  end subroutine run_cptr_solver

  subroutine run_rci(t, y, h, idid, stats, err_label)
    real(dp),         intent(inout) :: t, y(neqn), h
    integer,          intent(out)   :: idid
    type(rk_stats),   intent(out)   :: stats
    character(len=*), intent(in), optional :: err_label

    integer  :: stage
    real(dp) :: t_eval, y_eval(neqn)
    logical  :: step_rejected

    work = 0.0_dp
    idid = 0
    stats = rk_stats()

    ! Evaluate the initial k1 directly into column 1
    call rob_direct(neqn, t, y, work(:, 1))
    stats%nfev = stats%nfev + 1
    stage = 1
    step_rejected = .false.

    rci_loop: do
      call rk23_rci(stage, neqn, t, y, t_end, h, atol, rtol, work, &
                    t_eval, y_eval, idid, stats, step_rejected)
      select case(stage)
      case(2:4)
        call rob_direct(neqn, t_eval, y_eval, work(:, stage))
      case(6)
        if (idid == -1) then
          if (present(err_label)) then
            write(error_unit,'(A,A)') err_label, ": step-size underflow!"
          else
            write(error_unit,'(A)') "rk23_rci: step-size underflow!"
          end if
        end if
        exit rci_loop
      case default
        write(error_unit,'(A,I0)') "rk23_rci unexpected stage = ", stage
        error stop 1
      end select
    end do rci_loop
  end subroutine run_rci

end program rk_benchmark
