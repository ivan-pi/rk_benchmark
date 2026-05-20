program euler_benchmark
  use rk_kinds
  use rk_types
  use robertson_models
  use iso_fortran_env, only: error_unit
  implicit none

  integer,  parameter :: neqn = 3
  integer             :: N_runs, N_steps
  real(dp), parameter :: y_init(neqn) = [1.0_dp, 0.0_dp, 0.0_dp]
  real(dp), parameter :: t_start = 0.0_dp
  ! h_euler is chosen so the fixed-step count (N_steps) is well above 20 000;
  ! at h=5e-8 and t_end=1e-2 we take 200 000 steps, comparable to the ~137 000
  ! accepted steps the adaptive RK23 solver uses over [0, 100].
  real(dp), parameter :: h_euler = 5.0e-8_dp
  real(dp), parameter :: t_end   = 1.0e-2_dp   ! 200 000 steps at h=5e-8

  character(len=30), parameter :: plot_labels(6) = [ &
    "Direct RHS (rob_direct)        ", &
    "F77 Ext. (implicit iface)      ", &
    "Callback with RPAR/IPAR        ", &
    "Callback C-Style (ctx)         ", &
    "Functor Method (OOP)           ", &
    "Class(*) Select Type           "  ]
  character(len=*), parameter :: plot_data_file = "build/euler_mean_time_per_step.dat"

  real(dp) :: y(neqn), t
  integer  :: i

  ! Callback specifics
  real(dp) :: rpar(3)
  integer  :: ipar(1)
  type(robertson_ctx), target :: c_data
  type(robertson_params) :: params
  type(c_ptr) :: p_data
  type(robertson_functor) :: sys

  integer(8) :: t1, t2, count_rate
  real(dp)   :: elapsed, mean_elapsed
  real(dp)   :: elapsed_all(6)
  integer    :: plot_unit, io_stat
  logical    :: plot_data_enabled
  character(len=20) :: cli_arg
  integer    :: cli_stat
  real(dp)   :: y_ref(neqn)
  logical    :: val_ok
  integer    :: n_fail
  character(len=8) :: vstatus
  real(dp), parameter :: val_atol = 1.0e-6_dp, val_rtol = 1.0e-6_dp

  N_runs  = 500
  N_steps = nint((t_end - t_start) / h_euler)

  if (command_argument_count() >= 1) then
    call get_command_argument(1, cli_arg)
    read(cli_arg, *, iostat=cli_stat) N_runs
    if (cli_stat /= 0 .or. N_runs < 1) then
      write(error_unit,'(A)') "Usage: euler_benchmark [N_runs]"
      error stop 1
    end if
  end if

  call system_clock(count_rate=count_rate)

  ! Initialise callback context data once
  rpar = [0.04_dp, 1.0e4_dp, 3.0e7_dp]
  c_data%k1 = 0.04_dp; c_data%k2 = 1.0e4_dp; c_data%k3 = 3.0e7_dp
  p_data = c_loc(c_data)

  plot_data_enabled = .false.
  open(newunit=plot_unit, file=plot_data_file, status='replace', action='write', iostat=io_stat)
  if (io_stat /= 0) then
    write(error_unit,'(A,A)') "Warning: failed to open plot data file: ", plot_data_file
  else
    plot_data_enabled = .true.
    write(plot_unit,'(A)') "# id us_per_step label"
  end if

  write(*,'(A)') "Fixed-Step Forward Euler Benchmark (Minimal Solver Overhead)"
  write(*,'(A,I0)') "Integrations per test: ", N_runs
  write(*,'(A,G0.5)') "Euler step size h    : ", h_euler
  write(*,'(A,I0)')     "Steps per integration: ", N_steps

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

  ! 1 – reference run: direct inline
  t = t_start; y = y_init
  call run_euler_direct(t, y)
  y_ref = y
  write(*,'(I2,A2,A30,A8,3ES12.4)') 1, ". ", plot_labels(1), "REF     ", y_ref

  ! 2 – F77-style external callback
  t = t_start; y = y_init
  call run_euler_f77(t, y)
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 2, ". ", plot_labels(2), vstatus, y

  ! 3 – RPAR/IPAR callback
  t = t_start; y = y_init
  call run_euler_par(t, y)
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 3, ". ", plot_labels(3), vstatus, y

  ! 4 – C-style pointer
  t = t_start; y = y_init
  call run_euler_cptr(t, y)
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 4, ". ", plot_labels(4), vstatus, y

  ! 5 – Functor OOP
  t = t_start; y = y_init
  call run_euler_tb(t, y)
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 5, ". ", plot_labels(5), vstatus, y

  ! 6 – Class(*) Select Type
  t = t_start; y = y_init
  call run_euler_class_star(t, y)
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
  write(*,'(A4,A30,A10,A12)') "", "Interface", "Mean(s)", "us/step"
  write(*,'(A)') repeat("-", 80)

  ! 1. Direct RHS call (no callback overhead)
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    call run_euler_direct(t, y)
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(1) = mean_elapsed
  call print_row(1, plot_labels(1), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  ! 2. F77-Style External Callback (implicit interface)
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    call run_euler_f77(t, y)
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(2) = mean_elapsed
  call print_row(2, plot_labels(2), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  ! 3. SLATEC-like Callback (RPAR/IPAR)
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    call run_euler_par(t, y)
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(3) = mean_elapsed
  call print_row(3, plot_labels(3), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  ! 4. C-Style Pointer
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    call run_euler_cptr(t, y)
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(4) = mean_elapsed
  call print_row(4, plot_labels(4), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  ! 5. Functor OOP
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    call run_euler_tb(t, y)
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(5) = mean_elapsed
  call print_row(5, plot_labels(5), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  ! 6. Class(*) / Unlimited Polymorphic Context
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    call run_euler_class_star(t, y)
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(6) = mean_elapsed
  call print_row(6, plot_labels(6), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  write(*,'(A)') repeat("-", 80)
  write(*,'(A)') "Notes: Mean(s) is the mean time for one integration over all runs."
  write(*,'(A)') "       Each integration uses a fixed step-size (no adaptivity)."
  if (plot_data_enabled) then
    close(plot_unit)
    write(*,'(A,A)') "Wrote machine-readable plot data: ", trim(plot_data_file)
  end if

  ! ----------------------------------------------------------------------------
  ! Callback overhead analysis vs. Direct RHS call (strategy 1)
  ! ----------------------------------------------------------------------------
  block
    real(dp) :: penalty(5), log_sum, geo_mean
    integer  :: j
    character(len=30), parameter :: cb_labels(5) = [ &
      "F77 Ext. (implicit iface)     ", &
      "Callback with RPAR/IPAR       ", &
      "Callback C-Style (ctx)        ", &
      "Functor Method (OOP)          ", &
      "Class(*) Select Type          "  ]

    write(*,'(A)') ""
    write(*,'(A)') repeat("-", 80)
    write(*,'(A)') "Callback overhead vs. Direct RHS call (Test 1, no callback):"
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
      "  score > 1.0 : callbacks slower  (overhead relative to direct RHS call)"
    write(*,'(A)') &
      "  score < 1.0 : callbacks faster than direct RHS call"
    write(*,'(A)') ""
    write(*,'(A)') "Note: scores may vary between runs due to runtime load, cache effects, etc."
    write(*,'(A)') "      Results must be interpreted with care."
  end block

contains

  pure logical function is_close(y, y_ref, atol_v, rtol_v)
    real(dp), intent(in) :: y(:), y_ref(:), atol_v, rtol_v
    is_close = all(abs(y - y_ref) <= atol_v + rtol_v * abs(y_ref))
  end function is_close

  ! ─── Pseudo-integrators: one per callback style ──────────────────────────
  ! Each subroutine advances (t, y) by N_steps fixed Euler steps of size h_euler.
  ! They are called identically from both the validation pass and the benchmark
  ! pass, eliminating duplication between the two phases.

  subroutine run_euler_direct(t, y)
    real(dp), intent(inout) :: t, y(neqn)
    real(dp) :: dy(neqn)
    integer  :: rep
    do rep = 1, N_steps
      call rob_direct(neqn, t, y, dy)
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end subroutine run_euler_direct

  subroutine run_euler_f77(t, y)
    real(dp), intent(inout) :: t, y(neqn)
    real(dp) :: dy(neqn)
    integer  :: rep
    do rep = 1, N_steps
      call rob_direct(neqn, t, y, dy)
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end subroutine run_euler_f77

  subroutine run_euler_par(t, y)
    real(dp), intent(inout) :: t, y(neqn)
    real(dp) :: dy(neqn)
    integer  :: rep
    do rep = 1, N_steps
      call rob_par(neqn, t, y, dy, rpar, ipar)
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end subroutine run_euler_par

  subroutine run_euler_cptr(t, y)
    real(dp), intent(inout) :: t, y(neqn)
    real(dp) :: dy(neqn)
    integer  :: rep
    do rep = 1, N_steps
      call rob_cptr(neqn, t, y, dy, p_data)
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end subroutine run_euler_cptr

  subroutine run_euler_tb(t, y)
    real(dp), intent(inout) :: t, y(neqn)
    real(dp) :: dy(neqn)
    integer  :: rep
    do rep = 1, N_steps
      call sys%eval(neqn, t, y, dy)
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end subroutine run_euler_tb

  subroutine run_euler_class_star(t, y)
    real(dp), intent(inout) :: t, y(neqn)
    real(dp) :: dy(neqn)
    integer  :: rep
    do rep = 1, N_steps
      call rob_class_star(neqn, t, y, dy, params)
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end subroutine run_euler_class_star

  subroutine print_row(id, label, mean_elapsed, nsteps, plot_unit_out, plot_data_enabled_out)
    integer,          intent(in) :: id
    character(len=*), intent(in) :: label
    real(dp),         intent(in) :: mean_elapsed
    integer,          intent(in) :: nsteps
    integer,          intent(in), optional :: plot_unit_out
    logical,          intent(in), optional :: plot_data_enabled_out

    real(dp) :: us_per_step

    if (nsteps > 0) then
      us_per_step = mean_elapsed * 1.0e6_dp / real(nsteps, dp)
    else
      us_per_step = 0.0_dp
    end if

    write(*,'(I2,A2,A30,1X,G12.5,1X,G12.5)') id, ". ", label, mean_elapsed, us_per_step
    if (present(plot_unit_out) .and. present(plot_data_enabled_out)) then
      if (plot_data_enabled_out) then
        write(plot_unit_out,'(I2,1X,G0.12,1X,A)') id, us_per_step, '"'//trim(label)//'"'
      end if
    end if
  end subroutine print_row

end program euler_benchmark
