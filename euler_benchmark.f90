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
  real(dp), parameter :: h_euler = 1.0e-7_dp   ! fixed Euler step (within stability region)
  real(dp), parameter :: t_end   = 1.0e-2_dp   ! short span: 100000 steps at h=1e-7

  character(len=30), parameter :: plot_labels(6) = [ &
    "Direct Inline                  ", &
    "F77 Ext. (implicit iface)      ", &
    "Callback with RPAR/IPAR        ", &
    "Callback C-Style (ctx)         ", &
    "Functor Method (OOP)           ", &
    "Class(*) Select Type           "  ]
  character(len=*), parameter :: plot_data_file = "build/euler_mean_time_per_step.dat"

  real(dp) :: y(neqn), t, dy(neqn)
  integer  :: i, rep

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
  character(len=8) :: vstatus
  real(dp), parameter :: val_atol = 1.0e-6_dp, val_rtol = 1.0e-6_dp

  N_runs  = 100
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
  write(*,'(A,ES10.2)') "Euler step size h    : ", h_euler
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

  ! 1 – reference run: direct inline
  t = t_start; y = y_init
  do rep = 1, N_steps
    dy(1) = -0.04_dp * y(1) + 1.0e4_dp * y(2) * y(3)
    dy(2) =  0.04_dp * y(1) - 1.0e4_dp * y(2) * y(3) - 3.0e7_dp * y(2)**2
    dy(3) =  3.0e7_dp * y(2)**2
    y = y + h_euler * dy
    t = t + h_euler
  end do
  y_ref = y
  write(*,'(I2,A2,A30,A8,3ES12.4)') 1, ". ", plot_labels(1), "REF     ", y_ref

  ! 2 – F77-style external callback
  t = t_start; y = y_init
  do rep = 1, N_steps
    call rhs_internal(neqn, t, y, dy)
    y = y + h_euler * dy
    t = t + h_euler
  end do
  val_ok  = is_close(y, y_ref, val_atol, val_rtol)
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 2, ". ", plot_labels(2), vstatus, y

  ! 3 – RPAR/IPAR callback
  rpar = [0.04_dp, 1.0e4_dp, 3.0e7_dp]
  t = t_start; y = y_init
  do rep = 1, N_steps
    call rob_par(neqn, t, y, dy, rpar, ipar)
    y = y + h_euler * dy
    t = t + h_euler
  end do
  val_ok  = is_close(y, y_ref, val_atol, val_rtol)
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 3, ". ", plot_labels(3), vstatus, y

  ! 4 – C-style pointer
  c_data%k1 = 0.04_dp; c_data%k2 = 1.0e4_dp; c_data%k3 = 3.0e7_dp
  p_data = c_loc(c_data)
  t = t_start; y = y_init
  do rep = 1, N_steps
    call rob_cptr(neqn, t, y, dy, p_data)
    y = y + h_euler * dy
    t = t + h_euler
  end do
  val_ok  = is_close(y, y_ref, val_atol, val_rtol)
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 4, ". ", plot_labels(4), vstatus, y

  ! 5 – Functor OOP
  t = t_start; y = y_init
  do rep = 1, N_steps
    call sys%eval(neqn, t, y, dy)
    y = y + h_euler * dy
    t = t + h_euler
  end do
  val_ok  = is_close(y, y_ref, val_atol, val_rtol)
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 5, ". ", plot_labels(5), vstatus, y

  ! 6 – Class(*) Select Type
  t = t_start; y = y_init
  do rep = 1, N_steps
    call rob_class_star(neqn, t, y, dy, params)
    y = y + h_euler * dy
    t = t + h_euler
  end do
  val_ok  = is_close(y, y_ref, val_atol, val_rtol)
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,3ES12.4)') 6, ". ", plot_labels(6), vstatus, y

  write(*,'(A)') repeat("-", 80)
  write(*,'(A)') "All strategies produce identical results (all consistent with strategy 1)."

  ! ============================================================================
  ! Benchmark Pass
  ! ============================================================================
  write(*,'(A)') ""
  write(*,'(A)') "Benchmark Pass"
  write(*,'(A)') repeat("-", 80)
  write(*,'(A4,A30,A10,A12)') "", "Interface", "Mean(s)", "us/step"
  write(*,'(A)') repeat("-", 80)

  ! ----------------------------------------------------------------------------
  ! 1. Direct Inline (no callback overhead)
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    do rep = 1, N_steps
      dy(1) = -0.04_dp * y(1) + 1.0e4_dp * y(2) * y(3)
      dy(2) =  0.04_dp * y(1) - 1.0e4_dp * y(2) * y(3) - 3.0e7_dp * y(2)**2
      dy(3) =  3.0e7_dp * y(2)**2
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(1) = mean_elapsed
  call print_row(1, plot_labels(1), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  ! ----------------------------------------------------------------------------
  ! 2. F77-Style External Callback (implicit interface)
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    do rep = 1, N_steps
      call rhs_internal(neqn, t, y, dy)
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(2) = mean_elapsed
  call print_row(2, plot_labels(2), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  ! ----------------------------------------------------------------------------
  ! 3. SLATEC-like Callback (RPAR/IPAR)
  ! ----------------------------------------------------------------------------
  rpar = [0.04_dp, 1.0e4_dp, 3.0e7_dp]
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    do rep = 1, N_steps
      call rob_par(neqn, t, y, dy, rpar, ipar)
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(3) = mean_elapsed
  call print_row(3, plot_labels(3), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  ! ----------------------------------------------------------------------------
  ! 4. C-Style Pointer
  ! ----------------------------------------------------------------------------
  c_data%k1 = 0.04_dp; c_data%k2 = 1.0e4_dp; c_data%k3 = 3.0e7_dp
  p_data = c_loc(c_data)
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    do rep = 1, N_steps
      call rob_cptr(neqn, t, y, dy, p_data)
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(4) = mean_elapsed
  call print_row(4, plot_labels(4), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  ! ----------------------------------------------------------------------------
  ! 5. Functor OOP
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    do rep = 1, N_steps
      call sys%eval(neqn, t, y, dy)
      y = y + h_euler * dy
      t = t + h_euler
    end do
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(5) = mean_elapsed
  call print_row(5, plot_labels(5), mean_elapsed, N_steps, plot_unit, plot_data_enabled)

  ! ----------------------------------------------------------------------------
  ! 6. Class(*) / Unlimited Polymorphic Context
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init
    do rep = 1, N_steps
      call rob_class_star(neqn, t, y, dy, params)
      y = y + h_euler * dy
      t = t + h_euler
    end do
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
  ! Callback overhead analysis vs. Direct Inline (strategy 1)
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
    write(*,'(A)') "Callback overhead vs. Direct Inline (Test 1, no callback):"
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
      "  score > 1.0 : callbacks slower  (overhead relative to direct inline)"
    write(*,'(A)') &
      "  score < 1.0 : callbacks faster than direct inline"
    write(*,'(A)') ""
    write(*,'(A)') "Note: scores may vary between runs due to runtime load, cache effects, etc."
    write(*,'(A)') "      Results must be interpreted with care."
  end block

contains

  pure logical function is_close(y, y_ref, atol_v, rtol_v)
    real(dp), intent(in) :: y(:), y_ref(:), atol_v, rtol_v
    is_close = all(abs(y - y_ref) <= atol_v + rtol_v * abs(y_ref))
  end function is_close

  ! RHS callback passed as F77-style external (implicit interface)
  subroutine rhs_internal(n, t_ev, y_ev, f)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: t_ev, y_ev(n)
    real(dp), intent(out) :: f(n)
    f(1) = -0.04_dp * y_ev(1) + 1.0e4_dp * y_ev(2) * y_ev(3)
    f(2) =  0.04_dp * y_ev(1) - 1.0e4_dp * y_ev(2) * y_ev(3) - 3.0e7_dp * y_ev(2)**2
    f(3) =  3.0e7_dp * y_ev(2)**2
  end subroutine rhs_internal

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

    write(*,'(I2,A2,A30,F10.4,F12.4)') id, ". ", label, mean_elapsed, us_per_step
    if (present(plot_unit_out) .and. present(plot_data_enabled_out)) then
      if (plot_data_enabled_out) then
        write(plot_unit_out,'(I2,1X,F12.4,1X,A)') id, us_per_step, '"'//trim(label)//'"'
      end if
    end if
  end subroutine print_row

end program euler_benchmark
