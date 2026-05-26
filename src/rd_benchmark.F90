! ==============================================================================
! RK23 Reaction-Diffusion Benchmark
!
! Solves the 1-D finite-difference reaction-diffusion ODE system
!
!   du_i/dt = D*(u_{i-1} - 2*u_i + u_{i+1})/dx^2 - lambda*u_i   i=1..N
!
! (N interior points, Dirichlet BCs, dx = 1/(N+1)) using adaptive RK23
! with six different callback strategies, and measures the mean time
! per accepted step.
!
! The initial condition u_i(0) = sin(pi*i/(N+1)) is the exact lowest
! eigenmode of the discrete operator, so the exact time-marching solution
! is known analytically and the true time-integration error can be measured.
!
! Usage:  rd_benchmark [N_runs [N_grid]]
!   N_runs  : number of repeated integrations per strategy (default 50)
!   N_grid  : number of interior grid points / ODE unknowns (default 100)
! ==============================================================================
program rd_benchmark
  use rk_kinds
  use rk_types
  use rk_solvers, only: rk23_simple, rk23_par, rk23_tb, rk23_rci, rk23_class_star, rk_stats
#ifdef USE_EXTERNAL_C_RK23
  use rk_solvers, only: rk23_cptr => rk23_cptr_external
#else
  use rk_solvers, only: rk23_cptr
#endif
  use reaction_diffusion_models
  use iso_c_binding, only: c_funloc
  use iso_fortran_env, only: error_unit
  implicit none

  integer  :: neqn             ! number of interior grid points (= ODE size)
  integer  :: N_runs
  real(dp) :: t_start, t_end, h_init
  real(dp) :: atol_val, rtol
  character(len=*), parameter :: plot_data_file = "build/rd_mean_time_per_step.dat"
  character(len=30), parameter :: plot_labels(6) = [ &
    "F77 Ext. (implicit iface)      ", &
    "Callback with RPAR/IPAR        ", &
    "Callback C-Style (ctx)         ", &
    "Functor Method (OOP)           ", &
    "Reverse Communication          ", &
    "Class(*) Select Type           "  ]

  ! Allocatable state arrays (sized by neqn at runtime)
  real(dp), allocatable :: y(:), work(:,:), atol(:), y_init_arr(:), y_ref(:), y_exact_arr(:)

  real(dp) :: t, h
  integer  :: idid, i
  type(rk_stats) :: stats

  ! Callback-specific contexts
  real(dp) :: rpar(2)
  integer  :: ipar(1)
  type(rd_ctx), target :: c_data
  type(rd_params) :: params
  type(c_ptr) :: p_data
  type(rd_functor) :: sys

  integer(8) :: t1, t2, count_rate
  real(dp)   :: elapsed, mean_elapsed
  real(dp)   :: elapsed_all(6)
  integer    :: plot_unit, io_stat
  logical    :: plot_data_enabled
  character(len=20) :: cli_arg
  integer    :: cli_stat
  logical    :: val_ok
  integer    :: n_fail
  character(len=8) :: vstatus
  real(dp), parameter :: val_atol = 1.0e-5_dp, val_rtol = 1.0e-5_dp

  ! ============================================================================
  ! Parse command-line arguments
  ! ============================================================================
  N_runs   = 50
  neqn     = 100
  t_start  = 0.0_dp
  t_end    = 0.5_dp
  atol_val = 1.0e-6_dp
  rtol     = 1.0e-6_dp
  h_init   = 1.0e-4_dp

  if (command_argument_count() >= 1) then
    call get_command_argument(1, cli_arg)
    read(cli_arg, *, iostat=cli_stat) N_runs
    if (cli_stat /= 0 .or. N_runs < 1) then
      write(error_unit,'(A)') "Usage: rd_benchmark [N_runs [N_grid]]"
      error stop 1
    end if
  end if

  if (command_argument_count() >= 2) then
    call get_command_argument(2, cli_arg)
    read(cli_arg, *, iostat=cli_stat) neqn
    if (cli_stat /= 0 .or. neqn < 1) then
      write(error_unit,'(A)') "Usage: rd_benchmark [N_runs [N_grid]]"
      error stop 1
    end if
  end if

  ! ============================================================================
  ! Allocate and initialise arrays
  ! ============================================================================
  allocate(y(neqn), work(neqn,5), atol(neqn), y_init_arr(neqn), &
           y_ref(neqn), y_exact_arr(neqn))

  atol = atol_val

  ! Initial condition: lowest eigenmode sin(pi * x_i), x_i = i/(neqn+1)
  block
    real(dp) :: pi
    integer  :: j
    pi = acos(-1.0_dp)
    do j = 1, neqn
      y_init_arr(j) = sin(pi * real(j, dp) / real(neqn+1, dp))
    end do
  end block

  ! Set up callback contexts (all use the same D and lambda as module params)
  rpar(1) = rd_D;    rpar(2) = rd_lambda
  c_data%D = rd_D;  c_data%lambda = rd_lambda
  p_data   = c_loc(c_data)
  params%D = rd_D;  params%lambda = rd_lambda
  sys%D    = rd_D;  sys%lambda    = rd_lambda

  call system_clock(count_rate=count_rate)

  ! ============================================================================
  ! Open plot data file
  ! ============================================================================
  plot_data_enabled = .false.
  open(newunit=plot_unit, file=plot_data_file, status='replace', action='write', iostat=io_stat)
  if (io_stat /= 0) then
    write(error_unit,'(A,A)') "Warning: failed to open plot data file: ", plot_data_file
  else
    plot_data_enabled = .true.
    write(plot_unit,'(A)') "# id us_per_step label"
  end if

  write(*,'(A)') "RK23 Reaction-Diffusion Benchmark"
  write(*,'(A,G0.3,A,G0.3)') &
    "Problem: 1-D FD reaction-diffusion, D = ", rd_D, "  lambda = ", rd_lambda
  write(*,'(A,I0,A,I0)') "Grid points N = ", neqn, "    Integrations per test: ", N_runs
  write(*,'(A,G0.5,A,G0.5)') &
    "Interval: [", t_start, ", ", t_end
  write(*,'(A,G0.5,A,G0.5)') "Tolerances: atol = ", atol_val, "  rtol = ", rtol

  ! ============================================================================
  ! Validation / Warm-Up Pass
  ! ============================================================================
  write(*,'(A)') ""
  write(*,'(A)') "Validation / Warm-Up Pass"
  write(*,'(A)') repeat("-", 80)
  write(*,'(A4,A30,A8,2A16)') "", "Interface", "Status", "Y(1)", "Y(N/2)"
  write(*,'(A)') repeat("-", 80)
  n_fail = 0

  ! 1 – reference run (F77-style)
  work = 0.0_dp
  t = t_start; y = y_init_arr; h = h_init
  call rk23_simple(neqn, rhs_internal, t, y, t_end, h, atol, rtol, work, idid, stats)
  if (idid == -1) write(error_unit,*) "warm-up 1: step-size underflow!"
  y_ref = y
  write(*,'(I2,A2,A30,A8,2ES16.7)') 1, ". ", plot_labels(1), "REF     ", y_ref(1), y_ref(neqn/2)

  ! 2 – RPAR/IPAR
  work = 0.0_dp
  t = t_start; y = y_init_arr; h = h_init
  call rk23_par(neqn, rd_par, t, y, t_end, h, atol, rtol, work, rpar, ipar, idid, stats)
  if (idid == -1) write(error_unit,*) "warm-up 2: step-size underflow!"
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,2ES16.7)') 2, ". ", plot_labels(2), vstatus, y(1), y(neqn/2)

  ! 3 – C-pointer
  work = 0.0_dp
  t = t_start; y = y_init_arr; h = h_init
  call rk23_cptr(neqn, c_funloc(rd_cptr), t, y, t_end, h, atol, rtol, work, p_data, idid, stats)
  if (idid == -1) write(error_unit,*) "warm-up 3: step-size underflow!"
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,2ES16.7)') 3, ". ", plot_labels(3), vstatus, y(1), y(neqn/2)

  ! 4 – OOP functor
  work = 0.0_dp
  t = t_start; y = y_init_arr; h = h_init
  call rk23_tb(neqn, sys, t, y, t_end, h, atol, rtol, work, idid, stats)
  if (idid == -1) write(error_unit,*) "warm-up 4: step-size underflow!"
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,2ES16.7)') 4, ". ", plot_labels(4), vstatus, y(1), y(neqn/2)

  ! 5 – RCI
  t = t_start; y = y_init_arr; h = h_init
  call run_rci(t, y, h, idid, stats, "warm-up 5")
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,2ES16.7)') 5, ". ", plot_labels(5), vstatus, y(1), y(neqn/2)

  ! 6 – class(*)
  work = 0.0_dp
  t = t_start; y = y_init_arr; h = h_init
  call rk23_class_star(neqn, rd_class_star, t, y, t_end, h, atol, rtol, work, params, idid, stats)
  if (idid == -1) write(error_unit,*) "warm-up 6: step-size underflow!"
  val_ok = is_close(y, y_ref, val_atol, val_rtol)
  if (.not. val_ok) n_fail = n_fail + 1
  vstatus = merge("PASS    ", "FAIL    ", val_ok)
  write(*,'(I2,A2,A30,A8,2ES16.7)') 6, ". ", plot_labels(6), vstatus, y(1), y(neqn/2)

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

  ! 1. F77-Style External Callback (implicit interface)
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init_arr; h = h_init
    call rk23_simple(neqn, rhs_internal, t, y, t_end, h, atol, rtol, work, idid, stats)
    if (idid == -1) write(error_unit,*) "1. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(1) = mean_elapsed
  call print_row(1, plot_labels(1), mean_elapsed, stats, plot_unit, plot_data_enabled)

  ! 2. SLATEC-like Callback (RPAR/IPAR)
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init_arr; h = h_init
    call rk23_par(neqn, rd_par, t, y, t_end, h, atol, rtol, work, rpar, ipar, idid, stats)
    if (idid == -1) write(error_unit,*) "2. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(2) = mean_elapsed
  call print_row(2, plot_labels(2), mean_elapsed, stats, plot_unit, plot_data_enabled)

  ! 3. C-Style Pointer
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init_arr; h = h_init
    call rk23_cptr(neqn, c_funloc(rd_cptr), t, y, t_end, h, atol, rtol, work, p_data, idid, stats)
    if (idid == -1) write(error_unit,*) "3. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(3) = mean_elapsed
  call print_row(3, plot_labels(3), mean_elapsed, stats, plot_unit, plot_data_enabled)

  ! 4. Functor OOP
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init_arr; h = h_init
    call rk23_tb(neqn, sys, t, y, t_end, h, atol, rtol, work, idid, stats)
    if (idid == -1) write(error_unit,*) "4. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(4) = mean_elapsed
  call print_row(4, plot_labels(4), mean_elapsed, stats, plot_unit, plot_data_enabled)

  ! 5. Full RCI State-Machine Loop
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init_arr; h = h_init
    call run_rci(t, y, h, idid, stats, "5")
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
  mean_elapsed = elapsed / real(N_runs, dp)
  elapsed_all(5) = mean_elapsed
  call print_row(5, plot_labels(5), mean_elapsed, stats, plot_unit, plot_data_enabled)

  ! 6. Class(*) / Unlimited Polymorphic Context
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init_arr; h = h_init
    call rk23_class_star(neqn, rd_class_star, t, y, t_end, h, atol, rtol, work, params, idid, stats)
    if (idid == -1) write(error_unit,*) "6. Step-size underflow!"
  end do
  call system_clock(t2)
  elapsed = real(t2-t1, dp) / real(count_rate, dp)
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

  ! ============================================================================
  ! True error vs. exact analytical solution (from the last strategy-1 run)
  ! ============================================================================
  call rd_exact(neqn, t_end, rd_D, rd_lambda, y_exact_arr)

  block
    real(dp) :: max_err, l2_err, diff_i
    integer  :: j

    max_err = 0.0_dp
    l2_err  = 0.0_dp
    do j = 1, neqn
      diff_i  = abs(y_ref(j) - y_exact_arr(j))
      max_err = max(max_err, diff_i)
      l2_err  = l2_err + diff_i**2
    end do
    l2_err = sqrt(l2_err / real(neqn, dp))

    write(*,'(A)') ""
    write(*,'(A)') repeat("-", 80)
    write(*,'(A)') "Time-integration error vs. exact eigenmode solution (strategy 1):"
    write(*,'(A,ES10.3)') "  Max-norm error : ", max_err
    write(*,'(A,ES10.3)') "  L2-norm  error : ", l2_err
    write(*,'(A)') &
      "  (Exact solution: u_i(T) = sin(pi*i/(N+1)) * exp(-(mu1+lambda)*T))"
  end block

  ! ============================================================================
  ! Callback overhead analysis
  ! ============================================================================
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

  ! RHS callback passed to rk23_simple (F77-style, implicit interface).
  ! Delegates to the module-level rd_direct which uses the fixed parameters.
  subroutine rhs_internal(n, t_ev, y_ev, f)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: t_ev, y_ev(n)
    real(dp), intent(out) :: f(n)
    call rd_direct(n, t_ev, y_ev, f)
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

    write(*,'(I2,A2,A30,1X,G12.5,I8,I8,I8,1X,G12.5)') &
      id, ". ", label, mean_elapsed, s%accepted, s%rejected, s%nfev, &
      us_per_step
    if (present(plot_unit_out) .and. present(plot_data_enabled_out)) then
      if (plot_data_enabled_out) then
        write(plot_unit_out,'(I2,1X,G0.12,1X,A)') id, us_per_step, '"'//trim(label)//'"'
      end if
    end if
  end subroutine print_row

  ! Runs a single RK23 integration using the reverse communication interface.
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
    call rd_direct(neqn, t, y, work(:, 1))
    stats%nfev = stats%nfev + 1
    stage = 1
    step_rejected = .false.

    rci_loop: do
      call rk23_rci(stage, neqn, t, y, t_end, h, atol, rtol, work, &
                    t_eval, y_eval, idid, stats, step_rejected)
      select case(stage)
      case(2:4)
        call rd_direct(neqn, t_eval, y_eval, work(:, stage))
      case(6)
        if (idid == -1) then
          if (present(err_label)) then
            write(error_unit,'(A,A)') trim(err_label), ": step-size underflow!"
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

end program rd_benchmark
