program rk_benchmark
  use rk_kinds
  use rk_types
  use rk_solvers
  use robertson_models
  use plugin_path_mod
  use iso_c_binding,   only: c_char, c_funptr, c_null_char, c_f_procpointer, c_associated
  use iso_fortran_env, only: error_unit
  implicit none

  ! Interfaces for the C dlopen shim (dlopen_iface.c)
  interface
    function rkb_dlopen(path) result(handle) bind(c, name="rkb_dlopen")
      import :: c_ptr, c_char
      character(kind=c_char), intent(in) :: path(*)
      type(c_ptr) :: handle
    end function
    function rkb_dlsym(handle, sym) result(fptr) bind(c, name="rkb_dlsym")
      import :: c_ptr, c_funptr, c_char
      type(c_ptr),            value      :: handle
      character(kind=c_char), intent(in) :: sym(*)
      type(c_funptr) :: fptr
    end function
    subroutine rkb_dlclose(handle) bind(c, name="rkb_dlclose")
      import :: c_ptr
      type(c_ptr), value :: handle
    end subroutine
  end interface

  integer,  parameter :: neqn = 3, N_runs = 100
  real(dp), parameter :: y_init(neqn) = [1.0_dp, 0.0_dp, 0.0_dp]
  real(dp), parameter :: t_start = 0.0_dp, t_end = 100.0_dp
  real(dp), parameter :: atol(neqn) = 1.0e-8_dp, rtol = 1.0e-8_dp

  real(dp) :: y(neqn), t, h
  real(dp), target :: work(neqn, 5)
  integer :: idid, i

  ! Callback specifics
  real(dp) :: rpar(3)
  integer  :: ipar(1)
  type(robertson_ctx), target :: c_data
  type(robertson_params) :: params
  type(c_ptr) :: p_data
  type(robertson_functor) :: sys
  integer(8) :: t1, t2, count_rate

  ! dlopen benchmark variables (case 7)
  type(c_ptr)  :: dl_handle
  type(c_funptr) :: raw_fp
  procedure(func_cptr), pointer :: dyn_rob => null()

  work = 0.0_dp

  call system_clock(count_rate=count_rate)

  write(*,'(A)') "RK23 Final Refactored Benchmark (Clean FSAL Property)"
  write(*,'(A,I0)') "Integrations per test: ", N_runs
  write(*,'(A)') "--------------------------------------------------------------------------------"

  ! ----------------------------------------------------------------------------
  ! 1. Internal Procedure (Capturing host data)
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = 1.0e-3_dp
    call rk23_simple(neqn, rhs_internal, t, y, t_end, h, atol, rtol, work, idid)
    if (idid == -1) write(error_unit,*) "1. Step-size underflow!"
  end do
  call system_clock(t2)
  print '(A, F10.4, A, 3F12.4)', "1. Internal Proc (Host-Data): ", &
        real(t2-t1, dp)/real(count_rate, dp), " s | Final Y:", y


  ! ----------------------------------------------------------------------------
  ! 2. SLATEC-like Callback (RPAR/IPAR)
  ! ----------------------------------------------------------------------------
  rpar = [0.04_dp, 1.0e4_dp, 3.0e7_dp]
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = 1.0e-3_dp
    call rk23_par(neqn, rob_par, t, y, t_end, h, atol, rtol, work, rpar, ipar, idid)
    if (idid == -1) write(error_unit,*) "2. Step-size underflow!"
  end do
  call system_clock(t2)
  print '(A, F10.4, A, 3F12.4)', "2. Callback with RPAR/IPAR:   ", &
        real(t2-t1, dp)/real(count_rate, dp), " s | Final Y:", y


  ! ----------------------------------------------------------------------------
  ! 3. C-Style Pointer
  ! ----------------------------------------------------------------------------
  c_data%k1 = 0.04_dp; c_data%k2 = 1.0e4_dp; c_data%k3 = 3.0e7_dp
  p_data = c_loc(c_data)
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = 1.0e-3_dp
    call rk23_cptr(neqn, rob_cptr, t, y, t_end, h, atol, rtol, work, p_data, idid)
    if (idid == -1) write(error_unit,*) "3. Step-size underflow!"
  end do
  call system_clock(t2)
  print '(A, F10.4, A, 3F12.4)', "3. Callback C-Style (ctx):    ", &
        real(t2-t1, dp)/real(count_rate, dp), " s | Final Y:", y


  ! ----------------------------------------------------------------------------
  ! 4. Functor OOP
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = 1.0e-3_dp
    call rk23_tb(neqn, sys, t, y, t_end, h, atol, rtol, work, idid)
    if (idid == -1) write(error_unit,*) "4. Step-size underflow!"
  end do
  call system_clock(t2)
  print '(A, F10.4, A, 3F12.4)', "4. Functor Method (OOP):      ", &
        real(t2-t1, dp)/real(count_rate, dp), " s | Final Y:", y


  ! ----------------------------------------------------------------------------
  ! 5. Full RCI State-Machine Loop
  ! ----------------------------------------------------------------------------
  block
    ! RCI State Variables
    integer  :: stage
    real(dp) :: t_eval, y_eval(neqn)

    call system_clock(t1)
    do i = 1, N_runs
      t = t_start; y = y_init; h = 1.0e-3_dp
      work = 0.0_dp
      idid = 0

      ! Evaluate the initial k1 directly into column 1
      call rob_direct(neqn, t, y, work(:, 1))
      stage = 1

      integrate: do
        call rk23_rci(stage, neqn, t, y, t_end, h, atol, rtol, work, &
                      t_eval, y_eval, idid)

        select case(stage)
        case(2:4)
          ! The stage requested perfectly maps to the workspace column
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
    print '(A, F10.4, A, 3F12.4)', "5. Reverse Communication:     ", &
          real(t2-t1, dp)/real(count_rate, dp), " s | Final Y:", y
  end block

! ----------------------------------------------------------------------------
  ! 6. Class(*) / Unlimited Polymorphic Context
  ! ----------------------------------------------------------------------------
  call system_clock(t1)
  do i = 1, N_runs
    t = t_start; y = y_init; h = 1.0e-3_dp
    call rk23_class_star(neqn, rob_class_star, t, y, t_end, h, atol, rtol, work, params, idid)
    if (idid == -1) write(error_unit,*) "6. Step-size underflow!"
  end do
  call system_clock(t2)
  print '(A, F10.4, A, 3F12.4)', "6. Class(*) Select Type:      ", &
        real(t2-t1, dp)/real(count_rate, dp), " s | Final Y:", y


  ! ----------------------------------------------------------------------------
  ! 7. dlopen / Dynamic Loading  (compare with case 3 – same RHS, same ctx)
  ! The function pointer is obtained at runtime via dlsym instead of being
  ! resolved at link time; this measures any overhead from indirect dispatch.
  ! ----------------------------------------------------------------------------
  dl_handle = rkb_dlopen(robertson_plugin_path // c_null_char)
  if (.not. c_associated(dl_handle)) then
    write(error_unit,'(A)') "7. dlopen failed - skipping"
  else
    raw_fp = rkb_dlsym(dl_handle, "rob_cptr" // c_null_char)
    call c_f_procpointer(raw_fp, dyn_rob)

    call system_clock(t1)
    do i = 1, N_runs
      t = t_start; y = y_init; h = 1.0e-3_dp
      call rk23_cptr(neqn, dyn_rob, t, y, t_end, h, atol, rtol, work, p_data, idid)
      if (idid == -1) write(error_unit,*) "7. Step-size underflow!"
    end do
    call system_clock(t2)
    print '(A, F10.4, A, 3F12.4)', "7. dlopen Dynamic Loading:    ", &
          real(t2-t1, dp)/real(count_rate, dp), " s | Final Y:", y

    call rkb_dlclose(dl_handle)
  end if

contains

  subroutine rhs_internal(n, t_ev, y_ev, f)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: t_ev, y_ev(n)
    real(dp), intent(out) :: f(n)
    f(1) = -0.04_dp * y_ev(1) + 1.0e4_dp * y_ev(2) * y_ev(3)
    f(2) =  0.04_dp * y_ev(1) - 1.0e4_dp * y_ev(2) * y_ev(3) - 3.0e7_dp * y_ev(2)**2
    f(3) =  3.0e7_dp * y_ev(2)**2
  end subroutine rhs_internal

end program rk_benchmark
