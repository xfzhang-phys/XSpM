MODULE spm_util
  USE lib_math, ONLY : DP
  IMPLICIT NONE

  !
  ! SpM Input Parameters
  !
  ! Statistic type: Fermion/Boson
  CHARACTER(LEN=1) :: statistics = 'F'
  ! Size of tau-grid on the imaginary-time axis
  INTEGER :: ntau = 4000
  ! Size of omega-grid on the real-frequency axis
  INTEGER :: nw = 1000
  ! Inverse temperature
  REAL(DP) :: beta = 100
  ! Boundaries of the energy window
  REAL(DP) :: wmin = -4.0D0
  REAL(DP) :: wmax = 4.0D0
  ! Threshold for cutting off the SVD basis size
  REAL(DP) :: svd_cut = 1.0D-10
  ! If auto update the penalties in ADMM algorithm
  LOGICAL :: admm_auto_penalty = .FALSE.
  ! Maximal interations for the convergence of ADMM algorithm
  INTEGER :: admm_maxiter = 1000
  ! Initial ADMM penalty
  REAL(DP) :: admm_penalty = 1.0D0
  ! Convergence criteria of ADMM algorithm
  REAL(DP) :: admm_tol = 1.0D-3
  ! Boundaries and step-length for searching optimal lambda
  REAL(DP) :: loglam_min = -6.0D0
  REAL(DP) :: loglam_max = 2.0D0
  REAL(DP) :: loglam_step = 1.0D0
  !
  ! Local Variables
  !
  ! Input Green's function on the imaginary-time axis
  REAL(DP), ALLOCATABLE :: Gtau(:)
  ! Output spectral function on the real-frequency axis
  REAL(DP), ALLOCATABLE :: Aw_opt(:)
  ! Rows of SVD kernel function (ntau+1)
  INTEGER :: svd_M = 1
  ! Columns of SVD kernel function (nw+1)
  INTEGER :: svd_N = 1
  ! Reduced SVD rank
  INTEGER :: svd_L = 1
  ! SVD U matrix (reduced)
  REAL(DP), ALLOCATABLE :: svd_U(:, :)
  ! SVD S matrix (reduced)
  REAL(DP), ALLOCATABLE :: svd_S(:, :)
  ! SVD V^T matrix (reduced)
  REAL(DP), ALLOCATABLE :: svd_VT(:, :)
  ! Number of lambda to find the optimal spectra
  INTEGER :: nlam = 1
  ! Square errors for different lambda
  REAL(DP), ALLOCATABLE :: chi2(:)
  ! Spectra obtained on different lambda
  REAL(DP), ALLOCATABLE :: Aw(:, :)

CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE init_spm
    !
    IMPLICIT NONE

    ! Namelist to read from input file
    NAMELIST /SPM_PARAMS/ statistics, ntau, beta, nw, wmin, wmax, svd_cut,  &
      admm_auto_penalty, admm_maxiter, admm_penalty, admm_tol, loglam_min,  &
      loglam_max, loglam_step

    ! Open input file
    OPEN(UNIT=11, FILE='params.in', STATUS='OLD')
    READ(11, NML=SPM_PARAMS)
    CLOSE(11)

    ! Initialize local variables
    svd_M = ntau + 1
    svd_N = nw + 1
    nlam = MAX(FLOOR((loglam_max - loglam_min) / loglam_step), 0) + 1

    ! Initialize arrays
    ALLOCATE(Gtau(ntau+1), Aw_opt(nw+1))
    ALLOCATE(chi2(nlam), Aw(nw+1, nlam))
    Gtau = 0.0D0 ; Aw_opt = 0.0D0
    chi2 = 0.0D0 ; Aw = 0.0D0

  END SUBROUTINE init_spm
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE fin_spm
    !
    IMPLICIT NONE

    DEALLOCATE(Gtau, Aw_opt)
    DEALLOCATE(chi2, Aw)

  END SUBROUTINE fin_spm
  !-----------------------------------------------------------------------------

END MODULE spm_util

