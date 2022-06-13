MODULE spm_admm
  !
  ! According to J. Otsuki et al. PRE 95, 061302(R) (2017),
  ! which is denoted as [Otsuki2017]. The corresponding C++
  ! implementation is given in https://github.com/SpM-lab/SpM .
  !
  USE lib_math, ONLY : DP, MINV, NORM_L2
  USE spm_util, ONLY : Gtau, Aw, chi2
  USE spm_util, ONLY : svd_M, svd_N, svd_L, svd_U, svd_S, svd_VT
  USE spm_util, ONLY : admm_auto_penalty, admm_maxiter, admm_tol
  USE spm_util, ONLY : admm_penalty
  USE spm_util, ONLY : loglam_min, loglam_max, loglam_step
  IMPLICIT NONE

  ! ADMM LOG FILE
  INTEGER, PARAMETER :: ULOG = 21
  !
  ! Hyper-parameter
  REAL(DP) :: lambda
  !
  ! Penalty factors: mu' and mu in [Otsuki2017], correspondingly
  REAL(DP) :: mu1, mu2
  !
  ! Pre-stored factors:
  REAL(DP) :: Pfac
  ! (S^T S / lambda + (mu' + mu) I)^{-1} in [Otsuki2017]
  REAL(DP), ALLOCATABLE :: Dinv(:, :)
  ! e vector in [Otsuki2017]
  REAL(DP), ALLOCATABLE :: evec(:)
  ! V^T e in [Otsuki2017]
  REAL(DP), ALLOCATABLE :: Vrow(:)
  ! S^T / lambda in [Otsuki2017]
  REAL(DP), ALLOCATABLE :: Sl(:, :)
  ! D^{-1} Sl, D^{-1} mu', D^{-1} mu V^T & D^{-1} V^T e
  REAL(DP), ALLOCATABLE :: Pmat1(:, :), Pmat2(:, :), Pmat3(:, :), Pmat4(:)
  !
  ! If the minimization is converged
  LOGICAL :: cflag
  ! Residual errors
  REAL(DP) :: res1_prim, res1_dual, res2_prim, res2_dual
  ! x' (=rho') and y' (=Gtau') in [Otsuki2017]
  REAL(DP), ALLOCATABLE :: x_(:), y_(:)

CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE admm_minimizer(ilambda)
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ilambda
    !
    ! Local variables
    !
    ! If penalty has been update
    LOGICAL :: uflag
    ! Loop variables
    INTEGER :: iter, uskip, ucount
    ! Number of particles
    REAL(DP) :: nParticle
    ! Dummy factors
    REAL(DP) :: nu
    ! Auxiliary vectors: z', u', z & u in [Otsuki2017], correspondingly
    REAL(DP), ALLOCATABLE :: z1(:), u1(:), z2(:), u2(:)
    REAL(DP), ALLOCATABLE :: z1_old(:), z2_old(:)
    ! Vmat = VT^T, Vx = V x'
    REAL(DP), ALLOCATABLE :: Vmat(:, :), Vx(:)

    !
    ! Initialization
    !
    lambda = 10.0D0 ** (loglam_min+(ilambda-1)*loglam_step)
    nParticle = Gtau(1) + Gtau(svd_M)
    mu1 = admm_penalty ; mu2 = admm_penalty
    uskip = 20 ; ucount = 0
    !
    ALLOCATE(Dinv(svd_L, svd_L), evec(svd_N), Vrow(svd_L), Sl(svd_L, svd_L))
    ALLOCATE(Pmat1(svd_L, svd_L), Pmat2(svd_L, svd_L))
    ALLOCATE(Pmat3(svd_L, svd_N), Pmat4(svd_L))
    ALLOCATE(z1(svd_L), u1(svd_L), z2(svd_N), u2(svd_N))
    ALLOCATE(z1_old(svd_L), z2_old(svd_N))
    ALLOCATE(x_(svd_L), y_(svd_L), Vmat(svd_N, svd_L), Vx(svd_N))
    !
    CALL init_prefact
    !
    x_ = 0.0D0 ; z1 = 0.0D0 ; u1 = 0.0D0
    z2 = 0.0D0 ; u2 = 0.0D0
    y_ = MATMUL(TRANSPOSE(svd_U), Gtau)
    Vmat = TRANSPOSE(svd_VT)

    !
    ! Minimization
    !
    DO iter = 1, admm_maxiter

      !
      ! Record the previous results
      !
      z1_old = z1
      z2_old = z2

      !
      ! Update for minimization
      !
      ! x' <- Dinv * (Sl y' + mu' (z'-u') + mu VT (z-u) + nu VT e)
      x_ = MATMUL(Pmat1, y_)
      x_ = x_ + MATMUL(Pmat2, z1-u1)
      x_ = x_ + MATMUL(Pmat3, z2-u2)
      ! nu = (1 - < V xi1 >) / < V xi2 >
      nu = (nParticle - DOT_PRODUCT(Vrow, x_)) / Pfac
      ! x' = xi1 + nu xi2
      x_ = x_ + nu * Pmat4

      ! z' <- S_{1/mu'}(x' + u')
      CALL soft_thresholding(x_+u1, svd_L, 1.0D0/mu1, z1)
      !
      ! u' <- u' + x' - z'
      u1 = u1 + x_ - z1
      !
      ! z <- P+(V x' + u)
      Vx = MATMUL(Vmat, x_)
      CALL positive(Vx+u2, svd_N, z2)
      !
      ! u <- u + V x' - z
      u2 = u2 + Vx - z2

      !
      ! Check convergence
      !
      ! Residual errors
      res1_prim = NORM_L2(x_-z1, svd_L)
      res1_dual = NORM_L2(z1-z1_old, svd_L) * mu1
      res2_prim = NORM_L2(Vx-z2, svd_N)
      res2_dual = NORM_L2(z2-z2_old, svd_N) * mu2
      !
      ! Criteria: absolute tolerance (primary residual error per element)
      cflag = (res1_prim / SQRT(1.D0*svd_L)) < admm_tol
      cflag = cflag .AND. ((res2_prim / SQRT(1.D0*svd_N)) < admm_tol)
      IF (cflag) EXIT
      
      !
      ! Update pre-factors if auto_penalty
      !
      IF (admm_auto_penalty) THEN
        !
        ucount = ucount + 1
        !
        IF (MOD(ucount, uskip) == 0) THEN
          !
          uflag = .FALSE.
          ! mu'
          uflag = uflag .OR. update_penalty(res1_prim,res1_dual,mu1,u1,svd_L)
          ! mu
          uflag = uflag .OR. update_penalty(res2_prim,res2_dual,mu2,u2,svd_N)
          !
          ! Update the pre-factors
          IF (uflag) THEN
            CALL init_prefact
            ucount = 0
            uskip = uskip * 2
          END IF
          !
        END IF
      END IF

    END DO  ! iter

    !
    ! Record chi2, Aw and write log file
    !
    ! rho = V x'
    Aw(1:svd_N, ilambda) = z2
    !
    ! chi2 = 1/2 \sum_l (G'_l - s_l rho'_l)^2
    y_(1:svd_L) = y_(1:svd_L) - MATMUL(svd_S, x_)
    chi2(ilambda) = 0.5D0 * NORM_L2(y_, svd_L)**2
    !
    CALL write_log(ilambda)

    !
    ! Deallocation
    !
    DEALLOCATE(Dinv, evec, Vrow, Sl)
    DEALLOCATE(Pmat1, Pmat2, Pmat3, Pmat4)
    DEALLOCATE(z1, u1, z2, u2, z1_old, z2_old)
    DEALLOCATE(x_, y_, Vmat, Vx)

  END SUBROUTINE admm_minimizer
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE init_prefact
    IMPLICIT NONE
    !
    INTEGER :: i, L

    L = svd_L
    !
    ! Dinv = (S^T S / lambda + (mu' + mu) I)^{-1}
    Dinv = 0.0D0
    DO i = 1, L
      Dinv(i, i) = svd_S(i, i)**2 / lambda + (mu1 + mu2)
    END DO  ! i
    CALL MINV(Dinv, L)
    !
    ! e_j = 1
    evec = 1.0D0
    !
    ! Vrow = e^T V
    Vrow = MATMUL(evec, TRANSPOSE(svd_VT))
    !
    ! Sl = S^T / lambda
    Sl = svd_S / lambda
    !
    ! Pmat1 = Dinv * Sl
    Pmat1 = MATMUL(Dinv, Sl)
    !
    ! Pmat2 = Dinv mu'
    Pmat2 = Dinv * mu1
    !
    ! Pmat3 = Dinv * mu * V^T
    Pmat3 = MATMUL(Dinv, svd_VT) * mu2
    !
    ! Pmat4 = Dinv * V^T e
    Pmat4 = MATMUL(Dinv, MATMUL(svd_VT, evec))
    !
    ! Pfac = Vrow * Pmat4 = < V xi2 >
    Pfac = DOT_PRODUCT(Vrow, Pmat4)

  END SUBROUTINE init_prefact
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE soft_thresholding(vin, nv, alpha, vout)
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nv
    REAL(DP), INTENT(IN) :: alpha, vin(nv)
    REAL(DP), INTENT(OUT) :: vout(nv)
    !
    INTEGER :: iv
    REAL(DP) :: x

    DO iv = 1, nv
      x = vin(iv)
      IF (x >= -alpha .AND. x <= alpha) THEN
        vout(iv) = 0.0D0
      ELSE IF (x > alpha) THEN
        vout(iv) = x - alpha
      ELSE
        vout(iv) = x + alpha
      END IF
    END DO  ! iv

  END SUBROUTINE soft_thresholding
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE positive(vin, nv, vout)
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nv
    REAL(DP), INTENT(IN) :: vin(nv)
    REAL(DP), INTENT(OUT) :: vout(nv)
    !
    INTEGER :: iv
    REAL(DP) :: x

    DO iv = 1, nv
      x = vin(iv)
      IF (x > 0.0D0) THEN
        vout(iv) = x
      ELSE
        vout(iv) = 0.0D0
      END IF
    END DO  ! iv

  END SUBROUTINE positive
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  FUNCTION update_penalty(rp, rd, penalty, uvec, n)
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n
    ! primary residual and dual residual
    REAL(DP), INTENT(IN) :: rp, rd
    REAL(DP), INTENT(INOUT) :: penalty, uvec(n)
    LOGICAL :: update_penalty
    !
    REAL(DP), PARAMETER :: MAX_RESIDUAL_RATIO = 10.0D0
    REAL(DP), PARAMETER :: PENALTY_INCREMENT = 2.0D0
    REAL(DP) :: fac = PENALTY_INCREMENT

    !
    update_penalty = .TRUE.
    !
    ! Update penalties according to Eq.(3.13) in Boyd et al., Foundations and
    ! Trends in Machine Learning Vol. 3, No. 1 (2010) 1–122
    IF (rp > rd * MAX_RESIDUAL_RATIO) THEN
      fac = fac
    ELSE IF (rd > rp * MAX_RESIDUAL_RATIO) THEN
      fac = 1.0D0 / fac
    ELSE
      update_penalty = .FALSE.
      RETURN
    END IF

    penalty = penalty * fac
    uvec = uvec / fac
    RETURN

  END FUNCTION update_penalty
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE write_log(ilam)
    USE lib_math, ONLY : NORM_L1
    IMPLICIT NONE

    !
    INTEGER, INTENT(IN) :: ilam

    ! Open ADMM log file
    IF (ilam == 1) THEN
      OPEN(UNIT=ULOG, FILE='admm.out', STATUS='UNKNOWN')
    ELSE
      OPEN(UNIT=ULOG, FILE='admm.out', ACCESS='APPEND')
    END IF

    WRITE(ULOG, "(A)") "==============================================="
    WRITE(ULOG, "(X,A,E10.3)") "lambda = ", lambda
    IF (cflag) THEN
      WRITE(ULOG, "(X,A)") "converged"
    ELSE
      WRITE(ULOG, "(X,A)") "not converged"
    END IF
    WRITE(ULOG, "(X,A)") "|"
    WRITE(ULOG, "(X,A)") "|-- Residual Errors"
    WRITE(ULOG, "(X,A,'=',E10.3)") "| L1 Regular. (prim) ", res1_prim
    WRITE(ULOG, "(X,A,'=',E10.3)") "|             (dual) ", res1_dual
    WRITE(ULOG, "(X,A,'=',E10.3)") "| x >= 0      (prim) ", res2_prim
    WRITE(ULOG, "(X,A,'=',E10.3)") "|             (dual) ", res2_dual
    WRITE(ULOG, "(X,A)") "|"
    WRITE(ULOG, "(X,A)") "|-- Target Quantities"
    WRITE(ULOG, "(X,A,'=',E10.3)") "| || y'-Sx' ||_2     ", NORM_L2(y_, svd_L)
    WRITE(ULOG, "(X,A,'=',E10.3)") "| ||   x'   ||_1     ", NORM_L1(x_, svd_L)
    WRITE(ULOG, "(X,A)") "|"
    WRITE(ULOG, "(X,A)") "|-- Penalty Parameters"
    WRITE(ULOG, "(X,A,'=',E10.3)") "| for L1 Regular.    ", mu1
    WRITE(ULOG, "(X,A,'=',E10.3)") "| for x >= 0         ", mu2
    WRITE(ULOG, "(X,A)") "|"

    CLOSE(ULOG)

  END SUBROUTINE write_log
  !-----------------------------------------------------------------------------

END MODULE spm_admm
