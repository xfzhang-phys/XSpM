PROGRAM SpM
  USE lib_math, ONLY : DP
  USE spm_util, ONLY : ntau, nw, wmin, wmax, Gtau, Aw_opt
  USE spm_util, ONLY : init_spm, fin_spm
  IMPLICIT NONE

  !
  ! Local variables
  !
  INTEGER :: itau, iw
  REAL(DP) :: t_, g_, w_, dw_
  REAL(DP), ALLOCATABLE :: Gfs(:), Aws(:)

  !
  ! Initialize SpM
  !
  CALL init_spm
  !
  ALLOCATE(Gfs(ntau+1), Aws(nw+1))

  !
  ! Read input Green's function on the imaginary-time axis
  !
  OPEN(UNIT=11, FILE='Gt.in', STATUS='OLD')
  DO itau = 1, ntau + 1
    READ(11, *) t_, g_
    Gfs(itau) = -g_
  END DO  ! itau
  CLOSE(11)
  !
  Gtau(1:ntau+1) = Gfs(1:ntau+1)

  !
  ! Perform SpM analytical continuation
  !
  CALL run_spm
  !
  Aws(1:nw+1) = Aw_opt(1:nw+1)
  
  !
  ! Write the obtained spectra on the real-frequency axis
  !
  dw_ = (wmax - wmin) / DBLE(nw)
  OPEN(UNIT=11, FILE='Aw.dat', STATUS='UNKNOWN')
  DO iw = 1, nw + 1
    w_ = wmin + (iw - 1) * dw_
    WRITE(11, '(F12.6, F20.13)') w_, Aws(iw) / dw_
  END DO  ! iw
  CLOSE(11)

  !
  ! Finalize SpM
  !
  CALL fin_spm
  !
  DEALLOCATE(Gfs, Aws)

CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE run_spm
    USE spm_util, ONLY : statistics, Aw, chi2
    USE spm_util, ONLY : loglam_min, loglam_max, loglam_step, nlam
    USE spm_svd,  ONLY : calc_svd_fermion
    USE spm_admm, ONLY : admm_minimizer
    IMPLICIT NONE

    !
    INTEGER :: ilambda, ilam_opt
    REAL(DP) :: loga, b, loglam, max_, f_x2(nlam)

    !
    ! Initialize SVD
    !
    IF (statistics == 'F' .OR. statistics == 'f') THEN
      CALL calc_svd_fermion
    END IF

    !
    ! Solve the LASSO problem with different lambda by ADMM algorithm
    !
    DO ilambda = 1, nlam
      CALL admm_minimizer(ilambda)
    END DO  ! ilambda

    !
    ! Choose the best spectral function with optimal lambda
    !
    ! f = a lambda^b ==> log(f) = b log(lambda) + log(a)
    b = (LOG10(chi2(nlam)) - LOG10(chi2(1))) / (loglam_max - loglam_min)
    loga = LOG10(chi2(nlam)) - b * loglam_max
    !
    ! Find the optimal lambda
    OPEN(UNIT=11, FILE='find_optimal_lambda.out', STATUS='UNKNOWN')
    WRITE(11, "(A)") "#   log(lambda)   f(lambda)/chi2(lambda)    chi2(lambda)"
    DO ilambda = 1, nlam
      !
      loglam = loglam_min + (ilambda - 1) * loglam_step
      f_x2(ilambda) = b * loglam + loga - LOG10(chi2(ilambda))
      !
      WRITE(11, "(3F13.8)") loglam, f_x2(ilambda), chi2(ilambda)
    END DO  ! ilambda
    !
    max_ = -1.0D0
    DO ilambda = 1, nlam
      IF (f_x2(ilambda) > max_) THEN
        ilam_opt = ilambda
        max_ = f_x2(ilambda)
      END IF
    END DO
    WRITE(11, *)
    WRITE(11, "(A, F6.2)") "* lam_opt = 1.0E^", loglam_min+(ilam_opt-1)*loglam_step
    CLOSE(11)

    !
    ! Set the 'best' spectral function
    !
    Aw_opt(1:nw+1) = Aw(1:nw+1, ilam_opt)

  END SUBROUTINE run_spm
  !-----------------------------------------------------------------------------

END PROGRAM SpM
