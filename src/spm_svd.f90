MODULE spm_svd
  USE lib_math, ONLY : DP, XSVD
  USE spm_util, ONLY : beta, wmin, wmax, svd_cut
  USE spm_util, ONLY : ntau, nw, svd_M, svd_N, svd_L
  USE spm_util, ONLY : svd_U, svd_S, svd_VT
  IMPLICIT NONE

CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE calc_svd_fermion
    !
    IMPLICIT NONE

    !
    ! Local variables
    !
    ! Rank of S matrix
    INTEGER :: svd_L_orig = 1
    ! Kernel function
    REAL(DP), ALLOCATABLE :: Kmat(:, :)
    ! SVD arrays
    REAL(DP), ALLOCATABLE :: Umat(:, :), Svec(:), VTmat(:, :)
    ! Dummy variables
    INTEGER :: itau, iw, il
    REAL(DP) :: tau, w

    !
    ! Setup the kernel function
    !
    ALLOCATE(Kmat(svd_M, svd_N))

    DO itau = 1, ntau+1
      !
      tau = (itau - 1) * beta / DBLE(ntau)
      DO iw = 1, nw+1
        !
        w = wmin + (iw - 1) * (wmax - wmin) / DBLE(nw)
        ! Fermi kernel
        IF (w >= 0) THEN
          Kmat(itau, iw) = EXP(-tau * w) / (1.0D0 + EXP(-beta * w))
        ELSE
          Kmat(itau, iw) = EXP((beta - tau) * w) / (1.0D0 + EXP(beta * w))
        END IF
        !
      END DO  ! iw
    END DO  ! itau

    !
    ! SVD
    !
    svd_L_orig = MIN(svd_M, svd_N)
    ALLOCATE(Umat(svd_M, svd_M))
    ALLOCATE(Svec(svd_L_orig))
    ALLOCATE(VTmat(svd_N, svd_N))
    !
    CALL XSVD(Kmat, svd_M, svd_N, Umat, Svec, VTmat)

    !
    ! Reduce the dimension of SVD arrays by a cut off
    !
    DO il = 1, svd_L_orig
      IF (ABS(Svec(il)) < svd_cut) THEN
        svd_L = il - 1
        EXIT
      END IF
    END DO  ! il
    ! Allocate reduced SVD matrices
    ALLOCATE(svd_U(svd_M, svd_L))
    ALLOCATE(svd_S(svd_L, svd_L))
    ALLOCATE(svd_VT(svd_L, svd_N))
    ! Reduced arrays
    svd_U = 0.0D0 ; svd_VT = 0.0D0
    svd_U(1:svd_M, 1:svd_L) = Umat(1:svd_M, 1:svd_L)
    svd_VT(1:svd_L, 1:svd_N) = VTmat(1:svd_L, 1:svd_N)
    svd_S = 0.0D0
    DO il = 1, svd_L
      svd_S(il, il) = Svec(il)
    END DO  ! il

    !
    ! Deallocate local arrays
    !
    DEALLOCATE(Kmat, Umat, Svec, VTmat)

  END SUBROUTINE calc_svd_fermion
  !-----------------------------------------------------------------------------

END MODULE spm_svd
