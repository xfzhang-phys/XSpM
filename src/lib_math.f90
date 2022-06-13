MODULE lib_math
  USE, INTRINSIC :: iso_fortran_env
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = REAL64
  REAL(DP), PARAMETER :: PI = ACOS(-1.0D0)

CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE MINV(mat, n)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(INOUT) :: mat(n, n)
    INTEGER :: ipiv(n), info
    REAL(DP) :: work(n)

    CALL DGETRF(n, n, mat, n, ipiv, info)
    CALL DGETRI(n, mat, n, ipiv, work, n, info)

  END SUBROUTINE MINV
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE XSVD(a, m, n, u, s, vt)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: m, n
    REAL(DP), INTENT(INOUT) :: a(m, n)
    REAL(DP), INTENT(OUT) :: u(m, m), s(MIN(m,n)), vt(n, n)
    ! local variables
    INTEGER :: lda, ldu, ldvt, lwork, info
    REAL(DP), ALLOCATABLE :: work(:)

    lda = m ; ldu = m ; ldvt = n
    lwork = 4*MAX(3*MIN(m,n)+MAX(m,n), 5*MIN(m,n))

    ALLOCATE(work(lwork))
    !
    CALL DGESVD('A', 'A',  m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
    !
    DEALLOCATE(work)

    IF (info /= 0) THEN
      PRINT*, "FAILED TO PERFORM SVD!"
      STOP
    END IF

  END SUBROUTINE XSVD
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  FUNCTION NORM_L1(vec, n)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: vec(n)
    REAL(DP) :: NORM_L1

    NORM_L1 = SUM(ABS(vec))

  END FUNCTION NORM_L1
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  FUNCTION NORM_L2(vec, n)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: vec(n)
    REAL(DP) :: NORM_L2

    NORM_L2 = SQRT(SUM(vec**2))

  END FUNCTION NORM_L2
  !-----------------------------------------------------------------------------

END MODULE lib_math
