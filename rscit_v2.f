      subroutine iterate(fg,fgv,f,rho,f_new,rho_new,n,mask)
C       IMPLICIT NONE
C 
C === begin external variable declarations, including f2py statements ===
      integer :: n, ix, ig 
Cf2py integer intent(hide),depend(f) :: n=len(f)
      double precision fg(n,n), fgv(n,n), mask(n,n)
Cf2py double precision intent(in) :: fg, fgv, mask
      double precision f(n)
Cf2py double precision intent(in) :: f
      double precision rho(n)
Cf2py double precision intent(in) :: rho
      double precision f_new(n)
Cf2py double precision intent(out) :: f_new
      double precision rho_new(n)
Cf2py double precision intent(out) :: rho_new
C === end external variable declarations ===
C === begin internal variable declarations ===
      double precision s(n), a(n), b(n)
C === end internal variable declarations ===
C
C === start main program ===
      n = size(f)
      do ix = 1, n
        s(ix) = 0
        a(ix) = 0
        b(ix) = 0
      	do ig = 1, n
          
          s(ix) = s(ix) + ! remember to include mask
        enddo
      enddo
C
C       write(6,*) 'Hello, World!'
C       do i = 1, n
C C         write(*,"(10(F2.4,1X))") fg(i,1:10)
C           write(*,'(A,D5.3)') 'Hello ', f(i)
C C   50      format(A,F1.0)
C       enddo
      end subroutine