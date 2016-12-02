C       IMPLICIT NONE
      subroutine iterate(fg,fgv,f,rho,f_new,rho_new,n,mask,it)
C 
C === begin external variable declarations, including f2py statements ===
      integer :: n 
Cf2py integer intent(hide),depend(f) :: n=len(f)
      integer it
Cf2py integer intent(in) :: it
      double precision fg(n,n)
Cf2py double precision intent(in) :: fg
      double precision fgv(n,n), mask(n,n)
Cf2py double precision intent(in) :: fgv, mask
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
      integer :: ix, ig, iu
      double precision fun1(n), fun2(n), fun3(n)
      double precision nom(n,n), denom(n,n)
      double precision up, down, var
C === end internal variable declarations ===
C
C ====== start main program ======
      n = size(f)
C     zeroing everything that can be zeroed:
      do i = 1, n
        fun1(i) = 0.
        fun2(i) = 0.
        fun3(i) = 0.
        f_new(i) = 0.
        rho_new(i) = 0.
        do j = 1, n
          nom(i,j) = 0.
          denom(i,j) = 0.
        enddo
      enddo
C     calculate var, the allowed change in f&rho, dep on iteration
C     (assuming it starts at 0)
      if (it.le.4) then
        var = 1.2
      else if (it.le.11) then
        var = 1.1
      else if (it.le.20) then
        var = 1.05
      else if (it.le.29) then
        var = 1.025
      else
        var = 1.01
      endif
C ==  calculate helping functions: ==
      do ix = 1, n
        fun1(ix) = 0.
        fun2(ix) = 0.
        fun3(ix) = 0.
      	do ig = 1, ix-1 
C         Should it be ix-1 to avoid ig=0?
          iu = ix - ig
C         TODO Check whether this masking is right:
          fun1(ix) = fun1(ix) + f(ig)*rho(iu)
C      &        *mask(ix,ig) 
          if (fgv(ix,ig).gt.0.) then
           fun2(ix) = fun2(ix) + (f(ig)*rho(iu)/fgv(ix,ig))**2. 
C      &        *mask(ix,ig)
            fun3(ix) = fun3(ix) + f(ig)*rho(iu)*fg(ix,ig)
     &        /(fgv(ix,ig)**2.)
C      &        *mask(ix,ig)
          endif
        enddo 
C       end ig loop
        if (fun1(ix).gt.0.) then
          fun2(ix) = fun2(ix)/(fun1(ix)**3.)
          fun3(ix) = fun3(ix)/(fun1(ix)**2.)
        else
          fun2(ix) = 0.
          fun3(ix) = 0.
        endif
        do ig = 1, ix
          if ((fun1(ix)*fgv(ix,ig)).gt.0.) then
            nom(ix,ig) = fun2(ix) - fun3(ix) + 
     &        fg(ix,ig)/(fun1(ix)*(fgv(ix,ig)**2.))
            denom(ix,ig) = 1./(fun1(ix)*fgv(ix,ig))**2.
          else
            nom(ix,ig) = 0.
            denom(ix,ig) = 0.
          endif
        enddo
      enddo
C ==  calculate updated f: ==
      do ig = 1, n 
C     NB original uses igmin,igmax for limits above--is it same?
C       Nominator and denominator of current f(ig) point, to sum:
        up = 0.
        down = 0.
        do ix = ig+1, n 
C       Should it be ig+1 to not get iu=0?
          iu = ix - ig
          up = up + rho(iu)*nom(ix,ig)*mask(ix,ig) 
          down = down + (rho(iu)**2.)*denom(ix,ig)
C      &     *mask(ix,ig)
        enddo
        if (down.gt.0) then
          if ((up/down).gt.(var*f(ig))) then
            f_new(ig) = var*f(ig)
          else if ((up/down).lt.(f(ig)/var)) then
            f_new(ig) = f(ig)/var
          else
            f_new(ig) = up/down
          endif
        else
          f_new(ig) = 0.
        endif
      enddo 
C     end ig loop
C ==  calculate updated rho: ==
      do iu = 0, n
C       Nominator and denominator of current f(ig) point, to sum:
        up = 0.
        down = 0.
        do ix = iu, n
          ig = ix - iu
          up = up + f(ig)*nom(ix,ig)*mask(ix,ig)
          down = down + (f(ig)**2.)*denom(ix,ig)
C     &       *mask(ix,ig)
        enddo
        if (down.gt.0) then 
          if ((up/down).gt.var*rho(iu)) then
            rho_new(iu) = var*rho(iu)
          else if ((up/down).lt.rho(iu)/var) then
            rho_new(iu) = rho(iu)/var
          else
            rho_new(iu) = up/down
            if (iu.eq.5) then
              write(6,*)'rho_new(5) = ',rho_new(iu)
            endif
          endif
        else
          rho_new(iu) = 0.
        endif
      enddo
C     end iu loop

C       Test matrix ordering: Turns out to be (row,col).
C       do i = 1, n
C         fg(1,i) = 0
C       enddo
C
C       write(6,*) 'Hello, World!'
C       do i = 1, n
C C         write(*,"(10(F2.4,1X))") fg(i,1:10)
C           write(*,'(A,D5.3)') 'Hello ', f(i)
C C   50      format(A,F1.0)
C       enddo
      end subroutine