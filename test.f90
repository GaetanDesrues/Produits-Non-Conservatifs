module Solveur
  implicit none

contains


  subroutine ConditionInitiale(u, dx)
    implicit none
    real(kind=8), intent(inout), dimension(:) :: u
    real(kind=8), intent(in) :: dx
    integer :: i
    real(kind=8) :: x, sig, mu

    ! sig = 0.08
    ! mu = 0.5
    !
    ! do i = 1, size(u)
    !   x = i*dx
    !   if ((x<0.4) .or. x>(0.6)) then
    !     u(i) = 0
    !   else
    !     u(i) = sin((x-0.4)*3.1415*5)
    !   endif
    !   ! u(i) = 1/(sig*sqrt(2*3.1415))*exp(-((x-mu)/sig)**2/2)
    ! enddo


    do i = 1, size(u)
      x = i*dx
      if (x<0.5) then
        u(i) = 1
      else
        u(i) = 0
      endif
    enddo
  end subroutine



  real(kind=8) function minmod(a,b) result(c)
    implicit none
    real(kind=8), intent(in) :: a, b
    real(kind=8) :: aa, bb

    aa = a
    bb = b

    if (aa*bb <= 0) then
      c = 0
    else
      aa = abs(aa)
      bb = abs(bb)
      if (aa<bb) then
        c = aa
      else
        c = bb
      endif
    endif

  end function minmod


  subroutine CFLSub(cfl, dt, dx, a)
    implicit none

    real(kind=8), intent(in) :: cfl, dx, a
    real(kind=8), intent(inout) :: dt
    real(kind=8) :: cflX

    cflX = cfl*dx/a
    if (cflX < dt) then
      dt = cflX
      print*, "Attention CFL non respectée, dt =", dt
    endif

  end subroutine CFLSub




  subroutine Iteration(u, sigma, phi)
    implicit none
    real(kind=8), dimension(:), intent(inout) :: u
    real(kind=8), intent(in) :: sigma, phi
    real(kind=8), dimension(size(u)) :: ubis
    real(kind=8) :: i

    ubis = u
    u(1) = 0
    do i=2, size(u)
      u(i) = ubis(i) - sigma*phi*(ubis(i)-ubis(i-1))
    enddo


  end subroutine Iteration





  subroutine SaveSol(u, t, dx)
    implicit none
    real(kind=8), dimension(:), intent(in) :: u
    real(kind=8), intent(in) :: t, dx
    integer :: i
    character(len=5) :: ui
    character(len=4) :: temps
    character(len=6) :: x

    write(temps, '(F4.2)') t

    open(unit=15, file="Output/Sol_t=" // temps // ".txt", status="unknown")

    do i=1, size(u)
      write(x,'(F6.3)') i*dx
      write(ui, '(F5.2)') u(i)
      write(15, *) x // "  " // ui
    enddo

    close(15)
  end subroutine SaveSol



  subroutine SolExacte(xmax, dx, t, u, a)
    implicit none
    real(kind=8), intent(inout), dimension(:) :: u
    real(kind=8), intent(in) :: dx, xmax, t, a
    real(kind=8) :: x, sig, mu
    integer :: i

    ! sig = 0.08
    ! mu = 0.5

    ! do i = 1, size(u)
    !   x = i*dx - a*t
    !   if ((x<0.4) .or. x>(0.6)) then
    !     u(i) = 0
    !   else
    !     u(i) = sin((x-0.4)*3.1415*5)
    !   endif
    !   ! u(i) = 1/(sig*sqrt(2*3.1415))*exp(-((x-mu)/sig)**2/2)
    ! enddo

    do i = 1, size(u)
      x = i*dx - a*t
      if (x<0.5) then
        u(i) = 1
      else
        u(i) = 0
      endif
    enddo


  end subroutine SolExacte



end module Solveur
