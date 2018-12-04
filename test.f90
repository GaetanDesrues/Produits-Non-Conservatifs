module Solveur
  implicit none

contains


  subroutine ConditionInitiale(u, dx)
    implicit none
    real(kind=8), intent(inout), dimension(:,:) :: u
    real(kind=8), intent(in) :: dx
    integer :: i
    real(kind=8) :: x, sig, mu
    integer, dimension(2) :: shapeArray

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
    shapeArray = shape(u)
    u = 0
    do i = 1, shapeArray(1)
      x = i*dx

      ! Sinus
      ! if ((x>0.4) .and. x<(0.6)) then
      !   u(i,1) = sin((x-0.4)*3.1415*5)+0.5
      ! else
      !   u(i,1) = 0.5
      ! endif

      ! Gaussienne
      ! sig = 0.08
      ! mu = 1
      ! u(i,1) = 1/(sig*sqrt(2*3.1415))*exp(-((x-mu)/sig)**2/2)/100 + 1.2

      ! Marche
      if (x<20) then
        u(i,1) = 1
      else
        u(i,1) = 1.5
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
    real(kind=8), dimension(:,:), intent(inout) :: u
    real(kind=8), intent(in) :: sigma, phi
    real(kind=8), dimension(:,:), allocatable :: ubis, Fu
    real(kind=8) :: g, bip, bim, l1, l2, l3, l4, l5, l6
    integer :: i
    integer, dimension(2) :: shapeArray
    real(kind=8), dimension(2) :: fluxp, fluxm

    g=9.81

    shapeArray = shape(u)
    allocate(ubis(1:shapeArray(1), 1:shapeArray(2)))
    allocate(Fu(1:shapeArray(1), 1:shapeArray(2)))

    ! Système de Saint Venant
    do i=1, shapeArray(1)
      Fu(i,1) = u(i,2)
      Fu(i,2) = u(i,2)*u(i,2)/u(i,1)+0.5*g*u(i,1)*u(i,1)
    enddo

    ubis = u


    ! Flux de Lax-F
    ! do i=2, shapeArray(1)-1
    !   fluxm = 0.5*(Fu(i,:)+Fu(i-1,:))-0.5/sigma*(ubis(i,:)-ubis(i-1,:))
    !   fluxp = 0.5*(Fu(i+1,:)+Fu(i,:))-0.5/sigma*(ubis(i+1,:)-ubis(i,:))
    !   u(i,:) = ubis(i,:) - sigma*phi*(fluxp - fluxm)
    ! enddo


    ! Flux de Rusanov
    do i=2, shapeArray(1)-1
      l1 = abs(ubis(i-1,2)/ubis(i-1,1)+sqrt(g*ubis(i-1,1))) ! valeurs propres
      l2 = abs(ubis(i-1,2)/ubis(i-1,1)-sqrt(g*ubis(i-1,1)))
      l3 = abs(ubis(i,2)/ubis(i,1)+sqrt(g*ubis(i,1)))
      l4 = abs(ubis(i,2)/ubis(i,1)-sqrt(g*ubis(i,1)))
      bim = max(l1, l2, l3, l4)
      l5 = abs(ubis(i+1,2)/ubis(i+1,1)+sqrt(g*ubis(i+1,1)))
      l6 = abs(ubis(i+1,2)/ubis(i+1,1)-sqrt(g*ubis(i+1,1)))
      bip = max(l3, l4, l5, l6)

      fluxm = 0.5*(Fu(i,:)+Fu(i-1,:))-bim*0.5*(ubis(i,:)-ubis(i-1,:))
      fluxp = 0.5*(Fu(i+1,:)+Fu(i,:))-bip*0.5*(ubis(i+1,:)-ubis(i,:))
      u(i,:) = ubis(i,:) - sigma*phi*(fluxp - fluxm)
    enddo

    deallocate(ubis, Fu)
  end subroutine Iteration





  subroutine SaveSol(u, it, dx)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: u
    real(kind=8), intent(in) :: dx
    integer, intent(in) :: it
    integer :: i
    character(len=10) :: hi, qi, temps, x
    integer, dimension(2) :: shapeArray

    shapeArray = shape(u)

    write(temps, '(I6)') it

    open(unit=15, file="Output/Sol_it=" // trim(adjustl(temps)) // ".txt", status="unknown")

    do i=1, shapeArray(1)
      write(x,'(F10.6)') i*dx
      write(hi, '(F10.6)') u(i,1)
      write(qi, '(F10.6)') u(i,2)
      write(15, *) x // "  " // hi // "  " // qi
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
