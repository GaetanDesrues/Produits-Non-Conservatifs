module Solveur
  implicit none

contains


  subroutine ConditionInitiale(u, dx, type)
    implicit none
    real(kind=8), intent(inout), dimension(:,:) :: u
    real(kind=8), intent(in) :: dx
    integer, intent(in) :: type
    integer :: i
    real(kind=8) :: x, sig, mu
    integer, dimension(2) :: shapeArray

    shapeArray = shape(u)
    u = 0

    select case(type)
    case(0) ! Marche
      do i = 1, shapeArray(1)
        x = i*dx
        if (x<0.5) then
          u(i,1) = 2.1!1
          u(i,2) = 2.3
          u(i,3) = 0.84!3.485!0.4
        else
          u(i,1) = 8
          u(i,2) = 0
          u(i,3) = 0.32
        endif
      enddo

    case(1) ! Gaussienne
      sig = 0.08
      mu = 0.5
      do i = 1, shapeArray(1)
        x = i*dx
        u(i,1) = 1/(sig*sqrt(2*3.1415))*exp(-((x-mu)/sig)**2/2)/100
      enddo

    case(2) ! Sinus
      do i = 1, shapeArray(1)
        x = i*dx
        if ((x>0.4) .and. x<(0.6)) then
          u(i,1) = sin((x-0.4)*3.1415*5)+0.5
        else
          u(i,1) = 0.5
        endif
      enddo

    case default
      write(6,*) "Attention à la condition initiale"

    end select
  end subroutine




  subroutine Iteration(u, sigma, flux)
    implicit none
    real(kind=8), dimension(:,:), intent(inout) :: u
    real(kind=8), intent(in) :: sigma
    integer, intent(in) :: flux
    real(kind=8), dimension(:,:), allocatable :: ubis, Fu, dp, dm
    real(kind=8) :: g, bip, bim, l1, l2, l3, l4, l5, l6, v1, v2, sigSpeed
    integer :: i, j
    integer, dimension(2) :: shapeArray
    real(kind=8), dimension(:), allocatable :: fluxp, fluxm
    character(len=10) :: xx, xxx
    real(kind=8), dimension(3) :: d

    g=9.81

    shapeArray = shape(u)
    allocate(ubis(1:shapeArray(1), 1:shapeArray(2)))
    allocate(Fu(1:shapeArray(1), 1:shapeArray(2)))
    allocate(fluxm(1:shapeArray(2)), fluxp(1:shapeArray(2)))

    ! Système de Saint Venant - flux physique
    ! do i=1, shapeArray(1)
    !   Fu(i,1) = u(i,2)
    !   Fu(i,2) = u(i,2)*u(i,2)/u(i,1)+0.5*g*u(i,1)*u(i,1)
    ! enddo

    ! Système d'Euler coordonnées Lagrangiennes (conservatif)
    ! do i=1, shapeArray(1)
    !   Fu(i,1) = u(i,2)
    !   Fu(i,2) = (u(i,3)-0.5*u(i,2)**2)/(0.4*u(i,1))
    !   Fu(i,3) = u(i,2)*(u(i,3)-0.5*u(i,2)**2)/(0.4*u(i,1))
    ! enddo


    ubis = u

    select case (flux)
    case (0) ! Flux de Lax-Friedrichs
      open(unit=22, file="Output/LF.txt", status="unknown")
      do i=2, shapeArray(1)-1
        !!
        call racines(u(i-1,:), v1, v2)
        l1 = abs(v1)
        call racines(u(i,:), v1, v2)
        l3 = abs(v1)
        bim = max(l1, l3)
        call racines(u(i+1,:), v1, v2)
        l5 = abs(v1)
        bip = max(l3, l5)
        !!
        write(xx,'(F10.6)') bim
        write(xxx,'(F10.6)') bip
        write(22,*) xx // "  " // xxx
        fluxm = 0.5*(Fu(i,:)+Fu(i-1,:))-0.5/sigma*(ubis(i,:)-ubis(i-1,:))
        fluxp = 0.5*(Fu(i+1,:)+Fu(i,:))-0.5/sigma*(ubis(i+1,:)-ubis(i,:))
        u(i,:) = ubis(i,:) - sigma*(fluxp - fluxm)
      enddo
      close(22)


    case (1) ! Flux de Rusanov
      do i=2, shapeArray(1)-1
        call racines(u(i-1,:), v1, v2)
        l1 = abs(v1)
        call racines(u(i,:), v1, v2)
        l3 = abs(v1)
        bim = max(l1, l3)
        call racines(u(i+1,:), v1, v2)
        l5 = abs(v1)
        bip = max(l3, l5)

        bim=0.1
        bip=0.2

        ! if (isNaN(bim)) stop "Erreur : Valeur propre  =  NaN"
        ! write(6,*) bim, bip

        fluxm = 0.5*(Fu(i,:)+Fu(i-1,:))-bim*0.5*(ubis(i,:)-ubis(i-1,:))
        fluxp = 0.5*(Fu(i+1,:)+Fu(i,:))-bip*0.5*(ubis(i+1,:)-ubis(i,:))

        u(i,:) = ubis(i,:) - sigma*(fluxp - fluxm)
      enddo





    case (2) ! Schéma path-conservative
      allocate(dp(1:shapeArray(1)-1,1:3), dm(1:shapeArray(1)-1,1:3))
      dp = 0
      dm = 0

      do i=1, shapeArray(1)-1
        ! D'après la relation (23) de Pares

        ! if ((u(i+1,1)-u(i,1))/=0) write(6,*) -(u(i+1,2)-u(i,2))/(u(i+1,1)-u(i,1))

        d(1) = u(i,2)-u(i+1,2)
        d(2) = 2.5*(u(i+1,3)/u(i+1,1)-u(i,3)/u(i,1)) !u(i,2)-u(i+1,2)
        d(3) = 1.25*(u(i,3)/u(i,1)+u(i+1,3)/u(i+1,1))*(u(i+1,2)-u(i,2))

        do j=1,3 ! Sépare positif/négatif
          sigSpeed = (u(i+1,j)-u(i,j))
          if (sigSpeed/=0) sigSpeed = d(j)/(u(i+1,j)-u(i,j))
          if (sigSpeed>0) then
            dp(i,j) = d(j)
          else
            dm(i,j) = d(j)
          endif
        enddo

        ! dp contient sur chaque maille (i) un vecteur de dimension 3
        ! égal à D^+_{i+1/2}
        ! dm contient sur chaque maille (i) un vecteur de dimension 3
        ! égal à D^-_{i+1/2}
      enddo

      do i=2, shapeArray(1)-1
        u(i,:) = ubis(i,:) - sigma*(dp(i-1,:) + dm(i,:))
      enddo
      deallocate(dp, dm)

    case default
      write(6,*) "Attention au flux choisi"

    end select
    deallocate(ubis, Fu, fluxm, fluxp)
  end subroutine Iteration


  subroutine racines(u,v1,v2)
    implicit none
    real(kind=8), intent(in), dimension(1:3) :: u
    real(kind=8), intent(inout) :: v1, v2
    real(kind=8) :: a,b,c,d,g

    g = 1.4 - 1
    ! if (isNaN(u(1))) stop "Erreur : Valeur propre  =  NaN1"
    ! if (isNaN(u(2))) stop "Erreur : Valeur propre  =  NaN2"
    ! if (isNaN(u(3))) stop "Erreur : Valeur propre  =  NaN3"
    a = g**2*u(1)**2
    b = g*(u(3)-0.5*u(2)**2-2*u(2)*u(1))
    c = 2.5*u(2)**2 - u(3)

    d = b**2-4*a*c
    if (d<=0) write(6,*) "Attention au déterminant ! d = ", d
    v1 = (-b + sqrt(d))/(2*a)
    v2 = (-b - sqrt(d))/(2*a)
  end subroutine racines



  subroutine SaveSol(u, it, dx)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: u
    real(kind=8), intent(in) :: dx
    integer, intent(in) :: it
    integer :: i
    character(len=10) :: hi, qi, temps, x, ei
    integer, dimension(2) :: shapeArray

    shapeArray = shape(u)

    write(temps, '(I6)') it

    open(unit=15, file="Output/Sol/Sol_it=" // trim(adjustl(temps)) // ".txt", status="unknown")

    do i=1, shapeArray(1)
      write(x,'(F10.6)') i*dx
      write(hi, '(F10.6)') 1./u(i,1) ! Densité
      write(qi, '(F10.6)') u(i,2) ! Vitesse
      ! write(ei, '(F10.6)') (u(i,3)-0.5*u(i,2)**2)/(0.4*u(i,1)) ! Pression
      write(ei, '(F10.6)') u(i,3)/(0.4*u(i,1)) ! Pression
      write(15, *) x // "  " // hi // "  " // qi // "  " // ei

      ! Plan u-p
      ! write(hi, '(F10.6)') u(i,2) ! Vitesse
      ! write(ei, '(F10.6)') (u(i,3)-0.5*u(i,2)**2)/(0.4*u(i,1)) ! Pression
      ! write(15, *) hi // "  " // ei
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




  ! subroutine CFLSub(cfl, dt, dx, a)
  !   implicit none
  !
  !   real(kind=8), intent(in) :: cfl, dx, a
  !   real(kind=8), intent(inout) :: dt
  !   real(kind=8) :: cflX
  !
  !   cflX = cfl*dx/a
  !   if (cflX < dt) then
  !     dt = cflX
  !     print*, "Attention CFL non respectée, dt =", dt
  !   endif
  !
  ! end subroutine CFLSub


  !
  ! real(kind=8) function minmod(a,b) result(c)
  !   implicit none
  !   real(kind=8), intent(in) :: a, b
  !   real(kind=8) :: aa, bb
  !
  !   aa = a
  !   bb = b
  !
  !   if (aa*bb <= 0) then
  !     c = 0
  !   else
  !     aa = abs(aa)
  !     bb = abs(bb)
  !     if (aa<bb) then
  !       c = aa
  !     else
  !       c = bb
  !     endif
  !   endif
  !
  ! end function minmod



end module Solveur
