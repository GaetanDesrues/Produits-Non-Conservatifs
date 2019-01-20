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
          u(i,1) = 2.1
          u(i,2) = 2.3
          u(i,3) = 0.84 ! non-conservative
          ! u(i,3) = 3.485 ! conservative
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

    real(kind=8), dimension(:,:), allocatable :: ubis, Fu, dp, dm, dp2, dm2
    real(kind=8), dimension(:), allocatable :: fluxp, fluxm, fluxp2, fluxm2

    real(kind=8) :: bim, bi, bip, sigSpeed, pm, p, pp, thetap, thetam, limm, limp
    integer :: i, j
    integer, dimension(2) :: shapeArray
    real(kind=8), dimension(3) :: d, d2

    shapeArray = shape(u)
    allocate(ubis(1:shapeArray(1), 1:shapeArray(2)))
    allocate(Fu(1:shapeArray(1), 1:shapeArray(2)))
    allocate(fluxm(1:shapeArray(2)), fluxp(1:shapeArray(2)))
    allocate(fluxm2(1:shapeArray(2)), fluxp2(1:shapeArray(2)))

    ! Système de Saint Venant - flux physique
    ! do i=1, shapeArray(1)
    !   Fu(i,1) = u(i,2)
    !   Fu(i,2) = u(i,2)*u(i,2)/u(i,1)+0.5*9.81*u(i,1)*u(i,1)
    ! enddo

    ! Système d'Euler coordonnées Lagrangiennes (conservatif)
    ! do i=1, shapeArray(1)
    !   p = (u(i,3)-0.5*u(i,2)**2)/(0.4*u(i,1))
    !   Fu(i,1) = -u(i,2)
    !   Fu(i,2) = p
    !   Fu(i,3) = p*u(i,2)
    ! enddo


    ubis = u

    select case (flux)
    case (0) ! Flux de Lax-Friedrichs
      do i=2, shapeArray(1)-1
        fluxm = 0.5*(Fu(i,:)+Fu(i-1,:))-0.5/sigma*(ubis(i,:)-ubis(i-1,:))
        fluxp = 0.5*(Fu(i+1,:)+Fu(i,:))-0.5/sigma*(ubis(i+1,:)-ubis(i,:))
        u(i,:) = ubis(i,:) - sigma*(fluxp - fluxm)
      enddo


    case (1) ! Flux de Rusanov
      do i=3, shapeArray(1)-2
        pm = (u(i-1,3)-0.5*u(i-1,2)**2)/(0.4*u(i-1,1))
        p = (u(i,3)-0.5*u(i,2)**2)/(0.4*u(i,1))
        pp = (u(i+1,3)-0.5*u(i+1,2)**2)/(0.4*u(i+1,1))

        bim = sqrt(3.5*pm/u(i-1,1))
        bi = sqrt(3.5*p/u(i,1))
        bip = sqrt(3.5*pp/u(i+1,1))

        bim = max(bim, bi)
        bip = max(bi, bip)

        fluxm = 0.5*(Fu(i,:)+Fu(i-1,:))-bim*0.5*(ubis(i,:)-ubis(i-1,:))
        fluxp = 0.5*(Fu(i+1,:)+Fu(i,:))-bip*0.5*(ubis(i+1,:)-ubis(i,:))

        ! ordre 2 et limiteur de pente
        do j=1,3
          thetam = u(i,j)-u(i-1,j)
          if (thetam/=0) thetam = (u(i-1,j)-u(i-2,j))/thetam
          thetap = u(i+1,j)-u(i,j)
          if (thetap/=0) thetap = (u(i,j)-u(i-1,j))/thetap

          fluxm2(j) = 0.25*(2*Fu(i,j)+Fu(i+1,j)+Fu(i-1,j))-bim*0.25*(ubis(i+1,j)-ubis(i-1,j))
          fluxp2(j) = 0.25*(2*Fu(i+1,j)+Fu(i+2,j)+Fu(i,j))-bip*0.25*(ubis(i+2,j)-ubis(i,j))

          limm = 0 ! ordre 1
          limp = 0 ! ordre 1
          ! limm = minmod(1.d0,thetam)
          ! limp = minmod(1.d0,thetap)
          ! limm = superbee(thetam)
          ! limp = superbee(thetap)
          ! limm = vonLeer(thetam)
          ! limp = vonLeer(thetap)

          fluxm(j) = fluxm(j) + limm*(fluxm2(j) - fluxm(j))
          fluxp(j) = fluxp(j) + limp*(fluxp2(j) - fluxp(j))
        enddo

        u(i,:) = ubis(i,:) - sigma*(fluxp - fluxm)
      enddo





    case (2) ! Schéma path-conservative
      allocate(dp(1:shapeArray(1)-1,1:3), dm(1:shapeArray(1)-1,1:3))
      allocate(dp2(1:shapeArray(1)-1,1:3), dm2(1:shapeArray(1)-1,1:3))
      dp = 0
      dm = 0
      dp2 = 0
      dm2 = 0

      do i=3, shapeArray(1)-2
        ! D'après la relation (23) de Pares et (8) de Abgrall
        d(1) = u(i,2)-u(i+1,2)
        d(2) = 2.5*(u(i+1,3)/u(i+1,1)-u(i,3)/u(i,1))
        d(3) = 1.25*(u(i,3)/u(i,1)+u(i+1,3)/u(i+1,1))*(u(i+1,2)-u(i,2))

        ! ordre 2
        d2(1) = 0.5*(u(i,2)-u(i+2,2))
        d2(2) = 1.25*(u(i+2,3)/u(i+2,1)-u(i,3)/u(i,1))
        d2(3) = 0.5*(u(i+2,2)-u(i,2))*3.5*0.25*(u(i,3)/u(i,1)+2*u(i+1,3)/u(i+1,1)+u(i+2,3)/u(i+2,1))

        do j=1,3
          ! Lax Friedrichs
          ! dp(i,j) = 0.5*d(j)+0.5/sigma*(u(i+1,j)-u(i,j))
          ! dm(i,j) = 0.5*d(j)-0.5/sigma*(u(i+1,j)-u(i,j))

          ! Rusanov
          p = u(i,3)/(0.4*u(i,1))
          pp = u(i+1,3)/(0.4*u(i+1,1))

          bi = sqrt(3.5*p/u(i,1))
          bip = sqrt(3.5*pp/u(i+1,1))

          bip = max(bi, bip)
          ! write(6,*) bip

          ! limiteur de pente
          thetam = u(i,j)-u(i-1,j)
          if (thetam/=0) thetam = (u(i-1,j)-u(i-2,j))/thetam
          thetap = u(i+1,j)-u(i,j)
          if (thetap/=0) thetap = (u(i,j)-u(i-1,j))/thetap

          ! limm = 0 ! 0 = ordre 1, 1 = ordre 2
          ! limp = 0 ! 0 = ordre 1, 1 = ordre 2
          ! limm = minmod(1.d0,thetam)
          ! limp = minmod(1.d0,thetap)
          ! limm = superbee(thetam)
          ! limp = superbee(thetap)
          limm = vonLeer(thetam)
          limp = vonLeer(thetap)

          dp(i,j) = 0.5*d(j)+0.5*bip*(u(i+1,j)-u(i,j))
          dm(i,j) = 0.5*d(j)-0.5*bip*(u(i+1,j)-u(i,j))

          dp2(i,j) = 0.5*d2(j)+0.25*bip*(u(i+2,j)-u(i,j))
          dm2(i,j) = 0.5*d2(j)-0.25*bip*(u(i+2,j)-u(i,j))

          dp(i,j) = dp(i,j) + limp*(dp2(i,j) - dp(i,j))
          dm(i,j) = dm(i,j) + limm*(dm2(i,j) - dm(i,j))
        enddo
      enddo

      do i=2, shapeArray(1)-1
        u(i,:) = ubis(i,:) - sigma*(dp(i-1,:) + dm(i,:))
      enddo
      deallocate(dp, dm, dp2, dm2)

    case default
      write(6,*) "Attention au flux choisi"

    end select
    deallocate(ubis, Fu, fluxm, fluxp, fluxm2, fluxp2)
  end subroutine Iteration



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

    open(unit=15, file="Output/Sol/PCR2_VonLeer_it=" // trim(adjustl(temps)) // ".txt", status="unknown")

    do i=1, shapeArray(1)
      write(x,'(F10.6)') i*dx
      write(hi, '(F10.6)') 1./u(i,1) ! Densité
      write(qi, '(F10.6)') u(i,2) ! Vitesse
      ! write(ei, '(F10.6)') (u(i,3)-0.5*u(i,2)**2)/(0.4*u(i,1)) ! conservative
      write(ei, '(F10.6)') u(i,3)/(0.4*u(i,1)) ! non-conservative
      write(15, *) x // "  " // hi // "  " // qi // "  " // ei

      ! Plan u-p
      ! write(hi, '(F10.6)') u(i,2) ! Vitesse
      ! write(ei, '(F10.6)') (u(i,3)-0.5*u(i,2)**2)/(0.4*u(i,1)) ! Pression
      ! write(15, *) hi // "  " // ei
    enddo

    close(15)
  end subroutine SaveSol








  ! subroutine racines(u,v1,v2)
  !   implicit none
  !   real(kind=8), intent(in), dimension(1:3) :: u
  !   real(kind=8), intent(inout) :: v1, v2
  !   real(kind=8) :: a,b,c,d,g
  !
  !   g = 1.4 - 1
  !   ! if (isNaN(u(1))) stop "Erreur : Valeur propre  =  NaN1"
  !   ! if (isNaN(u(2))) stop "Erreur : Valeur propre  =  NaN2"
  !   ! if (isNaN(u(3))) stop "Erreur : Valeur propre  =  NaN3"
  !   a = g**2*u(1)**2
  !   b = g*(u(3)-0.5*u(2)**2-2*u(2)*u(1))
  !   c = 2.5*u(2)**2 - u(3)
  !
  !   d = b**2-4*a*c
  !   if (d<=0) write(6,*) "Attention au déterminant ! d = ", d
  !   v1 = (-b + sqrt(d))/(2*a)
  !   v2 = (-b - sqrt(d))/(2*a)
  ! end subroutine racines





  ! subroutine SolExacte(xmax, dx, t, u, a)
  !   implicit none
  !   real(kind=8), intent(inout), dimension(:) :: u
  !   real(kind=8), intent(in) :: dx, xmax, t, a
  !   real(kind=8) :: x, sig, mu
  !   integer :: i
  !
  !   ! sig = 0.08
  !   ! mu = 0.5
  !
  !   ! do i = 1, size(u)
  !   !   x = i*dx - a*t
  !   !   if ((x<0.4) .or. x>(0.6)) then
  !   !     u(i) = 0
  !   !   else
  !   !     u(i) = sin((x-0.4)*3.1415*5)
  !   !   endif
  !   !   ! u(i) = 1/(sig*sqrt(2*3.1415))*exp(-((x-mu)/sig)**2/2)
  !   ! enddo
  !
  !   do i = 1, size(u)
  !     x = i*dx - a*t
  !     if (x<0.5) then
  !       u(i) = 1
  !     else
  !       u(i) = 0
  !     endif
  !   enddo
  !
  !
  ! end subroutine SolExacte




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



  real(kind=8) function minmod(a,b) result(c)
    implicit none
    real(kind=8), intent(in) :: a, b

    if (a*b <= 0) then
      c = 0
    else
      if (abs(a)<abs(b)) then
        c = abs(a)
      else
        c = abs(b)
      endif
    endif

  end function minmod



  real(kind=8) function superbee(theta) result(c)
    implicit none
    real(kind=8), intent(in) :: theta

    c = max(0., min(1.,2*theta), min(2.,theta))

  end function superbee


  real(kind=8) function vonLeer(theta) result(c)
    implicit none
    real(kind=8), intent(in) :: theta

    c = (theta+abs(theta))/(1+abs(theta))

  end function vonLeer



end module Solveur
