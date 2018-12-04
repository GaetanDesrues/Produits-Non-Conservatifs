program Main
  use Solveur
  implicit none


  real(kind=8) :: a, x_min, x_max, nbMailles, dt, tF, CFL, phi, dx, sigma, t
  real(kind=8), dimension(:), allocatable :: u, ue

  real(kind=8) :: gg

  ! phi égal à : 0=minmod, 1=minmssod
  a = 1
  x_min = 0
  x_max = 1
  nbMailles = 150
  dt = 0.01
  tF = 1
  CFL = 1
  phi = 1

  dx = (x_max - x_min)/nbMailles

  call CFLSub(CFL, dt, dx, a)
  sigma = a*dt/dx
  ! print*, minmod(a,b)

  allocate(u(1:nbMailles)) ! Initialisation
  call ConditionInitiale(u, dx)

  t=0
  call SaveSol(u, t, dx) ! Save la sol initiale


  do while (t<tF)
    t = t+dt
    call Iteration(u, sigma, phi) ! Calcul de la sol à chaque pas de temps

    call SaveSol(u, t, dx)
  enddo

  ! gg = 0.10
  allocate(ue(1:nbMailles))
  ! ! Solution exacte
  ! call SolExacte(x_max, dx, gg, ue, a)
  ! call SaveSol(ue, tF+gg, dx)
  !
  ! gg = 0.20
  ! call SolExacte(x_max, dx, gg, ue, a)
  ! call SaveSol(ue, tF+gg, dx)
  !
  ! gg = 0.50
  ! call SolExacte(x_max, dx, gg, ue, a)
  ! call SaveSol(ue, tF+gg, dx)


  deallocate(u, ue)

end program Main
