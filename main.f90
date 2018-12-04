program Main
  use Solveur
  implicit none


  real(kind=8) :: x_min, x_max, nbMailles, dt, tF, CFL, phi, dx, sigma, t
  real(kind=8), dimension(:,:), allocatable :: u, ue
  integer :: nbEq, it
  integer, dimension(2) :: shapeArray ! Taille de la solution

  integer, parameter :: LAX_FRIEDRICHS=0, RUSANOV=1 ! Flux
  integer, parameter :: MARCHE=0, GAUSSIENNE=1, SINUS=2 ! Condition Initiale

  real(kind=8) :: gg

  ! phi égal à : 0=minmod, 1=minmssod

  x_min = 0
  x_max = 1
  nbMailles = 100
  dx = (x_max - x_min)/nbMailles
  t = 0
  tF = 0.5
  dt = 0.001
  nbEq = 2
  sigma = dt/dx

  allocate(u(1:nbMailles, 1:nbEq))
  shapeArray = shape(u)

  call ConditionInitiale(u, dx, GAUSSIENNE) ! Initialisation
  call SaveSol(u, 0, dx) ! Save la sol initiale

  it = 0
  do while (t<tF)
    t = t+dt
    it = it + 1
    call Iteration(u, sigma, LAX_FRIEDRICHS) ! Calcul de la sol à chaque pas de temps
    call SaveSol(u, it, dx)
  enddo
  !
  ! ! gg = 0.10
  ! allocate(ue(1:nbMailles))
  ! ! ! Solution exacte
  ! ! call SolExacte(x_max, dx, gg, ue, a)
  ! ! call SaveSol(ue, tF+gg, dx)

  deallocate(u)

end program Main















!
