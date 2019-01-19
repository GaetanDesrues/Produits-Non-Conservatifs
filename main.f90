program Main
  use Solveur
  implicit none


  real(kind=8) :: x_min, x_max, nbMailles, dt, tF, CFL, phi, dx, sigma, t
  real(kind=8), dimension(:,:), allocatable :: u, ue
  integer :: nbEq, it, itSave
  integer, dimension(2) :: shapeArray ! Taille de la solution

  integer, parameter :: LAX_FRIEDRICHS=0, RUSANOV=1, PARES=2 ! Flux
  integer, parameter :: MARCHE=0, GAUSSIENNE=1, SINUS=2 ! Condition Initiale

  real(kind=8) :: gg

  x_min = 0
  x_max = 1
  nbMailles = 1500
  dx = (x_max - x_min)/nbMailles
  t = 0
  tF = 0.2
  dt = 0.0001
  nbEq = 3
  sigma = dt/dx

  allocate(u(1:nbMailles, 1:nbEq))
  shapeArray = shape(u)

  call ConditionInitiale(u, dx, MARCHE) ! Initialisation
  call SaveSol(u, 0, dx) ! Save la sol initiale

  it = 0
  itSave = 0
  do while (t<tF)
    t = t+dt
    it = it + 1
    call Iteration(u, sigma, RUSANOV) ! Calcul de la sol à chaque pas de temps

    if(mod(it,20)==0) then
      itSave = itSave + 1
      call SaveSol(u, itSave, dx)
    endif
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
