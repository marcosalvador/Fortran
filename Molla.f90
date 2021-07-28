program harm
implicit none
real    :: massa, dt,ekin,epot
real, dimension (:), allocatable :: pos, vel, vel_parziale, f, kappa
integer :: nstep,it, i, N
write(unit=*,fmt="(a)",advance="no")"N : " 
read*,N
allocate (pos(N), vel(N), vel_parziale(N), f(N), kappa(N))
write(unit=*,fmt="(a)",advance="no")"delta t : "
write(unit=*,fmt="(a)",advance="no")"n.step: "
read*,nstep
write(unit=*,fmt="(a)",advance="no")"massa: "
read*,massa
write(unit=*,fmt="(a)",advance="no")"kappa: "
read*,kappa
write(unit=*,fmt="(a)",advance="no")"pos(0): "
read*,pos
write(unit=*,fmt="(a)",advance="no")"vel(0): "
read*,vel

it=0
write(unit=1,fmt=*)it,it*dt,pos,vel
epot=0
ekin=0
do i=1,N
epot =  0.5 * kappa(i) * pos(i)**2
ekin =  0.5 * massa * vel(i)**2
f(i) = - kappa(i) * pos(i)
end do
write(unit=2,fmt=*)it,dt*it,ekin,epot,ekin+epot

do it = 1,nstep										                                         !algoritmo di Verlet
	epot=0
	do i=1,N
   pos(i) = pos(i) + vel(i) * dt + 0.5* f(i)/massa * dt**2
   vel_parziale(i) = vel(i) + 0.5 * dt * f(i)/massa
   f(i) = - kappa(i) * pos(i)
   epot =  0.5 * kappa(i) * pos(i)**2
   vel(i) = vel_parziale(i) + 0.5 * dt * f(i)/massa
   ekin = 0.5 * massa * vel(i)**2
   end do
   write(unit=1,fmt=*)it,it*dt,pos,vel
   write(unit=2,fmt=*)it,it*dt,ekin,epot,ekin+epot
end do
deallocate(pos, vel, vel_parziale, f, kappa)
end program harm
