module gravitazione
implicit none
integer,parameter :: kr=selected_real_kind(15), nbody=3
real(kind=kr),parameter :: G=6.673E-11_kr
real(kind=kr),dimension(nbody) :: mass

contains
  subroutine interazione(pos,f,epot)
   real(kind=kr), intent(in), dimension(:,:) :: pos
   real(kind=kr), intent(out) :: epot
   real(kind=kr), intent(out), dimension(:,:) :: f
   real(kind=kr), dimension(size(pos,1)) :: posij
   real(kind=kr) :: rij
   integer ::i,j
   epot = 0
   f    = 0
   do i=1,nbody
      do j=1,nbody
         if( i==j ) cycle
         posij(:)  = pos(:,j)-pos(:,i)
         rij=sqrt( dot_product(posij,posij) )
         epot   = epot   -(G*mass(i)*mass(j))/rij
         f(:,i) = f(:,i) +(G*mass(i)*mass(j)) * posij(:)/rij**3
      end do
   end do
   epot = epot/2    
   
  end subroutine interazione
end module gravitazione


program principale
use gravitazione
implicit none

real(kind=kr) :: dt,ekin,epot
real(kind=kr), dimension(3,nbody) :: pos,vel,f
real, dimension(3) :: cdm,vcm
integer :: z,nstep,it,i

print*, "scegli il set di dati"
read*, z

if (z==0) then
	
	write(unit=*,fmt="(a)",advance="no")"delta t : "
	read*,dt
	write(unit=*,fmt="(a)",advance="no")"n.step: "
	read*,nstep
	write(unit=*,fmt="(a)",advance="no")"pos(0): "
	read*,pos
	write(unit=*,fmt="(a)",advance="no")"vel(0): "
	read*,vel
	write(unit=*,fmt="(a)",advance="no")"mass: "
	read*,mass

else
	
	dt=86400_kr								!Sole, Terra, Giove
	nstep=3650_kr
	pos(:,1)=[0.0_kr,0.0_kr,0.0_kr]
	pos(:,2)=[152.1e9_kr,0.0_kr,0.0_kr]
 	pos(:,3)=[816.62e9_kr,0.0_kr,0.0_kr]
 	vel(:,1)=[0.0_kr,0.0_kr,0.0_kr]
	vel(:,2)=[0.0_kr,29.3e3_kr,0.0_kr]
 	vel(:,3)=[0.0_kr,12.4e3_kr,0.0_kr]
	mass(1)=1.989e30_kr
	mass(2)=5.972e24_kr
	mass(3)=1.e30_kr
	
end if
	
	
do i=1,3									!posizione relativa rispetto a quella del centro di massa
cdm(i) = sum(mass*pos(i,:)/sum(mass))
end do
do i=1,nbody
pos(:,i) = pos(:,i) - cdm(:)
end do


do i=1,3									!velocità relativa rispetto a quella del centro di massa
vcm(i) = sum(mass*vel(i,:)/sum(mass))
end do
do i=1,nbody
vel(:,i) = vel(:,i) - vcm(:)
end do
	

it=0
write(unit=1,fmt=*)it,it*dt,pos,vel
ekin = 0.5 * sum(spread(mass,1,3)*vel**2)
call interazione(pos,f,epot)
write(unit=2,fmt=*)it,dt*it,ekin,epot,ekin+epot

do it = 1,nstep
	pos = pos + vel * dt + 0.5 * f/spread(mass,1,3) * dt**2
    vel = vel + 0.5 * dt * f/spread(mass,1,3)
    call interazione(pos,f,epot)
    vel = vel + 0.5 * dt * f/spread(mass,1,3)
    write(unit=1,fmt=*)it,it*dt,pos,vel
    ekin = 0.5 * sum(spread(mass,1,3)*vel**2)
    write(unit=2,fmt=*)it,it*dt,ekin,epot,ekin+epot
end do
end program principale
