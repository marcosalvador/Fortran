!Marco Salvador

module prec
   integer, parameter :: rk=selected_real_kind(8)
end module prec



module funzioni
use prec
implicit none

contains
	
	function f(x,y) result(z)						!funzione da integrare
	real(kind=rk), intent(in) :: x,y
	real(kind=rk) :: z
	z = sqrt(x)/(exp(x-y)-1.0_rk)					!x>=0 ; y<0
	end function



	function trap(fun,j) result (res)				!integrale con il metodo dei trapezi
		integer, parameter :: n = 1000
		real(kind=rk), parameter :: a = 0.0_rk , b = 10.0_rk
		real(kind=rk), intent(in) :: j
		integer :: i
		real(kind=rk) :: res,h
			interface
				function fun(x,y) result(z)
				use prec
				real(kind=rk), intent(in) :: x,y
				real(kind=rk) :: z
				end function fun
			end interface
			
    res = 0
	h = ((b-a)/n)
	do i=1,n-1
	res = res + fun(a+i*h,j)
	end do
	res = h * (fun(a,j)/2 + res + fun(b,j)/2)	
    end function trap
    
    

    function cavsim(fun,j) result(res)				!integrale con il metodo di Cavalieri-Simpson
		integer, parameter :: n = 1000				!n deve essere pari
		real(kind=rk), parameter :: a = 0.0_rk , b = 10.0_rk
		real(kind=rk), intent(in) :: j
		integer :: i
		real(kind=rk) :: res,h
			interface
				function fun(x,y) result(z)
				use prec
				real(kind=rk), intent(in) :: x,y
				real(kind=rk) :: z
				end function fun
			end interface		
	
	res = 0
	h = ((b-a)/n)
	do i=1,n-1
	if (mod(i,2)==0) then
		res = res + 2*fun(a+i*h,j)
	else
		res = res + 4*fun(a+i*h,j)
	end if
	end do
	res = h/3 * (fun(a,j) + res + fun(b,j))
	end function cavsim
    
end module funzioni



module bisezione									!metodo di bisezione per il calcolo degli zeri
use prec											!dell'equazione integrale
use funzioni
implicit none

contains

  	subroutine bis(area,rho,res)					!"area" è un argomento simbolico al posto 
  	integer, parameter :: nmax = 20					!del quale si inserirà "trap" o "cavsim"
  	real(kind=rk), parameter :: tol = 0.01_rk
  	real(kind=rk) :: a,b							![a,b] è l'intervallo dove si cercano le
  	real(kind=rk), intent(in) :: rho				!soluzioni y dell'equazione integrale
  	real(kind=rk), intent(out) :: res 
  	real(kind=rk) :: c
  	integer :: i
  		interface
  			function area(fun,j) result (integ)
  			use prec
  			real(kind=rk), intent(in) :: j
			real(kind=rk) :: integ
				interface
					function fun(x,y) result(z)
					use prec
					real(kind=rk), intent(in) :: x,y
					real(kind=rk) :: z
					end function fun
				end interface
			end function area
  		end interface
  	
  		a = -50.0_rk								!estremi di [a,b]
		b = -0.01_rk
		if ((area(f,a)-rho)*(area(f,b)-rho)>0) then 
			print*, "Non ci sono zeri nell'intervallo [a,b]"
			res = 1/(rho-rho)
		elseif ((area(f,a)-rho) == 0.0_rk) then 
			print*, "y=a è soluzione"
			res = a
		elseif ((area(f,b)-rho) == 0.0_rk) then 
			print*, "y=b è soluzione"
			res = b
		else
			do i=1,nmax
			c=(a+b)/2
			if (abs(b-a)<tol) exit
			if ((area(f,c)-rho) == 0.0_rk) exit
				if ((area(f,a)-rho)*(area(f,c)-rho)<0) then
				b=c
				else
				a=c
				end if
			end do
			res=c
			end if
	
	end subroutine bis

end module bisezione





program integrale
use prec
use funzioni
use bisezione
implicit none
integer :: s
real(kind=rk) :: y,w,k,rho,res

print*, "Inserire il valore di y"					!calcolo dell'integrale generalizzato
read*, y

w = trap(f,y)
print*, "Integrazione trapezoidale", w

k = cavsim(f,y)
print*, "Integrazione Cavalieri-Simpson", k


print*, "Inserire il valore di rho"				!inversione numerica
read*, rho
call bis(trap,rho,res)
print*, "y1 =", res
call bis(cavsim,rho,res)
print*, "y2 =", res


!print*, "Dati bisezione con metodo trapezi"		!dati per graficare l'inversione numerica
!rho = 0.2_rk
!do s = 1, 20
!	rho = rho / s
!	call bis(trap,rho,res)
!	print*, "y =", res
!	write(unit=1,fmt=*) rho,res	
!end do


!print*, "Dati bisezione con metodo Cavalieri-Simpson"
!rho = 0.2_rk
!do s = 1, 20
!	rho = rho / s
!	call bis(cavsim,rho,res)
!	print*, "y =", res
!	write(unit=2,fmt=*) rho,res
!end do



end program integrale