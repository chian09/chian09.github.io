program cylinder
real, dimension (:,:), allocatable :: A,B,C,T_tk,temp,ti,temp_sol, T_tk1, temp1year, temp10year, temp50year, temp100year, r, D
real, dimension (:,:), allocatable ::source_vector,source
double precision  dr, dt, kappa, s, isol,r1,r2,rc,a_source, Trod, Tau_0,t1
integer T, n, m,i, j,k,p,x,y, IER, PGBEG
r1=0
rc=100
r2=rc
T=100
n=9
m=100
dr= real(rc)/real(n+1)
dt= real(T)/real(m)
t1=0.0
kappa= 2*10**(7)
s= real(kappa*dt)/real(dr**2)
print *, kappa,s, dt

allocate(A(n,n))
A=0.0
do i=1,n
do j=1,n
if (i .eq. j) then
A(i,j) = 1 + 2*s
else if (i .eq. j+1) then
A(i,j) = -s +real(s)/real(2*i)
else if (i .eq. j-1) then
A(i,j) = -s - real(s)/real(2*i)
end if
end do
end do
allocate(D(n,1))
D=0
allocate(T_tk(n,1))
T_tk = 300.0
!print*, T_tk

allocate(T_tk1(n,1))
T_tk1=0
i=1
D(i,1) = -s+(s/2)*T_tk1(1,1)!neumann T_0=T_1
i=n
D(i,1)= -s-(s/(2*i))*300


allocate(temp1year(n,1))
temp1year=0
allocate(temp10year(n,1))
temp10year=0
allocate(temp50year(n,1))
temp50year=0
allocate(temp100year(n,1))
temp100year=0

!r1= r1+dr
!r2=r2-dr
!x= (r2-r1)/n
!print*, r1
!allocate(r(x+1,1))
!do y=1,x+1
!r(y,1) = r1
!r1=r1+dr
!print*, r1, r(y,1)
!end do
!print *,D
allocate(source_vector(n,1))
source_vector = 0;
allocate(source(n,1))
source= 0;

a_source = 25;   !a = 25cm radius of rod is 25cm
Trod = 1;  !initial temperature change due to nuclear rod is 1K
tau_0 = 100;  !half-life of rod is 100 years

do i = 1,n
    if (r(i,1).lt. a_source) then
       source_vector(i,1) = 1;
    end if
end do

isol = 1
allocate(temp_sol(4,1))
temp_sol = 0 


do k= 1,m
t1= t1 + dt
do i=1,n
source(i,1) = (real(Trod*exp(-real(t1)/real(tau_0)))/real(a_source**2))*source_vector(i,1);
end do
do i=1,n
B(i,1)= D(i,1)+T_tk(i,1)+kappa*dt*source(i,1)
end do
call gauss_1(A,B,T_tk1,n)
!if (k.eq.10)then
!print*,T_tk1
!end if

T_tk= T_tk1
!ti(k,1) = t1

if (k.eq.10) then
       temp_sol(1,1) = t1
       temp1year = T_tk1
elseif (k.eq.40) then
       temp_sol(2,1) = t1
       
       temp10year = T_tk1
elseif (k.eq.80) then
       temp_sol(3,1) = t1
      
       temp50year = T_tk1
elseif (k.eq.100) then
       temp_sol(4,1) = t1
       
       temp100year = T_tk1
end if
end do
do i= 1,n

print*,temp1year(i,1)
end do
!IER = PGBEG(0,'?',1,1)
!IF (IER.NE.1) STOP
!CALL PGENV(0.,100.,300.,301.,0,1) 
!CALL PGLAB('(x)', '(y)', 'A Simple Graph')
!CALL PGLINE(9,r,temp1year)
!CALL PGLINE(9,r,temp10year)
!CALL PGLINE(9,r,temp50year)
!CALL PGLINE(9,r,temp100year)
!CALL PGEND
end



 subroutine gauss_1(a,b,x,n)
!============================================================
! Solutions to a system of linear equations A*x=b
! Method: the basic elimination (simple Gauss elimination)
! Alex G. November 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! b(n)   - vector of the right hand coefficients b
! n      - number of equations
! output ...
! x(n)   - solutions
! comments ...
! the original arrays a(n,n) and b(n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real a(n,n), b(n), x(n)
real c
integer i, j, k

!step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      c=a(i,k)/a(k,k)
      a(i,k) = 0.0
      b(i)=b(i)- c*b(k)
      do j=k+1,n
         a(i,j) = a(i,j)-c*a(k,j)
      end do
   end do
end do

!step 2: back substitution
x(n) = b(n)/a(n,n)
do i=n-1,1,-1
   c=0.0
   do j=i+1,n
     c= c + a(i,j)*x(j)
   end do 
   x(i) = (b(i)- c)/a(i,i)
end do
end subroutine gauss_1
