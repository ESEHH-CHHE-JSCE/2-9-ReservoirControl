! mod_variables
module mod_variables
implicit none

integer maxhs, dt
real, allocatable:: st(:), ht(:)
real ca, a

end module mod_variables

! runge.f90
program runge
use mod_variables
implicit none

integer maxt, t, i, maxtt, tt
real, allocatable:: s(:)
real, allocatable:: inflow(:), outflow(:)
real k0, k1, k2, k3, in
real f, q, h, spre, spos

! STEP 0: Reading Parameters

open(10, file = "param.txt", status = 'old')
read(10, *) maxt
read(10, *) dt
read(10, *) maxhs
read(10, *) ca
read(10, *) a
close(10)

allocate( s(0:maxt), inflow(0:maxt), outflow(0:maxt) )
allocate( ht(0:maxhs), st(0:maxhs) )

! STEP 1: Reading Data
open(15, file = "inflow.txt", status = 'old')
do t = 1, maxt
 read(15, *) inflow(t) ! Hourly inflow [m3/s]
enddo
close(15)

open(20, file = "HS_table.txt", status = 'old')
do i = 1, maxhs
 read(20, *) ht(i), st(i) ! Water level [m], Storage [m3]
enddo
close(20)
ht(0) = 0.
st(0) = 0.

! STEP 2: Runge-Kutta
s(0) = 0
maxtt = 3600 / dt
do t = 1, maxt
 spre = s(t-1)
 do tt = 1, maxtt
  in = inflow(t)
  k0 = f(in, spre)
  k1 = f(in, spre + k0 / 2.0 * dt)
  k2 = f(in, spre + k1 / 2.0 * dt)
  k3 = f(in, spre + k2 * dt)
  spos = spre + dt * (k0 + 2 * k1 + 2 * k2 + k3) / 6.
  if(spos .lt. 0) spos = 0.
  spre = spos
 enddo
 s(t) = spos
enddo

! STEP 3: Output
open(30, file = "out.txt")
write(*,'("    t         i(t)         s(t)         h(t)         q(t)")')
do t = 1, maxt
 write(30, '(i5, 4f13.3)') t, inflow(t), s(t), h(s(t)), q(s(t))
 write(*, '(i5, 4f13.3)') t, inflow(t), s(t), h(s(t)), q(s(t))
enddo
close(30)

end program runge

! FUNCTION f
function f(in, s)
use mod_variables
implicit none
integer i
real in, s, f, q
f = in - q(s)
return
end

! FUNCTION q
function q(s)
use mod_variables
implicit none
real s, q, h, hh
hh = h(s)
if(hh.lt.0) hh = 0.
q = a * Ca * sqrt(2.0 * 9.81 * hh)
return
end

! FUNCTION h
function h(s)
use mod_variables
implicit none
real s, h
integer i
if(s .gt. st(maxhs))then
 write(*,*) s, st(maxhs)
 stop "error."
endif
do i = 1, maxhs
 if(st(i) .gt. s ) then
  h = ht(i-1) + (ht(i) - ht(i-1)) * (s - st(i-1)) / (st(i) - st(i-1))
  exit
 endif
enddo
return
end
