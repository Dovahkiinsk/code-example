module head
parameter(im=461,jm=421,kk=8,ncy=144,nfl=40)
double precision rho0,U,Re,msp,dt,pi,tol,omg
double precision xc,yc,rp,xycy(0:ncy,2),cd_cy,cl_cy
double precision ms_rto,U_red_cy,vecy1,vecy2,dyc,voc,acc_y,vor_up,vor_dn,pow(2),spow
double precision Lf,rho_f,ks,kb,x_ld(2),y_ld(2),xyfl(0:nfl,2,2),vefl(0:nfl,2,2),fbfl(0:nfl,2,2),fbcfl(0:nfl,2,2),cd_f(2),cl_f(2),idfl(0:nfl,2),jdfl(0:nfl,2)
double precision x(im),y(jm),xl,yl
integer iac_cy,jac_cy,nij_cy,nij_f,iac_f(2),jac_f(2),idcy(0:ncy),jdcy(0:ncy)
double precision ux(im,jm),uy(im,jm),rho(im,jm),vor(im,jm)
double precision f(im,jm,0:kk),feq(im,jm,0:kk),ex(0:kk),ey(0:kk)
double precision AA_cy(ncy,ncy),AA_fl(nfl+1,nfl+1,nfl)
double precision ffx(im,jm),ffy(im,jm),xyofl(0:nfl,2,2)

double precision AA_cy1(ncy,25),AA_cy2(ncy,ncy,25),AAx(ncy,25),AAy(ncy,25)  !!  LXJ
integer nn                    !!LXJ
!stream!!!!!!!!!!!!
double precision fm
dimension fm(im,jm,0:kk)
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!velocity-modified
dimension bu_cy(ncy),bv_cy(ncy),fbcy(0:ncy,2)
double precision bu_cy,bv_cy,fbcy
end module
    
program main
use head
use openacc

character (len=20) name2
character*(5) name1
nf=10001
call cpu_time(start_time)
call parainitial
call initialrouandu
call coefficient_cy
call invam(AA_cy,ncy)
call cpu_time(end_time)
print *, "running time is:", end_time-start_time

end
!
subroutine parainitial
use head

double precision eta,disl,dism,disr,dis
double precision St,wr,wn,kn

pi=4.0*atan(1.0)

xl=40.0; yl=32.0
rp=0.5;
xc=0.0; yc=0.0

x_ld(1)=xc+rp*cos(pi*45./180.); y_ld(1)=yc+rp*sin(pi*45./180.)
x_ld(2)=xc+rp*cos(pi*45./180.); y_ld(2)=yc-rp*sin(pi*45./180.)

U=0.1
rho0=1.0
Re=200.0

Lf=1.0
rho_f=0.1*rho0*Lf
ks=2000.*rho0*U**2*Lf

St=0.193		!viv St number
wr=0.5
wn=2.*pi*St*U/wr
kn=1.8751
kb=rho_f*wn**2*(Lf/kn)**4

ms_rto=2.0
U_red_cy=5.0

dyc=0.0
voc=0.0
acc_y=0.0

msp=1./100.
nij_cy=50
nij_f=2*nij_cy

ni1_cy=60; ni2_cy=60+nij_cy+50
nj1_cy=160

iac_cy=140+1!!!!!!!!!!!!圆心坐标
jac_cy=(jm-1)/2+1!!!!!!!圆心坐标
iac_f(1)=iac_cy+35
jac_f(1)=jac_cy+35
iac_f(2)=iac_cy+35
jac_f(2)=jac_cy-35

dt=msp*0.8
tol=3.0*U/(Re*dt)+0.5
omg=1.0/tol
print *, 'single relaxation time:',tol
print *, 'Kb is:',kb
!------------- generate computational mesh -------------!
do i=iac_cy-ni1_cy,iac_cy+ni2_cy
x(i)=2.2*real(i-iac_cy+ni1_cy)/real(ni1_cy+ni2_cy)+xc-0.6
end do
do j=jac_cy-nj1_cy,jac_cy+nj1_cy
y(j)=3.2*real(j-jac_cy+nj1_cy)/real(2*nj1_cy)+yc-1.6
end do

eta=0.85
disl=xl/4.-0.6
disr=3.*xl/4.-1.5
dis=yl/2.-1.6

do i=1,iac_cy-ni1_cy
x(i)=disl*(eta*atan((real(i-1)/real(iac_cy-ni1_cy-1))*tan(1./eta)))-xl/4.
end do
do i=iac_cy+ni2_cy,im
x(i)=x(iac_cy+ni2_cy)+disr*(1.-eta*atan((1.-real(i-iac_cy-ni2_cy)&
    /real(im-iac_cy-ni2_cy))*tan(1./eta)))
end do

do j=1,jac_cy-nj1_cy
y(j)=dis*(eta*atan((real(j-1)/real(jac_cy-nj1_cy-1))*tan(1./eta)))-yl/2.
end do
do j=jac_cy+nj1_cy,jm
y(j)=y(jac_cy+nj1_cy)+dis*(1.-eta*atan((1.-real(j-jac_cy-nj1_cy)&
    /real(jm-jac_cy-nj1_cy))*tan(1./eta)))
end do
!------------- generate cylinder points -------------!
do ip=0,ncy
xycy(ip,1)=xc+rp*cos(2.*pi*real(ip)/real(ncy))
xycy(ip,2)=yc+rp*sin(2.*pi*real(ip)/real(ncy))

idcy(ip)=int((xycy(ip,1)-xc)/msp)+iac_cy
jdcy(ip)=int((xycy(ip,2)-yc)/msp)+jac_cy
end do
!------------- generate filament points -------------!
!do np=1,2
!do ip=0,nfl
!xyfl(ip,1,np)=x_ld(np)+Lf*real(ip)/real(nfl)*cos(0.*pi/180.)
!xyfl(ip,2,np)=y_ld(np)+Lf*real(ip)/real(nfl)*sin(0.*pi/180.)

!idfl(ip,np)=int((xyfl(ip,1,np)-x_ld(np))/msp)+iac_f(np)
!jdfl(ip,np)=int((xyfl(ip,2,np)-y_ld(np))/msp)+jac_f(np)
!end do
!end do

!nij_f=4*nij_cy
!ni1_cy=60; nj1_cy=180
!ni1_f=10; nj1_f=180
!
!iac_cy=140+1
!jac_cy=(jm-1)/2+1
!iac_f=iac_cy+nij_cy
!jac_f=(jm-1)/2+1

!dt=msp
!tol=3.0*U/(Re*dt)+0.5
!omg=1.0/tol
!print *, 'single relaxation time:',tol
!print *, 'Kb is:',kb
!!------------- generate computational mesh -------------!
!do i=iac_cy-ni1_cy,iac_f+ni1_f+nij_f
!x(i)=(0.2+2.*rp+Lf)*real(i-iac_cy+ni1_cy)/&
!    real(iac_f+ni1_f+nij_f-iac_cy+ni1_cy)+xc-0.6
!end do
!do j=jac_cy-nj1_cy,jac_cy+nj1_cy
!y(j)=3.6*real(j-jac_cy+nj1_cy)/real(2*nj1_cy)+yc-1.8
!end do
!
!eta=0.85
!disl=xl/4.-rp-0.6
!disr=3.*xl/4.-2.1
!dis=yl/2.-1.8
!
!do i=1,iac_cy-ni1_cy
!x(i)=disl*(eta*atan((real(i-1)/real(iac_cy-ni1_cy-1))*tan(1./eta)))-xl/4.
!end do
!do i=iac_f+ni1_f+nij_f,im
!x(i)=x(iac_f+ni1_f+nij_f)+disr*(1.-eta*atan((1.-real(i-iac_f-ni1_f-nij_f)&
!    /real(im-iac_f-ni1_f-nij_f))*tan(1./eta)))
!end do
!
!do j=1,jac_f-nj1_f
!y(j)=dis*(eta*atan((real(j-1)/real(jac_f-nj1_f-1))*tan(1./eta)))-yl/2.
!end do
!do j=jac_f+nj1_f,jm
!y(j)=y(jac_f+nj1_f)+dis*(1.-eta*atan((1.-real(j-jac_f-nj1_f)&
!    /real(jm-jac_f-nj1_f))*tan(1./eta)))
!end do
!!------------- generate cylinder points -------------!
!do ip=0,ncy
!xycy(ip,1)=xc+rp*cos(2.*pi*real(ip)/real(ncy))
!xycy(ip,2)=yc+rp*sin(2.*pi*real(ip)/real(ncy))
!
!idcy(ip)=int((xycy(ip,1)-xc)/msp)+iac_cy
!jdcy(ip)=int((xycy(ip,2)-yc)/msp)+jac_cy
!end do
!!------------- generate filament points -------------!
!do ip=0,nfl
!xyfl(ip,1)=x_ld+Lf*real(ip)/real(nfl)*cos(0.*pi/180.)
!xyfl(ip,2)=y_ld+Lf*real(ip)/real(nfl)*sin(0.*pi/180.)
!
!idfl(ip)=int((xyfl(ip,1)-x_ld)/msp)+iac_f
!jdfl(ip)=int((xyfl(ip,2)-y_ld)/msp)+jac_f
!end do

ex(0)= 0.0; ey(0)= 0.0
ex(1)= 1.0; ey(1)= 0.0
ex(2)= 1.0; ey(2)= 1.0
ex(3)= 0.0; ey(3)= 1.0
ex(4)=-1.0; ey(4)= 1.0
ex(5)=-1.0; ey(5)= 0.0
ex(6)=-1.0; ey(6)=-1.0
ex(7)= 0.0; ey(7)=-1.0
ex(8)= 1.0; ey(8)=-1.0
end
!
subroutine initialrouandu
use head

open(2,file='wholefield.dat')
do i=1,im
do j=1,jm
read(2,*) ux(i,j),uy(i,j),rho(i,j)
end do
end do
close(2)

!open(2,file='info_filament.dat')
!do np=1,2
!do ip=0,nfl
!read(2,*) xyfl(ip,1,np),xyfl(ip,2,np),vefl(ip,1,np),vefl(ip,2,np),fbfl(ip,1,np),fbfl(ip,2,np)
           !xyfl(ip,1),xyfl(ip,2),vefl(ip,1),vefl(ip,2),fbfl(ip,1),fbfl(ip,2)

!xyofl(ip,1,np)=xyfl(ip,1,np); xyofl(ip,2,np)=xyfl(ip,2,np)
!end do
!end do
!close(2)
vecy1=0.0; vecy2=0.0
end


subroutine coefficient_cy
use head

double precision temc,dh
double precision drx1,dry1,detax1,detay1,drx2,dry2,detax2,detay2

n=ncy
dh=2.5
 AA_cy1=0.0
 AA_cy2=0.0
!$acc kernels
do ip=1,ncy
   nn=0
!$acc loop seq
do i=idcy(ip)-5,idcy(ip)+5
!$acc loop seq
do j=jdcy(ip)-5,jdcy(ip)+5
    
drx1=(x(i)-xycy(ip,1))/msp
dry1=(y(j)-xycy(ip,2))/msp

if(abs(drx1).le.dh.and.abs(dry1).le.dh) then
    nn=nn+1
    AAx(ip,nn)=x(i)
    AAy(ip,nn)=y(j)
	if(abs(drx1).le.0.5) then
		detax1=3./8.+pi/32.-drx1**2/4.
	elseif(abs(drx1).le.1.5) then
		detax1=1./4.+(1.-abs(drx1))*sqrt(-2.+8.*abs(drx1)&
     		  -4.*drx1**2)/8.-asin(sqrt(2.)*(abs(drx1)-1.))/8.
	else
		detax1=17./16.-pi/64.-3.*abs(drx1)/4.+drx1**2/8.&
     		  +(abs(drx1)-2.)*sqrt(-14.+16*abs(drx1)-4.*drx1**2)/16.&
     		  +asin(sqrt(2.)*(abs(drx1)-2.))/16.
	endif
	if(abs(dry1).le.0.5) then
		detay1=3./8.+pi/32.-dry1**2/4.
	elseif(abs(dry1).le.1.5) then
		detay1=1./4.+(1.-abs(dry1))*sqrt(-2.+8.*abs(dry1)&
     		  -4.*dry1**2)/8.-asin(sqrt(2.)*(abs(dry1)-1.))/8.
	else
		detay1=17./16.-pi/64.-3.*abs(dry1)/4.+dry1**2/8.&
     		  +(abs(dry1)-2.)*sqrt(-14.+16*abs(dry1)-4.*dry1**2)/16.&
     		  +asin(sqrt(2.)*(abs(dry1)-2.))/16.
    endif
    AA_cy1(ip,nn)=detax1*detay1
end if
end do
end do
end do

do ip=1,ncy
    do i=1,25    
	do jp=1,ncy
        
	drx2=(AAx(ip,i)-xycy(jp,1))/msp
	dry2=(AAy(ip,i)-xycy(jp,2))/msp
    
	if(abs(drx2).le.dh.and.abs(dry2).le.dh) then
		if(abs(drx2).le.0.5) then
			detax2=3./8.+pi/32.-drx2**2/4.
		elseif(abs(drx2).le.1.5) then
			detax2=1./4.+(1.-abs(drx2))*sqrt(-2.+8.*abs(drx2)&
     			  -4.*drx2**2)/8.-asin(sqrt(2.)*(abs(drx2)-1.))/8.
		else
			detax2=17./16.-pi/64.-3.*abs(drx2)/4.+drx2**2/8.&
     			  +(abs(drx2)-2.)*sqrt(-14.+16*abs(drx2)-4.*drx2**2)/16.&
     			  +asin(sqrt(2.)*(abs(drx2)-2.))/16.
		endif
		if(abs(dry2).le.0.5) then
			detay2=3./8.+pi/32.-dry2**2/4.
		elseif(abs(dry2).le.1.5) then
			detay2=1./4.+(1.-abs(dry2))*sqrt(-2.+8.*abs(dry2)&
     			  -4.*dry2**2)/8.-asin(sqrt(2.)*(abs(dry2)-1.))/8.
		else
			detay2=17./16.-pi/64.-3.*abs(dry2)/4.+dry2**2/8.&
     			  +(abs(dry2)-2.)*sqrt(-14.+16*abs(dry2)-4.*dry2**2)/16.&
     			  +asin(sqrt(2.)*(abs(dry2)-2.))/16.
		endif
		AA_cy2(ip,jp,i)=0.5*detax2*detay2*dt/msp/msp
    end if
    end do 
    end do 
end do

do i=1,ncy
do j=1,ncy
AA_cy(i,j)=0.0
end do
end do

do ip=1,ncy
   do jp=1,ncy
    !$acc loop seq
     do i=1,25
        AA_cy(ip,jp)= AA_cy(ip,jp)+AA_cy1(ip,i)*AA_cy2(ip,jp,i)
end do
end do
end do
!$acc end kernels

end subroutine coefficient_cy

subroutine invam(a,n)
use head

dimension a(n,n),is(n),js(n)
double precision a,dd,t

l=1

do 100 k=1,n
dd=0.
do 10 i=k,n
do 10 j=k,n
if(abs(a(i,j)).gt.dd) then
	dd=abs(a(i,j))
	is(k)=i
	js(k)=j
endif
10 continue
if(dd+1.0.eq.1.0) then
	l=0
	print *, 'err, not inv'
!	return
endif
do 30 j=1,n
t=a(k,j)
a(k,j)=a(is(k),j)
a(is(k),j)=t
30 continue
do 40 i=1,n
t=a(i,k)
a(i,k)=a(i,js(k))
a(i,js(k))=t
40 continue
a(k,k)=1./a(k,k)
do 50 j=1,n
if(j.ne.k) then
	a(k,j)=a(k,j)*a(k,k)
endif
50 continue
do 70 i=1,n
if(i.ne.k) then
	do 60 j=1,n
	if(j.ne.k) then
		a(i,j)=a(i,j)-a(i,k)*a(k,j)
	endif
60	continue
endif
70 continue
do 80 i=1,n
if(i.ne.k) then
	a(i,k)=-a(i,k)*a(k,k)
endif
80 continue
100 continue
do 130 k=n,1,-1
do 110 j=1,n
t=a(k,j)
a(k,j)=a(js(k),j)
a(js(k),j)=t
110 continue
do 120 i=1,n
t=a(i,k)
a(i,k)=a(i,is(k))
a(i,is(k))=t
120 continue
130 continue

return
end subroutine invam
!
!

