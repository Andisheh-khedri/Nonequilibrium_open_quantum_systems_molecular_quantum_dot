!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! **********************************************************************
! This program aims to calculate the Observables 
! for SRLM coupled to phononic bath at a given Temperature! 
! using Keldysh FRG and also the comparison to Perturbation Theory!
! Using the ODEPACK library to solve the flow equations.
! Compared to the previous version using dqawc to perform the hilbert transforms!
!##################################################################
program KeldyshFRG_SRLMCP
!------------------------------------------------------------------
implicit none
!------------------------------------------------------------------
common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon&
       ,dnespec,dneSigmaRre,dneSigmaRim,dneSigmaK
!------------------------------------------------------------------
!  Parameter declaration:
!------------------------------------------------------------------
INTEGER,parameter::iw1max=200,NEQ=4*iw1max,NTAB=20000

DOUBLE PRECISION::zero,pi,w0,lambda,gama,e0,ep,betaL,betaR&
                 ,gamaR,gamaL,muL,muR,beta,betaph,betaaux,voltage
INTEGER:: openstatus,Flag, ISTATE, IOPT, LRW, LIW, MF,MU, ML,k,i,j,ibeta,ie0int,ivoltage,NEQ2
DOUBLE PRECISION,DIMENSION(800)::YFRG,YDFRG,YFRGcommon
DOUBLE PRECISION::x, xout, RTOL, ATOL,JAC,xcommon,fcoef,y,yout
DOUBLE PRECISION,allocatable,dimension(:)::RWORK
INTEGER,allocatable,dimension(:)::IWORK
!------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(200)::w1,NdisFRG
DOUBLE PRECISION,DIMENSION(2*iw1max)::Ymat,YDmat
DOUBLE PRECISION,DIMENSION(4)::Ycur,YDcur
!------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(NEQ/4)::rananum,rnum,rFRG
COMPLEX*16,DIMENSION(NEQ/4)::sRana,sKana,ana,sRananon,sKananon
COMPLEX*16::ic,Intana,self,dp,dm,dtp,dtm,IntCaushy
DOUBLE PRECISION::ap,am,SHart,fermi,fermiaux
INTEGER*8::sgnn
!------------------------------------------------------------------
! Unnecessary:
!------------------------------------------------------------------
INTEGER,parameter::ITOL=1,nrowpd=2,NEQ1=1
DOUBLE PRECISION,DIMENSION(nrowpd,NEQ1)::pd1
!------------------------------------------------------------------
REAL*8, external::cpv
!------------------------------------------------------------------
! Parameters for Quadpack:
!------------------------------------------------------------------
DOUBLE PRECISION::wmin,wmax,wpole,abserr,epsabs,epsrel,resultquad
INTEGER,PARAMETER::limit=100,lenw=4*limit
INTEGER,DIMENSION(limit)::Iworkquad
INTEGER::last,neval,ier
DOUBLE PRECISION::workquad(lenw)
DOUBLE PRECISION, EXTERNAL::fKimcaushy,fRArecaushy,fKrecaushy,fRAimcaushy
!------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(0:200)::vol
DOUBLE PRECISION,DIMENSION(1:3)::Jch,JQR,JQL
DOUBLE PRECISION::G0,deltavoltage,deltadeltaTem,deltaTem&
                ,Gdif,Sdif,Kedif,PFdif,ZTdif,Ldif&
                ,dJchdv,dJchddeltaT,dJQdv,dJQddeltaT,dJQLdv,dJQLddeltaT,Jcharge,JheatR,JheatL&
                ,dJchdvdec,dJchddeltaTdec,dJQRdvdec,dJQRddeltaTdec,dJQLdvdec,dJQLddeltaTdec&
                ,Gdifdec,Sdifdec,Kedifdec,PFdifdec,slogGdifdec&
                ,ZTdifdec,Ldifdec,etaZT,etaZTdec,etaCA,eta,powerlinear,power,powerlineardec,fanoana
INTEGER::iparder,ideltaTem,ilambda,ipardermax
!------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:3)::YInL,YDINL
DOUBLE PRECISION,DIMENSION(1:4)::YInR,YDINR
INTEGER::NEQ4
DOUBLE PRECISION,DIMENSION(1:3,1:iw1max)::nespec,neSigmaRre,neSigmaRim,neSigmaK
DOUBLE PRECISION,DIMENSION(1:iw1max)::dnespec
DOUBLE PRECISION,DIMENSION(1:iw1max)::dneSigmaRre,dneSigmaRim,dneSigmaK
DOUBLE PRECISION,DIMENSION(1:3)::YIncor,YDIncor,YIncorvolt,YIncordeltaT,YInnedis,YDInnedis,YInnedisdeltaT,YInnedisvolt
!------------------------------------------------------------------
real :: start_time, stop_time
!------------------------------------------------------------------
call cpu_time(start_time)
!------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
open(unit=2,file="kk17datal0.5w020gTL0.1gdT2.0ge0epp2gVscanrest3Aug.txt",status="new"&
                               ,action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
!------------------------------------------------------------------
!  Parameter of the model:
!------------------------------------------------------------------
zero=0.0d+00
pi=atan2(zero,-gama)
ic=cmplx(0.0d+00,1.0d+00)
JAC=0.0d+00

w0=1.0d+00!gama*20.0d+00!20.0d+00
gama=w0/20.0d+00!2.0d+00*w0*(sqrt(0.5d+00))!w0/20.0d+00!0.0001!w0/1.5d+00
!w0=gama/10.0d+00

gamaR=gama/2.0d+00!10.0d+00*gama/11.0d+00
gamaL=gama/2.0d+00!1.0d+00*gama/11.0d+00

G0=4.0d+00*gamaL*gamaR/(gama**2)

fcoef=1.0d+00
!------------------------------------------------------------------
Do ilambda=1,1

deltavoltage=1.0d+00*gama*(10.0d+00**(-4.0d+00+(4.0d+00*(1-1)/40.0d+00)))!0.0005d+00!
deltadeltaTem=1.0d+00*gama*(10.0d+00**(-4.0d+00+(4.0d+00*(1-1)/40.0d+00)))!0.0005d+00

lambda=0.5d+00*w0
ep=(lambda**2)/w0
!------------------------------------------------------------------
! Gate Voltage:
!------------------------------------------------------------------
Do ie0int=1,1!01
!------------------------------------------------------------------
!e0=ep!+30.0d+00*gama!*(0.0d+00+10.0d+00*((ie0int-1)/80.0d+00))!0.0d+00!ep+1.0d+00*gama!ep!+(4.0d+00*gama*(ie0int-1))!*(0.0d+00+80.0d+00*((ie0int-1)/100.0d+00))!*(10.0d+00**(-2.0d+00+4.0d+00*((ie0int-1)/40.0)))!-gama!+6.0d+00*gama!-2.0d+00*gama!-2.0d+00*gama!!+16.0d+00*gama!(1.0d+00*w0)!+100.0d+00*gama!-w0
e0=ep+2.0d+00*gama
!------------------------------------------------------------------
! Temperature:
!------------------------------------------------------------------
Do ibeta=1,1
!------------------------------------------------------------------
!beta=(1.0d+00/((10.0d+00**(-2.0d+00+4.0d+00*((ibeta-1)/14.0d+00)))*gama))!*(10.0d+00**(-2.0d+00+(4.0d+00*(ibeta-1)/10.0d+00)))!100000.0d+00
beta=1.0d+00/(gama*0.1d+00)
betaaux=beta
betaph=beta

print*,"Temperatures",1.0d+00/(gama*beta),1.0d+00/(w0*beta)
!------------------------------------------------------------------
! Temperature Gradient:
!------------------------------------------------------------------
Do ideltaTem=1,1
!------------------------------------------------------------------
!------------------------------------------------------------------
! Bias Voltage:
!------------------------------------------------------------------
Do ivoltage=1,2
!------------------------------------------------------------------
if(ivoltage==1)then
ipardermax=1!2!2!3
deltaTem=2.0d+00*gama!*(0.001d+00+20.0d+00*((ideltaTem-1)/80.0d+00))!2.0d+00*gama!
else
ipardermax=1!2!2!1
deltaTem=2.0d+00*gama!*(0.001d+00+20.0d+00*((ideltaTem-1)/80.0d+00))!2.0d+00*gama!1.0d+00*gama*(0.001d+00+10.0d+00*((ideltaTem-1)/40.0d+00))!
endif
!------------------------------------------------------------------
if(ivoltage<11)then
!voltage=1.0d+00*gama*(10**(-1.0d+00+((2.0d+00)*(ivoltage-1)/(10.0d+00))))
voltage=+0.0d+00+3.0d+00*w0*(ivoltage-1.0d+00)/(10.0d+00)*exp((ivoltage-10.0d+00)/9.d+00)
!elseif(ivoltage<12)then
! voltage=w0*1.0d+00*((3.0d+00*gama/w0)+(3.0d+00-((3.0d+00*gama/w0)))*((ivoltage-5)/10.0d+00))
else
voltage=3.0d+00*w0+(6.0d+00*W0-3.0d+00*w0)*(1.0d+00*((ivoltage-11)/4.0d+00))
endif

voltage=gama*(2.5d+00+((ivoltage-1)*0.5/1.0d+00))

print*,"note","ie0int=",ie0int,"ivoltage=",ivoltage
!------------------------------------------------------------------
! Infinitesimal changes in Temperature gradient and bias Voltage 
!                to calculate partial derivative:
!------------------------------------------------------------------
Do iparder=1,ipardermax
!------------------------------------------------------------------
if(iparder==2)then
voltage=voltage+deltavoltage
elseif(iparder==3)then
deltaTem=deltaTem+deltadeltaTem
voltage=voltage-deltavoltage
endif
!------------------------------------------------------------------
betaL=1.0d+00/((1.0d+00/beta)) !+(deltaTem/2.0d+00)
betaR=1.0d+00/((1.0d+00/beta)+(deltaTem/1.0d+00)) 

print*,"Temperatures1",beta,betaL,betaR,deltaTem*betaR
!------------------------------------------------------------------
vol(ivoltage)=voltage

muL=voltage*(gamaR/gama)!/2.0d+00!(gamaR/gama)*-0.01*gama!0.0d+00
muR=-voltage*(gamaL/gama)!/2.0d+00!(gamaL/gama)*0.01*gama!0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
! Incoming energy grid:
!------------------------------------------------------------------
Do k=1,iw1max
if(k<(iw1max/2)+1)then

if(k<41)then
w1(k)=(e0-ep)+1.5d+00*w0-28.0d+00*w0*(((k-40.5)*exp((abs(40.5-k)-40.0d+00)/10.0d+00)/40.0d+00))
elseif(k<61)then
w1(k)=(e0-ep)+w0-(0.5d+00)*w0*(((k-50.5)*exp((abs(50.5-k)-10.0d+00)/8.0d+00)/10.0d+00))
else
w1(k)=(e0-ep)+(w0/2.0d+00)*(exp(-(k-61.0d+00)*0.18d+00))
endif

!if(k<((iw1max/6)+1))then
!w1(k)=1.5d+00*w0-28.0d+00*w0*(((k-50.5)*exp((abs(50.5-k)-50.0d+00)/10.0d+00)/50.0d+00))
!elseif(k<((iw1max/10)+(2*iw1max/6)+1))then
!w1(k)=w0-(0.5d+00)*w0*(((k-90.5)*exp((abs(90.5-k)-40.0d+00)/8.0d+00)/40.0d+00))
!else
!w1(k)=(w0/2.0d+00)*(exp(-(k-131.0d+00)*0.3d+00))
!endif
else                                     
w1(k)=(e0-ep)-(w1(iw1max-k+1)-(e0-ep))
endif
!write(1,*)k,w1(k)/w0,(w1(k)-w0)/w0
!------------------------------------------------------------------
end do
!------------------------------------------------------------------
!------------------------------------------------------------------
! FRG: Sloving Flow equations:
!------------------------------------------------------------------
!------------------------------------------------------------------
MF=10
!------------------------------------------------------------------
! Dimension Declaration of WORK:
!------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ+NEQ**2 
LIW=20+NEQ
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ+(2*ML+MU)*NEQ
LIW=20+NEQ
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=10000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-8.0d+00)!(-16.0d+00)   !
ATOL=w0*10.0d+00**(-8.0d+00)!(-16.0d+00)   !
ISTATE=1
x=100.0d+00*w0!w1(1)
xout=0.0d+00!0.01*w1(NEQ/12)
ml=1
mu=1

YFRG(1:NEQ/4)=0.0d+00                !Imaginary Part of Retarded
YFRG((NEQ/4)+1:(2*NEQ/4))=e0-ep      !Real Part of Retarded

YFRG(1+(2*NEQ/4):3*NEQ/4)=0.0d+00    !Imaginary Part of Keldysh
YFRG(1+(3*NEQ/4):(NEQ))=0.0d+00      !Real Part of Keldysh

YFRGcommon=YFRG
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL FFRG (NEQ, x, YFRG, YDFRG)

!CALL JAC (NEQ, x, Yfrg, ml, mu, pd, nrowpd)

CALL DLSODE (FFRG, NEQ, YFRG, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"FRG",ISTATE,xout/w0,RTOL/gama

deallocate(RWORK,IWORK)
!-----------------------------------------------------------------------
Do k=1,iw1max
IntCaushy=0.0d+00
!-----------------------------------------------------------------------
! Checking the Caushy Principal Value:
!-----------------------------------------------------------------------
epsabs=w0*(10.0d+00)**(-5.0d+00)
epsrel=w0*(10.0d+00)**(-5.0d+00)
wpole=w1(k)-w0
wmin=w1(k)-w0-20.0d+00*w0
wmax=w1(k)-w0+20.0d+00*w0

! CALL dqawc(fRArecaushy,wmin,wmax,wpole,epsabs,epsrel,resultquad,abserr,neval,ier,limit,lenw,last,iworkquad,workquad)

! IntCaushy=resultquad



epsabs=w0*(10.0d+00)**(-5.0d+00)
epsrel=w0*(10.0d+00)**(-5.0d+00)
wpole=w1(k)-w0
wmin=w1(k)-w0-20.0d+00*w0
wmax=w1(k)-w0+20.0d+00*w0

! CALL dqawc(fRAimcaushy,wmin,wmax,wpole,epsabs,epsrel,resultquad,abserr,neval,ier,limit,lenw,last,iworkquad,workquad)   

! IntCaushy=IntCaushy+ic*resultquad

! write(1,*)k,w1(k)/w0,real(IntCaushy),aimag(IntCaushy)&
!          ,real(-pi*((1.0d+00/((w1(k)-w0)+ic*gama))**2)),aimag(-pi*((1.0d+00/((w1(k)-w0)+ic*gama))**2))
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Perturbation Theory (zero T):
!-----------------------------------------------------------------------
SHart=0.0d+00!-ep
ap=1.0d+00/(gama**2+(e0+SHart-w0-w1(k))**2)
am=-1.0d+00/(gama**2+(e0+SHart+w0-w1(k))**2)


sRana(k)=ic*((lambda**2)/2.0d+00)*((gama*am*tanh(beta*(w1(k)-w0)/2.0d+00)+gama*ap*tanh(beta*(w1(k)+w0)/2.0d+00)))&
        -((lambda**2)/(2.0d+00*pi))*(-gama*(ap+am)*log((e0+SHart)**2+gama**2)&
                          +2*gama*(ap*log(abs(w1(k)+w0))+am*log(abs(w1(k)-w0)))&
                          +2*(ap*(w1(k)+w0-(e0+SHart))+am*(w1(k)-w0-(e0+SHart)))*atan2((e0+SHart),gama))&
                          +((lambda**2)/2.0d+00)*((1.0d+00/(w1(k)-w0-(e0+SHart)+ic*gama))&
                                             +(1.0d+00/(w1(k)+w0-(e0+SHart)+ic*gama)))&
         +e0+((-ep)*(1.0d+00-(2.0d+00/pi)*atan(e0/gama)))

!-----------------------------------------------------------------------

SRananon(k)=e0+((-ep)*(gamaL/gama)*(1.0d+00-(2.0d+00/pi)*atan((e0-muL)/gama)))&
            +((-ep)*(gamaR/gama)*(1.0d+00-(2.0d+00/pi)*atan((e0-muR)/gama)))&
            +(lambda**2)*(gamaL/(gama**2+(w1(k)+w0-e0)**2))&
                        *(((1.0d+00/pi)*(log((sqrt(gama**2+(e0-muL)**2))/(abs(w1(k)+w0-muL)))))&
                         +(((w1(k)+w0-e0)/(2.0d+00*gama))*(1+(-2.0d+00/pi)*atan((e0-muL)/gama))))&
            +(lambda**2)*(gamaL/(gama**2+(w1(k)-w0-e0)**2))&
                        *(((-1.0d+00/pi)*(log((sqrt(gama**2+(e0-muL)**2))/(abs(w1(k)-w0-muL)))))&
                         +(((w1(k)-w0-e0)/(2.0d+00*gama))*(1+(2.0d+00/pi)*atan((e0-muL)/gama))))&
            +(lambda**2)*(gamaR/(gama**2+(w1(k)+w0-e0)**2))&
                        *(((1.0d+00/pi)*(log((sqrt(gama**2+(e0-muR)**2))/(abs(w1(k)+w0-muR)))))&
                         +(((w1(k)+w0-e0)/(2.0d+00*gama))*(1+(-2.0d+00/pi)*atan((e0-muR)/gama))))&
            +(lambda**2)*(gamaR/(gama**2+(w1(k)-w0-e0)**2))&
                        *(((-1.0d+00/pi)*(log((sqrt(gama**2+(e0-muR)**2))/(abs(w1(k)-w0-muR)))))&
                         +(((w1(k)-w0-e0)/(2.0d+00*gama))*(1+(2.0d+00/pi)*atan((e0-muR)/gama))))&
            -ic*(lambda**2)*(gamaL/(gama**2+(w1(k)+w0-e0)**2))*((1.0d+00-tanh(beta*(w1(k)+w0-muL)/2.0d+00))/2.0d+00)&
            -ic*(lambda**2)*(gamaL/(gama**2+(w1(k)-w0-e0)**2))*((1.0d+00+tanh(beta*(w1(k)-w0-muL)/2.0d+00))/2.0d+00)&
            -ic*(lambda**2)*(gamaR/(gama**2+(w1(k)+w0-e0)**2))*((1.0d+00-tanh(beta*(w1(k)+w0-muR)/2.0d+00))/2.0d+00)&
            -ic*(lambda**2)*(gamaR/(gama**2+(w1(k)-w0-e0)**2))*((1.0d+00+tanh(beta*(w1(k)-w0-muR)/2.0d+00))/2.0d+00)
!-----------------------------------------------------------------------
SKananon(k)=-ic*(lambda**2)*(gamaL/(gama**2+(w1(k)-w0-e0)**2))&
                           *(1.0d+00+((tanh(beta*(w1(k)-w0-muL)/2.0d+00))/(tanh(beta*w0/2.0d+00))))&
            -ic*(lambda**2)*(gamaL/(gama**2+(w1(k)+w0-e0)**2))&
                           *(-1.0d+00+((tanh(beta*(w1(k)+w0-muL)/2.0d+00))/(tanh(beta*w0/2.0d+00))))&
            -ic*(lambda**2)*(gamaR/(gama**2+(w1(k)-w0-e0)**2))&
                           *(1.0d+00+((tanh(beta*(w1(k)-w0-muR)/2.0d+00))/(tanh(beta*w0/2.0d+00))))&
            -ic*(lambda**2)*(gamaR/(gama**2+(w1(k)+w0-e0)**2))&
                           *(-1.0d+00+((tanh(beta*(w1(k)+w0-muR)/2.0d+00))/(tanh(beta*w0/2.0d+00))))
!-----------------------------------------------------------------------
!This is only valid at the particle-hole symmetric point:
!-----------------------------------------------------------------------
sKana(k)=ic*(lambda**2)*gama*(((-1.0d+00-tanh(beta*(w1(k)-w0)))/((w1(k)-w0)**2+gama**2))&
                             +((1.0d+00-tanh(beta*(w1(k)+w0)))/((w1(k)+w0)**2+gama**2)))
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Intana=2.0d+00*(-gama*(ap+am)*log((e0+SHart)**2+gama**2)&
                          +2*gama*(ap*log(abs(w1(k)+w0))+am*log(abs(w1(k)-w0)))&
                          +2*(ap*(w1(k)+w0-(e0+SHart))+am*(w1(k)-w0-(e0+SHart)))*atan2((e0+SHart),gama))
                          
Intana=2.0d+00*gama*pi*((-1.0d+00/(gama**2+(w0-w1(k))**2))+(1.0d+00/(gama**2+(w0+w1(k))**2)))                          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Spectral Function:
!-----------------------------------------------------------------------
rananum(k)=(-1.0d+00/pi)*(-gama+(aimag(sRana(k))))/((w1(k)-real(sRana(k)))**2+(gama-(aimag(sRana(k))))**2)
rFRG(k)=(-1.0d+00/pi)*(-gama+(YFRG(k)))/((w1(k)-YFRG((NEQ/4)+k))**2+(gama-(YFRG(k)))**2)
!-----------------------------------------------------------------------
!NdisFRG(k)=(1.0d+00/2.0d+00)-((YFRG(k+(2*NEQ/4)))/(4.0d+00*YFRG(k)))
NdisFRG(k)=(0.5d+00)-((YFRG(k+(2*NEQ/4))-(2.0d+00*gama*(1.0d+00-2.0d+00*fermi(w1(k)))))/(4.0d+00*(-gama+YFRG(k))))
!-----------------------------------------------------------------------
!write(1,*)ivoltage,(muL-muR)/gama,k,w1(k)/w0,(e0-ep)/gama,pi*gama*rananum(k),pi*gama*rFRG(k)&
!                  ,YFRG(k)/gama,YFRG(k+(NEQ/4))/ep,YFRG(k+(2*NEQ/4))/gama,YFRG(k+(3*NEQ/4))/ep&
!                  ,NdisFRG(k)&
!                  ,aimag(sRana(k))/gama,real(sRana(k))/ep&
!                  ,2.0d+00*YFRG(k)*(1.0d+00-2.0d+00*fermi(w1(k)))/gama,voltage/gama,deltaTem/gama&
!                  ,ilambda,ep/w0,fermi(w1(k)),(lambda/w0)**2
!print*,(muL-muR)/gama

!write(3,*)w1(k)/gama,((1.0d+00-tanh(betaL*(w1(k)-muL)/2.0d+00))/2.0d+00)-((1.0d+00-tanh(betaR*(w1(k)-muR)/2.0d+00))/2.0d+00)&
!         ,pi*gama*rFRG(k),((1.0d+00-tanh(betaL*(w1(k)-muL)/2.0d+00))/2.0d+00),((1.0d+00-tanh(betaR*(w1(k)-muR)/2.0d+00))/2.0d+00)

!----------------------------------------------------------------------

!write(1,*)k,w1(k)/w0,(e0-ep)/gama,real(SRana(k))/gama,aimag(SRana(k))/gama,real(SRananon(k))/gama,aimag(SRananon(k))/gama&
!         ,real(SKananon(k))/gama,aimag(SKananon(k))/gama&
!         ,YFRG(k)/gama,YFRG(k+(NEQ/4))/gama,YFRG(k+(2*NEQ/4))/gama,YFRG(k+(3*NEQ/4))/gama&
!         ,aimag(SKananon(k)-(1.0d+00-2.0d+00*fermi(w1(k)))*(2.0d+00*ic*aimag(SRananon(k))))&
!         ,(aimag(SKananon(k)-(1.0d+00-2.0d+00*fermi(w1(k)))*(2.0d+00*ic*aimag(SRananon(k)))))&
!           /((w1(k)-e0-real(SKananon(k)))**2+(gama-aimag(SKananon(k)))**2)&
!         ,(aimag(SKananon(k)-(1.0d+00-2.0d+00*fermi(w1(k)))*(2.0d+00*ic*aimag(SRananon(k)))))&
!           /((w1(k)-e0-0.0d+00*real(SKananon(k)))**2+(gama-0.0d+00*aimag(SKananon(k)))**2)
           
!----------------------------------------------------------------------
nespec(iparder,k)=rFRG(k)
neSigmaRim(iparder,k)=YFRG(k)
neSigmaRre(iparder,k)=YFRG((NEQ/4)+k)
neSigmaK(iparder,k)=YFRG(2*(NEQ/4)+k)
!----------------------------------------------------------------------
! checking:
!----------------------------------------------------------------------
! YFRGcommon(k)=aimag(sRana(k))        !0.0d+00!
! YFRGcommon(iw1max+k)=real(sRana(k))  !0.0d+00!
!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Going to the imaginary axis:
!----------------------------------------------------------------------
!----------------------------------------------------------------------
MF=10
!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ+NEQ**2 
LIW=20+NEQ
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ+(2*ML+MU)*NEQ
LIW=20+NEQ
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=20000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-13.0d+00)!(-16.0d+00)   !
ATOL=w0*10.0d+00**(-13.0d+00)!(-16.0d+00)   !
ISTATE=1
y=-10000.0d+00*w0
yout=10000.0d+00*w0
ml=1
mu=1
!------------------------------------------------------------------
! Initial Conditions:
!------------------------------------------------------------------
Ymat(1:iw1max)=0.0d+00                    !Imaginary Part of Greens Function
Ymat(iw1max+1:2*iw1max)=0.0d+00           !Real Part of Greens Function
!------------------------------------------------------------------
!------------------------------------------------------------------
! CALL Fint (2*iw1max, y, Ymat, YDmat)

!CALL JAC (NEQ, x, Yfrg, ml, mu, pd, nrowpd)

! CALL DLSODE (Fint, 2*iw1max, Ymat, y, yout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

!print*,"Fint",ISTATE,yout/w0,RTOL/gama

deallocate(RWORK,IWORK)
!-----------------------------------------------------------------------
Do k=1,iw1max

!-----------------------------------------------------------------------------------------------------------------------------
! Analytical results for perturbation theory in matsubara space: (+self-consistent Hartree term)
!-----------------------------------------------------------------------------------------------------------------------------

dp=-1.0d+00/(e0+SHart+w0-ic*(w1(k)+gama))
dm=-1.0d+00/(e0+SHart-w0-ic*(w1(k)+gama))
dtp=-1.0d+00/(e0+SHart+w0-ic*(w1(k)-gama))
dtm=-1.0d+00/(e0+SHart-w0-ic*(w1(k)-gama))        

ana(k)=((ic*(lambda**2)/(2.0d+00*pi))*(0.5d+00*(dtp-dtm-dp+dm)*log(gama**2+(e0+SHart)**2)&
                             -0.5d+00*(dtp-dtm-dp+dm)*log(w0**2+w1(k)**2)&
                             +ic*(dp-dtp)*atan2(w0,w1(k))+ic*(dtm-dm)*atan2(-w0,w1(k))&
                             -ic*(dp+dm)*pi*sgnn(w0)+ic*(dtp-dtm)*atan2(e0+SHart,-gama)&
                             -ic*(dp-dm)*atan2(e0+SHart,gama)&
                             -ic*pi*(dtp-dtm)*sgnn(e0+SHart)))&
                             +e0+SHart
                   
!---------------------------------------------------------------------- 

!Ymat(k)=(aimag(1.0d+00/(ic*(w1(k)+gama*tanh(beta*w1(k)/2.0d+00)))))
!Ymat(k+iw1max)=(real(1.0d+00/(ic*(w1(k)+SIGN(gama,w1(k))))))

self=ic*(w1(k)+SIGN(gama,w1(k)))-e0-(1.0d+00/(ic*Ymat(k)+Ymat(k+iw1max)))

! write(2,*)k,w1(k)/w0,aimag(self)/gama,real(self)/ep,Ymat(k)/gama,Ymat(iw1max+k)/ep&
!             ,aimag(ana(k))/gama,real(ana(k))/ep,YFRGcommon(k)/gama,YFRGcommon(iw1max+k)/ep&
!             ,YFRG(k)/gama,YFRG(iw1max+k)/ep,aimag(sRana(k))/gama,real(sRana(k))/ep&
!             ,YFRG(k+2*iw1max),YFRG(k+3*iw1max),(tanh(beta*w1(k)/2.0d+00))*(2.0d+00*YFRG(k))

end do
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Charge and Heat Currents:
!----------------------------------------------------------------------
!----------------------------------------------------------------------
MF=10
NEQ2=4
!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ2
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ2+NEQ2**2 
LIW=20+NEQ2
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ2+(2*ML+MU)*NEQ2
LIW=20+NEQ2
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=20000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-15.0d+00)!(-16.0d+00)   !
ATOL=w0*10.0d+00**(-15.0d+00)!(-16.0d+00)   !
ISTATE=1
y=min(muR,muL)-(30.0d+00/betaL)
yout=max(muR,muL)+(30.0d+00/betaL)!10.0d+00*w0
ml=1
mu=1
!------------------------------------------------------------------
! Initial Conditions:
!------------------------------------------------------------------
Ycur(1:4)=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL Fcurrent (4, y, Ycur, YDcur)

CALL DLSODE (Fcurrent, 4, Ycur, y, yout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"Fcur",ISTATE,yout/w0,RTOL/gama,Ycur(1)/(gama),Ycur(2)/(gama**2),Ycur(4)&
             ,vol(ivoltage)/gama,deltaTem/gama,betaL*gama,betaR/gama

print*,"violation of current conservation:",Ycur(4)/gama

deallocate(RWORK,IWORK)
!-----------------------------------------------------------------------
Jch(iparder)=Ycur(1)
JQR(iparder)=Ycur(2)+((voltage/2.0d+00)*Ycur(1))
JQL(iparder)=Ycur(3)+((voltage/2.0d+00)*Ycur(1))
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Moments of nonequilibrium spectral function:
!----------------------------------------------------------------------
if(iparder==1)then! .and. ipardermax>1)then
!----------------------------------------------------------------------
NEQ4=3
MF=10
!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ4
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ4+NEQ4**2 
LIW=20+NEQ4
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ4+(2*ML+MU)*NEQ4
LIW=20+NEQ4
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=20000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-15.0d+00)!(-16.0d+00)   !
ATOL=w0*10.0d+00**(-15.0d+00)!(-16.0d+00)   !
ISTATE=1
y=(1.0d+00)*0.5d+00*(vol(ivoltage))-(40.0d+00/betaL)
yout=(1.0d+00)*0.5d+00*(vol(ivoltage))+(40.0d+00/betaL)
ml=1
mu=1
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initial Conditions:
!------------------------------------------------------------------
YInL(1:3)=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL FInL (3, y, YInL, YDInL)

CALL DLSODE (FInL, 3, YInL, y, yout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"FInL",ISTATE,yout/w0,RTOL/gama,YInL(1),YInL(2),YInL(3)
!-----------------------------------------------------------------------
deallocate(RWORK,IWORK)
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
NEQ4=4
MF=10
!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ4
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ4+NEQ4**2 
LIW=20+NEQ4
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ4+(2*ML+MU)*NEQ4
LIW=20+NEQ4
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=20000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-16.0d+00)!(-16.0d+00)   !
ATOL=w0*10.0d+00**(-16.0d+00)!(-16.0d+00)   !
ISTATE=1
y=-0.5d+00*(vol(ivoltage))-(40.0d+00/betaR)
yout=-0.5d+00*(vol(ivoltage))+(40.0d+00/betaR)
ml=1
mu=1
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initial Conditions:
!------------------------------------------------------------------
YInR(1:4)=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL FInR (4, y, YInR, YDInR)

CALL DLSODE (FInR, 4, YInR, y, yout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"FInR",ISTATE,yout/w0,RTOL/gama,YInR(1),YInR(2),YInR(3),YInR(4)
!-----------------------------------------------------------------------
deallocate(RWORK,IWORK)
!-----------------------------------------------------------------------
elseif(iparder==2)then
!----------------------------------------------------------------------
dnespec(1:iw1max)=(nespec(2,1:iw1max)-nespec(1,1:iw1max))/deltavoltage

dneSigmaRre(1:iw1max)=(neSigmaRre(2,1:iw1max)-neSigmaRre(1,1:iw1max))/deltavoltage
dneSigmaRim(1:iw1max)=(neSigmaRim(2,1:iw1max)-neSigmaRim(1,1:iw1max))/deltavoltage
dneSigmaK(1:iw1max)=(neSigmaK(2,1:iw1max)-neSigmaK(1,1:iw1max))/deltavoltage

NEQ4=3
MF=10
!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ4
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ4+NEQ4**2 
LIW=20+NEQ4
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ4+(2*ML+MU)*NEQ4
LIW=20+NEQ4
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=20000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-16.0d+00)!(-16.0d+00)   !
ATOL=w0*10.0d+00**(-16.0d+00)!(-16.0d+00)   !
ISTATE=1
y=-0.5d+00*(vol(ivoltage))-(40.0d+00/betaL)
yout=0.5d+00*(vol(ivoltage))+(40.0d+00/betaL)
ml=1
mu=1
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initial Conditions:
!------------------------------------------------------------------
YIncor(1:3)=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL FIncor (3, y, YIncor, YDIncor)

CALL DLSODE (FIncor, 3, YIncor, y, yout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"FIncorL-Voltage",ISTATE,yout/w0,RTOL/gama,YIncor(1),YIncor(2),YIncor(3)

YIncorvolt(1)=YIncor(1)
YIncorvolt(2)=YIncor(2) !#L
YIncorvolt(3)=YIncor(3) !#R
!-----------------------------------------------------------------------
deallocate(RWORK,IWORK)
!-----------------------------------------------------------------------
NEQ4=2
MF=10
!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ4
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ4+NEQ4**2 
LIW=20+NEQ4
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ4+(2*ML+MU)*NEQ4
LIW=20+NEQ4
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=20000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-16.0d+00)!(-16.0d+00)   !
ATOL=w0*10.0d+00**(-16.0d+00)!(-16.0d+00)   !
ISTATE=1
y=-0.5d+00*(vol(ivoltage))-(40.0d+00/betaL)
yout=0.5d+00*(vol(ivoltage))+(40.0d+00/betaL)
ml=1
mu=1
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initial Conditions:
!------------------------------------------------------------------
YInnedis(1:2)=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL FInnedis(2, y, YInnedis, YDInnedis)

CALL DLSODE (FInnedis, 2, YInnedis, y, yout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"FInnedis-Voltage",ISTATE,yout/w0,RTOL/gama,YInnedis(1),YInnedis(2)

YInnedisvolt(1)=YInnedis(1)
YInnedisvolt(2)=YInnedis(2)
!-----------------------------------------------------------------------
deallocate(RWORK,IWORK)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
elseif(iparder==3)then
!----------------------------------------------------------------------
dnespec(1:iw1max)=(nespec(3,1:iw1max)-nespec(1,1:iw1max))/deltadeltaTem

dneSigmaRre(1:iw1max)=(neSigmaRre(3,1:iw1max)-neSigmaRre(1,1:iw1max))/deltadeltaTem
dneSigmaRim(1:iw1max)=(neSigmaRim(3,1:iw1max)-neSigmaRim(1,1:iw1max))/deltadeltaTem
dneSigmaK(1:iw1max)=(neSigmaK(3,1:iw1max)-neSigmaK(1,1:iw1max))/deltadeltaTem


NEQ4=3
MF=10
!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ4
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ4+NEQ4**2 
LIW=20+NEQ4
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ4+(2*ML+MU)*NEQ4
LIW=20+NEQ4
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=20000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-16.0d+00)!(-16.0d+00)   !
ATOL=w0*10.0d+00**(-16.0d+00)!(-16.0d+00)   !
ISTATE=1
y=-0.5d+00*(vol(ivoltage))-(40.0d+00/betaL)
yout=0.5d+00*(vol(ivoltage))+(40.0d+00/betaL)
ml=1
mu=1
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initial Conditions:
!------------------------------------------------------------------
YIncor(1:3)=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL FIncor (3, y, YIncor, YDIncor)

CALL DLSODE (FIncor, 3, YIncor, y, yout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"FIncor-deltaT",ISTATE,yout/w0,RTOL/gama,YIncor(1),YIncor(2),YIncor(3)

YIncordeltaT(1)=YIncor(1)
YIncordeltaT(2)=YIncor(2)
YIncordeltaT(3)=YIncor(3)
!-----------------------------------------------------------------------
deallocate(RWORK,IWORK)
!-----------------------------------------------------------------------
NEQ4=2
MF=10
!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ4
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ4+NEQ4**2 
LIW=20+NEQ4
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ4+(2*ML+MU)*NEQ4
LIW=20+NEQ4
endif
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1
RWORK(5:10)=0.0d+00
IWORK(5:10)=0
IWORK(6)=20000
!------------------------------------------------------------------
RTOL=w0*10.0d+00**(-16.0d+00)!(-16.0d+00)   !
ATOL=w0*10.0d+00**(-16.0d+00)!(-16.0d+00)   !
ISTATE=1
y=-0.5d+00*(vol(ivoltage))-(40.0d+00/betaL)
yout=0.5d+00*(vol(ivoltage))+(40.0d+00/betaL)
ml=1
mu=1
!------------------------------------------------------------------
!------------------------------------------------------------------
! Initial Conditions:
!------------------------------------------------------------------
YInnedis(1:2)=0.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
CALL FInnedis(2, y, YInnedis, YDInnedis)

CALL DLSODE (FInnedis, 2, YInnedis, y, yout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)

print*,"FInnedis-deltaT",ISTATE,yout/w0,RTOL/gama,YInnedis(1),YInnedis(2)


YInnedisdeltaT(1)=YInnedis(1)
YInnedisdeltaT(2)=YInnedis(2)
!-----------------------------------------------------------------------
deallocate(RWORK,IWORK)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
endif
!-----------------------------------------------------------------------
end do !iparder
!-----------------------------------------------------------------------
Jcharge=Jch(1)
JheatR=JQR(1)
JheatL=JQL(1)
!-----------------------------------------------------------------------
! Linear response:
!-----------------------------------------------------------------------
!if(ivoltage==1)then

dJchdv=(Jch(2)-Jch(1))/deltavoltage
dJQdv=(JQL(2)-JQL(1))/deltavoltage

dJchddeltaT=(Jch(3)-Jch(1))/(deltadeltaTem)
dJQddeltaT=(JQL(3)-JQL(1))/(deltadeltaTem)

Gdif=dJchdv
Sdif=dJchddeltaT/dJchdv
Kedif=-(-dJQddeltaT+((dJQdv)*(dJchddeltaT/dJchdv)))

PFdif=(Sdif**2)*Gdif
ZTdif=PFdif/(beta*Kedif)
Ldif=beta*kedif/Gdif

!-----------------------------------------------------------------------
! Results from the decomposition:
!-----------------------------------------------------------------------
dJchdvdec=0.0d+00*YIncorvolt(1)+((YInL(1)+YInR(1))/2.0d+00)
dJQRdvdec=(Jch(1)/2.0d+00)-(muR*(dJchdvdec))+(((0.5d+00*YInR(2))+YInnedisvolt(2)-YIncorvolt(3)+YInnedisvolt(1))/(gamaL/gama))
dJQLdvdec=(Jch(1)/2.0d+00)+(muL*(dJchdvdec))-(((0.5d+00*YInL(2))-YInnedisvolt(2)+YIncorvolt(2)-YInnedisvolt(1))/(gamaR/gama))

dJchddeltaTdec=YIncordeltaT(1)+(((-betaR*(YInR(2)-((muR)*YInR(1))))))
dJQLddeltaTdec=(muL*(dJchddeltaTdec))+(((YInnedisdeltaT(2))-(YIncordeltaT(2)-YInnedisdeltaT(1)))/(gamaR/gama))
dJQRddeltaTdec=-(muR*(dJchddeltaTdec))-(((betaR*(YInR(3)-(muR*YInR(2))))&
                                      -(YInnedisdeltaT(2))+(YIncordeltaT(3)-YInnedisdeltaT(1)))/(gamaL/gama))

Gdifdec=dJchdvdec
Sdifdec=dJchddeltaTdec/dJchdvdec
Kedifdec=(-dJQRddeltaTdec+((dJQRdvdec)*(dJchddeltaTdec/dJchdvdec)))

PFdifdec=(Sdifdec**2)*Gdifdec
ZTdifdec=PFdifdec/(beta*Kedifdec)
Ldifdec=beta*kedifdec/Gdifdec

if(Gdifdec<0.5d+00)then
slogGdifdec=2.0d+00*Gdifdec
else
slogGdifdec=log10(20.0d+00*Gdifdec)
endif


print*,"Gdif",Gdif,Gdifdec
print*,"Sdif",Sdif,Sdifdec
print*,"Kedif",Kedif,Kedifdec

print*,"dJchdv",dJchdv,Gdifdec
print*,"dJQRdv",dJQdv,dJQRdvdec,dJQLdvdec
print*,"dJchddeltaT",dJchddeltaT,dJchddeltaTdec
print*,"dJQddeltaT",dJQddeltaT,dJQRddeltaTdec,dJQLddeltaTdec

!endif
!-----------------------------------------------------------------------
powerlinear=(voltage**2)*Gdif*(-1.0d+00+Sdif*(-deltaTem/voltage))
powerlineardec=(voltage**2)*Gdifdec*(-1.0d+00+Sdifdec*(-deltaTem/voltage))
power=-Jcharge*voltage

etaZT=(ZTdif/(ZTdif+2.0d+00))*0.5d+00 !*(0.1d+00*gama*beta)
etaZTdec=(ZTdifdec/(ZTdifdec+2.0d+00))*0.5d+00 !*(0.1d+00*gama*beta)
etaCA=1.0d+00-(betaR/betaL)

if(Jcharge<0.0d+00 .and. JheatR<0.0d+00)then

eta=(Jcharge*(muL-muR)/JheatR)

else

eta=-0.01d+00

endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
fanoana=-(G0/(2.0d+00*gama))*((1.0d+00/((((muL-muR)/(2.0d+00*gama))-((e0-ep)/gama))**3))&
                             +(1.0d+00/((((muL-muR)/(2.0d+00*gama))+((e0-ep)/gama))**3)))
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!if(Jcharge<0.0d+00 .and. JheatR<0.0d+00)then

write(2,*)ivoltage,(muL-muR)/gama,(e0-ep)/w0,1.0d+00/(beta*gama)&
            ,Jcharge/(gama*G0),JheatR/((gama**2)*G0),JheatL/((gama**2)*G0),Ycur(4),(lambda/w0)**2&
!            ,Jcharge*38.74044586573/(gama*G0),e0/w0&
!            ,atan((voltage-2.0d+00*(e0-ep))/(2.0d+00*gama))+atan((voltage+2.0d+00*(e0-ep))/(2.0d+00*gama))&
!            ,Gdifdec/G0,Gdif/G0,Jch(2)/gama&
!            ,0.5d+00*((1.0d+00/(1.0d+00+((voltage-2.0d+00*(e0-ep))/(2.0d+00*gama))**2))&
!            +(1.0d+00/(1.0d+00+((voltage+2.0d+00*(e0-ep))/(2.0d+00*gama))**2))),slogGdifdec!&            
             ,power/(gama**2),powerlinear/(gama**2),powerlineardec/(gama**2),eta/etaCA,eta,etaZT,etaZTdec,deltaTem/gama,betaL/betaR&
            ,Gdif/G0,Sdif,Kedif/gama,PFdif,ZTdif,Ldif&
            ,ilambda,ep/w0,deltavoltage/gama&
            ,Gdifdec/G0,Sdifdec,Kedifdec/gama,PFdifdec,ZTdifdec,Ldifdec&               
            ,Jch(1)/gama,Jch(2)/gama,Jch(3)/gama,JQL(1)/(gama**2),JQL(2)/(gama**2),JQL(3)/(gama**2)&
            ,deltavoltage/gama,deltadeltaTem/gama&
            ,dJchdv,dJchddeltaT,dJQdv,dJQddeltaT,dJchdvdec,dJchddeltaTdec&
            ,dJQRdvdec,dJQRddeltaTdec,dJQLdvdec!,dJQLddeltaTdec

!endif
! write(2,*)ivoltage,(muL-muR)/gama,(e0-ep)/gama,1.0d+00/(beta*gama),fanoana

print*,"Temperatures",betaL,betaR,deltaTem*betaR
!-----------------------------------------------------------------------
end do !ivoltage
!-----------------------------------------------------------------------
end do !ideltaTem
!-----------------------------------------------------------------------
end do !ibeta
!-----------------------------------------------------------------------
end do !ie0int
!-----------------------------------------------------------------------
end do !ilambda
!-----------------------------------------------------------------------
call cpu_time(stop_time)
print*,"total time",stop_time - start_time
!----------------------------------------------------------------------
close(1)
!----------------------------------------------------------------------
end program KeldyshFRG_SRLMCP
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!*****************FRG:*************************************************
!----------------------------------------------------------------------
subroutine Ffrg(NEQ, x, Yfrg, YDfrg)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon

INTEGER::k,NEQ,i
INTEGER,parameter::NEQ1=1,ITOLin=1
DOUBLE PRECISION::x,xcommon,w,wout
DOUBLE PRECISION,DIMENSION(800)::YFRG,YDFRG,SFRG,SDFRG,YFRGcommon
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,beta,fcoef,betaL,betaR,gamaR,gamaL&
                  ,muL,muR,fermi,fermiaux,bose,betaph,betaaux
COMPLEX*16::ic,GLRFRG,GLAFRG,GLKFRG,IntK,IntRA,RSelfinterp,ASelfinterp,KSelfinterp
COMPLEX*16::SLKplusw0,SLKminusw0,SLRplusw0,SLRminusw0,SLAplusw0,SLAminusw0,Intana&
            ,IntKpmw0,IntRApmw0
DOUBLE PRECISION::IntCaushy,ap,am,SHart
DOUBLE PRECISION,DIMENSION(200)::w1
!----------------------------------------------------------------------
! Parameters for Quadpack:
!----------------------------------------------------------------------
DOUBLE PRECISION::wmin,wmax,wpole,abserr,epsabs,epsrel,resultquad,result
INTEGER,PARAMETER::limit=400,lenw=4*limit
INTEGER,DIMENSION(limit)::Iworkquad,iord
INTEGER::last,neval,ier
DOUBLE PRECISION::workquad(lenw)
DOUBLE PRECISION, EXTERNAL::fKimcaushy,fRArecaushy,fKrecaushy,fRAimcaushy,funSRA,funSK
!----------------------------------------------------------------------
INTEGER,PARAMETER::pointsdim=6,limitdis=400,leniwdis=2*limitdis+pointsdim,lenwdis=2*leniwdis-pointsdim
DOUBLE PRECISION,DIMENSION(limit)::alist,blist,rlist,elist
DOUBLE PRECISION,DIMENSION(pointsdim)::points
INTEGER,DIMENSION(leniwdis)::Iworkquaddis
DOUBLE PRECISION::workquaddis(lenwdis)
!------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
integer*4, parameter::DIM_P=4          !the spatial dimension.
integer*4, parameter::DATA_NUM=200    !the number of data points
integer*4, parameter::INTERP_NUM=200     !the number of points at which interpolation is to be done

real*8::T_DATA(DATA_NUM)               !the value of the independent variable at the sample points
real*8::P_DATA1(DIM_P,DATA_NUM)         !the value of the dependent variables at the sample points
real*8::T_INTERP1(INTERP_NUM),T_INTERP2(INTERP_NUM)           !the value of the independent variable at the interpolation points
real*8:: P_INTERP1(DIM_P,INTERP_NUM),P_INTERP2(DIM_P,INTERP_NUM)      !the interpolated values of the dependent variables
!----------------------------------------------------------------------
! Parameters for ODEPACK library:
!----------------------------------------------------------------------
INTEGER:: Flagin, ISTATEin, IOPTin, LRWin, LIWin, MFin,MUin, MLin
DOUBLE PRECISION,DIMENSION(NEQ1)::YSRA,YDSRA
DOUBLE PRECISION::RTOLin, ATOLin,JACin
DOUBLE PRECISION,allocatable,dimension(:)::RWORKin
INTEGER,allocatable,dimension(:)::IWORKin
!------------------------------------------------------------------
!----------------------------------------------------------------------
ep=(lambda**2)/w0
SHart=-ep
ic=cmplx(0.0d+00,1.0d+00)
!---------------------------------------------------------------------- 
IntK=cmplx(0.0d+00,0.0d+00)
IntRA=cmplx(0.0d+00,0.0d+00)
YFRGcommon=YFRG
xcommon=x

!print*,x/w0

Do k=1,NEQ/4

GLRFRG=1.0d+00/(w1(k)+ic*(x+gama)-fcoef*ic*YFRG(k)-fcoef*YFRG((NEQ/4)+k)-(1.0d+00-fcoef)*e0)
GLAFRG=CONJG(GLRFRG)
GLKFRG=GLRFRG*GLAFRG*(fcoef*ic*YFRG((2*NEQ/4)+k)+fcoef*YFRG((3*NEQ/4)+k)-2.0d+00*ic*gama*(1.0d+00-2.0d+00*fermi(w1(k)))&
                                                            -2.0d+00*ic*x*((1.0d+00-2.0d+00*fermiaux(w1(k)))))
! Retarded
SFRG(k)=aimag(ic*((GLRFRG)**2))
SFRG((NEQ/4)+k)=real(ic*((GLRFRG)**2))

! Keldysh
SFRG((2*NEQ/4)+k)=aimag(ic*GLRFRG*GLKFRG-ic*GLAFRG*GLKFRG+GLRFRG*GLAFRG*(2.0d+00*ic*(1.0d+00-2.0d+00*fermiaux(w1(k)))))
SFRG((3*NEQ/4)+k)=real(ic*GLRFRG*GLKFRG-ic*GLAFRG*GLKFRG+GLRFRG*GLAFRG*(2.0d+00*ic*(1.0d+00-2.0d+00*fermiaux(w1(k)))))


!if(w1(k)>0)then

!IntK=IntK+1.0d+00*(ic*(((SFRG((2*NEQ/6)+k)+SFRG((2*NEQ/6)+k+1))/2.0d+00))&
!                      +((SFRG((3*NEQ/6)+k)+SFRG((3*NEQ/6)+k+1))/2.0d+00))*abs(w1(k+1)-w1(k))

!IntRA=IntRA+1.0d+00*(ic*(SFRG(k))+SFRG((1*NEQ/6)+k)+ic*(SFRG((4*NEQ/6)+k))+SFRG((5*NEQ/6)+k))*abs(w1(k+1)-w1(k))

!IntRA=IntRA+2.0d+00*(((SFRG((1*NEQ/6)+k)+SFRG((1*NEQ/6)+k+1))/2.0d+00))*abs(w1(k+1)-w1(k))

!else
!IntK=IntK+1.0d+00*(ic*(SFRG((2*NEQ/6)+k))+SFRG((3*NEQ/6)+k))*abs(w1(k)-w1(k-1))
!IntRA=IntRA+1.0d+00*(ic*(SFRG(k))+SFRG((1*NEQ/6)+k)+ic*(SFRG((4*NEQ/6)+k))+SFRG((5*NEQ/6)+k))*abs(w1(k)-w1(k-1))

!endif


T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA1(1,k)=YFRG((DATA_NUM+1)-k)  
P_DATA1(2,k)=YFRG((1*DATA_NUM)+(DATA_NUM+1)-k)
P_DATA1(3,k)=YFRG((2*DATA_NUM)+(DATA_NUM+1)-k)
P_DATA1(4,k)=YFRG((3*DATA_NUM)+(DATA_NUM+1)-k)


T_INTERP1(k)=w1(k)-w0
T_INTERP2(k)=w1(k)+w0
!----------------------------------------------------------------------
!IntRA=cmplx(0.0d+00,0.0d+00)
!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------
call interp_linear(DIM_P,NEQ/4,T_DATA,P_DATA1,INTERP_NUM,T_INTERP1,P_INTERP1)
call interp_linear(DIM_P,NEQ/4,T_DATA,P_DATA1,INTERP_NUM,T_INTERP2,P_INTERP2)
!----------------------------------------------------------------------
! write(1,*)x/w0,w1(24)/w0,T_DATA(24)/w0,P_INTERP1(2,24)/gama,YFRGcommon(300+24)/gama
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!------------------------------------------------------------------
! Sloving Integrals: IntRA
!------------------------------------------------------------------
!------------------------------------------------------------------
MFin=10
!------------------------------------------------------------------
! Dimension Declaration of WORK:
!------------------------------------------------------------------
If(MFin==10)then
LRWin=20+16*NEQ1
LIWin=20
endif
allocate(RWORKin(LRWin),IWORKin(LIWin))
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
Flagin=1
!------------------------------------------------------------------
IOPTin=1
RWORKin(5:10)=0.0d+00
IWORKin(5:10)=0
IWORKin(6)=1000
!------------------------------------------------------------------
RTOLin=w0*10.0d+00**(-4.0d+00)!(-16.0d+00)   !
ATOLin=w0*10.0d+00**(-4.0d+00)!(-16.0d+00)   !
ISTATEin=1
w=20.0d+00*w0
wout=-20.0d+00*w0
mlin=1
muin=1

YSRA(1)=0.0d+00                !Imaginary Part of Retarded
!------------------------------------------------------------------
!------------------------------------------------------------------
!CALL FintSRA (NEQ1, w, YSRA, YDSRA)


!CALL DLSODE (FintSRA, NEQ1, YSRA, w, wout, ITOLin, RTOLin, ATOLin, Flagin,ISTATEin, IOPTin&
!             , RWORKin, LRWin, IWORKin, LIWin, JACin, MFin)

!print*,"intRA",ISTATEin,wout/w0,RTOLin/gama

!deallocate(RWORKin,IWORKin)
!IntRA=YSRA(1)
!----------------------------------------------------------------------
!epsabs=w0*(10.0d+00)**(-4.0d+00)
!epsrel=w0*(10.0d+00)**(-4.0d+00)
!CALL dqagie(funSRA,0.0d+00,2,epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
IntRA=0.0d+00!result

!if(ier>0)then
!print*,"ERROR"
!endif
!----------------------------------------------------------------------
epsabs=w0*(10.0d+00)**(-8.0d+00)
epsrel=w0*(10.0d+00)**(-8.0d+00)
!wmin=-100.0d+00*w0
!wmax=100.0d+00*w0
!wpole=w0

CALL dqagie(funSK,0.0d+00,2,epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

IntK=ic*result

if(ier>0)then
!print*,"ERROR",ier,neval,abserr/w0
endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
Do k=1,DATA_NUM
!----------------------------------------------------------------------
! Initialization:
!----------------------------------------------------------------------
IntKpmw0=cmplx(0.0d+00,0.0d+00)
IntRApmw0=cmplx(0.0d+00,0.0d+00)
!----------------------------------------------------------------------
! Interpolation at w1(k)-w0:
!----------------------------------------------------------------------
GLRFRG=1.0d+00/(w1(k)-w0+ic*(x+gama)-fcoef*ic*P_INTERP1(1,k)-fcoef*P_INTERP1(2,k)-(1.0d+00-fcoef)*e0)
GLAFRG=CONJG(GLRFRG)
GLKFRG=GLRFRG*GLAFRG*(fcoef*ic*P_INTERP1(3,k)+fcoef*P_INTERP1(4,k)-2.0d+00*ic*gama*(1.0d+00-2.0d+00*fermi(w1(k)-w0))&
                                                      -2.0d+00*ic*x*((1.0d+00-2.0d+00*fermiaux(w1(k)-w0))))

SLRminusw0=ic*((GLRFRG)**2)
SLAminusw0=CONJG(SLRminusw0)!-ic*((GLAFRG)**2)
SLKminusw0=ic*GLRFRG*GLKFRG-ic*GLAFRG*GLKFRG+GLRFRG*GLAFRG*(2.0d+00*ic*(1.0d+00-2.0d+00*fermiaux(w1(k)-w0)))

!SLKminusw0=(tanh(beta*(w1(k)-w0)/2.0d+00))*(SLRminusw0-SLAminusw0)

!----------------------------------------------------------------------
! Interpolation at w1(k)+w0:
!----------------------------------------------------------------------
GLRFRG=1.0d+00/(w1(k)+w0+ic*(x+gama)-fcoef*ic*P_INTERP2(1,k)-fcoef*P_INTERP2(2,k)-(1.0d+00-fcoef)*e0)
GLAFRG=CONJG(GLRFRG)
GLKFRG=GLRFRG*GLAFRG*(fcoef*ic*P_INTERP2(3,k)+fcoef*P_INTERP2(4,k)-2.0d+00*ic*gama*(1.0d+00-2.0d+00*fermi(w1(k)+w0))&
                                                      -2.0d+00*ic*x*((1.0d+00-2.0d+00*fermiaux(w1(k)+w0))))

SLRplusw0=ic*((GLRFRG)**2)
SLAplusw0=CONJG(SLRplusw0)!-ic*((GLAFRG)**2)
SLKplusw0=ic*GLRFRG*GLKFRG-ic*GLAFRG*GLKFRG+GLRFRG*GLAFRG*(2.0d+00*ic*(1.0d+00-2.0d+00*fermiaux(w1(k)+w0)))

!SLKplusw0=(tanh(beta*(w1(k)+w0)/2.0d+00))*(SLRplusw0-SLAplusw0)

!----------------------------------------------------------------------
! Caushy principal value: 
!----------------------------------------------------------------------
epsabs=w0*(10.0d+00)**(-8.0d+00)
epsrel=w0*(10.0d+00)**(-8.0d+00)
wmin=-20.0d+00*w0
wmax=20.0d+00*w0
wpole=w1(k)-w0

CALL dqawc(fKimcaushy,wmin,wmax,wpole,epsabs,epsrel,resultquad,abserr,neval,ier,limit,lenw,last,iworkquad,workquad)

IntKpmw0=-ic*resultquad

if(ier>0)then
!print*,"ERROR-coushy",ier,neval,abserr/w0,w1(k)/w0
endif

epsabs=w0*(10.0d+00)**(-8.0d+00)
epsrel=w0*(10.0d+00)**(-8.0d+00)
wmin=-20.0d+00*w0
wmax=20.0d+00*w0
wpole=w1(k)+w0

CALL dqawc(fKimcaushy,wmin,wmax,wpole,epsabs,epsrel,resultquad,abserr,neval,ier,limit,lenw,last,iworkquad,workquad)   

IntKpmw0=IntKpmw0+ic*resultquad

if(ier>0)then
!print*,"ERROR--coushy",ier,neval,abserr/w0,w1(k)/w0
endif
!-----------------------------------------------------------------------
epsabs=w0*(10.0d+00)**(-8.0d+00)
epsrel=w0*(10.0d+00)**(-8.0d+00)
wmin=-20.0d+00*w0 !w1(k)-w0
wmax=20.0d+00*w0  !w1(k)-w0+
wpole=w1(k)-w0

!CALL dqawc(fRArecaushy,wmin,wmax,wpole,epsabs,epsrel,resultquad,abserr,neval,ier,limit,lenw,last,iworkquad,workquad)

!IntRApmw0=-resultquad

epsabs=w0*(10.0d+00)**(-8.0d+00)
epsrel=w0*(10.0d+00)**(-8.0d+00)
wmin=-20.0d+00*w0  !w1(k)+w0
wmax=20.0d+00*w0   !w1(k)+w0+
wpole=w1(k)+w0

!CALL dqawc(fRArecaushy,wmin,wmax,wpole,epsabs,epsrel,resultquad,abserr,neval,ier,limit,lenw,last,iworkquad,workquad)   

!IntRApmw0=IntRApmw0+resultquad

!write(2,*)k,(w1(k)/w0),aimag(IntKpmw0)/gama
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Checking:
!IntRApmw0=cmplx(0.0d+00,0.0d+00)
!IntKpmw0=cmplx(0.0d+00,0.0d+00)
!IntK=0.0d+00
!IntRA=0.0d+00
!----------------------------------------------------------------------
! Flow Equations:
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Retarded
!----------------------------------------------------------------------
YDFRG(k)=aimag(-ic*(ep/(2.0d+00*pi))*IntK-ic*((lambda**2)/(4.0d+00*pi))*IntKpmw0&
         -1.0d+00*((lambda**2)/4.0d+00)*(SLKminusw0-SLKplusw0)&
         -1.0d+00*((lambda**2)/2.0d+00)*(SLRminusw0+SLRplusw0)*(1.0d+00+2.0d+00*bose(w0)))

YDFRG((NEQ/4)+k)=real(-ic*(ep/(2.0d+00*pi))*IntK-ic*((lambda**2)/(4.0d+00*pi))*IntKpmw0&
         -1.0d+00*((lambda**2)/4.0d+00)*(SLKminusw0-SLKplusw0)&
         -1.0d+00*((lambda**2)/2.0d+00)*(SLRminusw0+SLRplusw0)*(1.0d+00+2.0d+00*bose(w0)))         
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!Keldysh  
!----------------------------------------------------------------------
YDFRG((2*NEQ/4)+k)=aimag(-1.0d+00*((lambda**2)/2.0d+00)*(SLRminusw0-SLRplusw0-SLAminusw0+SLAplusw0)&
                         -1.0d+00*((lambda**2)/2.0d+00)*(SLKminusw0+SLKplusw0)*(1.0d+00+2.0d+00*bose(w0)))

YDFRG((3*NEQ/4)+k)=real(-1.0d+00*((lambda**2)/2.0d+00)*(SLRminusw0-SLRplusw0-SLAminusw0+SLAplusw0)&
                        -1.0d+00*((lambda**2)/2.0d+00)*(SLKminusw0+SLKplusw0)*(1.0d+00+2.0d+00*bose(w0)))


!----------------------------------------------------------------------               
!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------  
end subroutine Ffrg
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!-------------------------------------------------------------
!*************************************************************
DOUBLE PRECISION function fKimcaushy(w2)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon

DOUBLE PRECISION,intent(in)::w2

!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
INTEGER*4, parameter::DIM_P=4          !the spatial dimension.
integer*4, parameter::DATA_NUM=200    !the number of data points
INTEGER*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

DOUBLE PRECISION::T_DATA(DATA_NUM)               !the value of the independent variable at the sample points
DOUBLE PRECISION::P_DATA1(DIM_P,DATA_NUM)         !the value of the dependent variables at the sample points
DOUBLE PRECISION::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
DOUBLE PRECISION:: P_INTERP1(DIM_P,INTERP_NUM)      !the interpolated values of the dependent variables
DOUBLE PRECISION::ep,lambda,w0,gama,e0,pi,beta,x,xcommon,fcoef&
                  ,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,fermi,fermiaux
DOUBLE PRECISION,DIMENSION(200)::w1
INTEGER::NEQ,k
DOUBLE PRECISION,DIMENSION(800)::YFRGcommon,YDFRG,SFRG,SDFRG
COMPLEX*16::ic,GLRFRG,GLAFRG,GLKFRG,RSelfinterp,ASelfinterp,KSelfinterp
COMPLEX*16::SLRinterp,SLAinterp,SLKinterp
!----------------------------------------------------------------------
!----------------------------------------------------------------------
ic=cmplx(0.0d+00,1.0d+00)
ep=(lambda**2)/w0

!print*,NEQ
!Do k=1,NEQ
Do k=1,DATA_NUM

T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA1(1,k)=YFRGcommon((DATA_NUM+1)-k)  
P_DATA1(2,k)=YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k)
P_DATA1(3,k)=YFRGcommon((2*DATA_NUM)+(DATA_NUM+1)-k)
P_DATA1(4,k)=YFRGcommon((3*DATA_NUM)+(DATA_NUM+1)-k)
!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------
T_INTERP(INTERP_NUM)=w2
call interp_linear(DIM_P,DATA_NUM,T_DATA,P_DATA1,INTERP_NUM,T_INTERP,P_INTERP1)
!----------------------------------------------------------------------
RSelfinterp=ic*P_INTERP1(1,1)+P_INTERP1(2,1)
ASelfinterp=CONJG(RSelfinterp)
KSelfinterp=ic*P_INTERP1(3,1)+P_INTERP1(4,1)

GLRFRG=1.0d+00/(w2+ic*(xcommon+gama)-fcoef*RSelfinterp-(1.0d+00-fcoef)*e0)
GLAFRG=1.0d+00/(w2-ic*(xcommon+gama)-fcoef*ASelfinterp-(1.0d+00-fcoef)*e0)
GLKFRG=GLRFRG*GLAFRG*(fcoef*KSelfinterp-2.0d+00*ic*gama*(1.0d+00-2.0d+00*fermi(w2))&
                       -2.0d+00*ic*xcommon*(1.0d+00-2.0d+00*fermiaux(w2)))

SLRinterp=ic*((GLRFRG)**2)
SLAinterp=CONJG(SLRinterp)!-ic*((GLAFRG)**2)
SLKinterp=ic*GLRFRG*GLKFRG-ic*GLAFRG*GLKFRG+GLRFRG*GLAFRG*(2.0d+00*ic*(1.0d+00-2.0d+00*fermiaux(w2)))


fKimcaushy=aimag(SLKinterp)


end function fKimcaushy
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***********************************************************************
DOUBLE PRECISION function fRArecaushy(w2)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon

DOUBLE PRECISION,INTENT(in)::w2

!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
INTEGER*4, parameter::DIM_P=4          !the spatial dimension.
integer*4, parameter::DATA_NUM=200    !the number of data points
INTEGER*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

DOUBLE PRECISION::T_DATA(DATA_NUM)               !the value of the independent variable at the sample points
DOUBLE PRECISION::P_DATA1(DIM_P,DATA_NUM)         !the value of the dependent variables at the sample points
DOUBLE PRECISION::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
DOUBLE PRECISION:: P_INTERP1(DIM_P,INTERP_NUM)      !the interpolated values of the dependent variables
DOUBLE PRECISION::ep,lambda,w0,gama,e0,pi,beta,x,xcommon,fcoef,betaL,betaR&
                  ,gamaR,gamaL,muL,muR,fermi,fermiaux,betaph,betaaux
DOUBLE PRECISION,DIMENSION(200)::w1
INTEGER::NEQ,k
DOUBLE PRECISION,DIMENSION(800)::YFRGcommon,YDFRG,SFRG,SDFRG
COMPLEX*16::ic,GLRFRG,GLAFRG,GLKFRG,RSelfinterp,ASelfinterp,KSelfinterp
COMPLEX*16::SLRinterp,SLAinterp,SLKinterp
!----------------------------------------------------------------------
!----------------------------------------------------------------------
ic=cmplx(0.0d+00,1.0d+00)
ep=(lambda**2)/w0

!print*,NEQ
!Do k=1,NEQ
Do k=1,DATA_NUM

T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA1(1,k)=YFRGcommon((DATA_NUM+1)-k)  
P_DATA1(2,k)=YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k)
P_DATA1(3,k)=YFRGcommon((2*DATA_NUM)+(DATA_NUM+1)-k)
P_DATA1(4,k)=YFRGcommon((3*DATA_NUM)+(DATA_NUM+1)-k)
!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------
T_INTERP(INTERP_NUM)=w2
call interp_linear(DIM_P,DATA_NUM,T_DATA,P_DATA1,INTERP_NUM,T_INTERP,P_INTERP1)
!----------------------------------------------------------------------
RSelfinterp=ic*P_INTERP1(1,1)+P_INTERP1(2,1)
ASelfinterp=CONJG(RSelfinterp)
KSelfinterp=ic*P_INTERP1(3,1)+P_INTERP1(4,1)

GLRFRG=1.0d+00/(w2+ic*(xcommon+gama)-fcoef*RSelfinterp-(1.0d+00-fcoef)*e0)
GLAFRG=1.0d+00/(w2-ic*(xcommon+gama)-fcoef*ASelfinterp-(1.0d+00-fcoef)*e0)
GLKFRG=GLRFRG*GLAFRG*(fcoef*KSelfinterp-2.0d+00*ic*gama*(1.0d+00-2.0d+00*fermi(w2))&
                       -2.0d+00*ic*xcommon*((1.0d+00-2.0d+00*fermiaux(w2))))

SLRinterp=ic*((GLRFRG)**2)
SLAinterp=CONJG(SLRinterp)!-ic*((GLAFRG)**2)
SLKinterp=ic*GLRFRG*GLKFRG-ic*GLAFRG*GLKFRG+GLRFRG*GLAFRG*(2.0d+00*ic*(1.0d+00-2.0d+00*fermiaux(w2)))


fRArecaushy=real(SLRinterp+SLAinterp)

end function fRArecaushy
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!*******Integral for going to the Imaginary Axis:***********************
!-----------------------------------------------------------------------
subroutine Fint(NEQ, y, Ymat, YDmat)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon

INTEGER::k,i,NEQ
DOUBLE PRECISION::x,xcommon,y,rhoy
DOUBLE PRECISION,DIMENSION(800)::YFRGcommon
DOUBLE PRECISION,DIMENSION(NEQ)::Ymat,YDmat
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,beta,fcoef,betaL,betaR&
                  ,gamaR,gamaL,muL,muR,betaph,betaaux
DOUBLE PRECISION::IntCaushy,ap,am,SHart
DOUBLE PRECISION,DIMENSION(200)::w1
COMPLEX*16::ic
!------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
INTEGER*4::DATA_NUM
INTEGER*4, parameter::DIM_P=2          !the spatial dimension.
!INTEGER*4, parameter::DATA_NUM=NEQ/2   !the number of data points
INTEGER*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

REAL*8::T_DATA(NEQ/2)               !the value of the independent variable at the sample points
REAL*8::P_DATA(DIM_P,NEQ/2)         !the value of the dependent variables at the sample points
REAL*8::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
REAL*8:: P_INTERP(DIM_P,INTERP_NUM)    !the interpolated values of the dependent variables
!----------------------------------------------------------------------
!----------------------------------------------------------------------
ep=(lambda**2)/w0
SHart=-ep
ic=cmplx(0.0d+00,1.0d+00)
DATA_NUM=200 
!----------------------------------------------------------------------
Do k=1,DATA_NUM
T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA(1,k)=YFRGcommon((DATA_NUM+1)-k)
P_DATA(2,k)=YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k)

! P_DATA(1,k)=(-1.0d+00/pi)*(-gama+(YFRGcommon((DATA_NUM+1)-k)))&
!              /((T_DATA(k)-YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k))**2+(gama-(YFRGcommon((DATA_NUM+1)-k)))**2) 


end do       
!----------------------------------------------------------------------
T_INTERP=y
call interp_linear(DIM_P,DATA_NUM,T_DATA,P_DATA,INTERP_NUM,T_INTERP,P_INTERP)
!---------------------------------------------------------------------- 
rhoy=(-1.0d+00/pi)*(-gama+(P_INTERP(1,1)))&
              /((y-(P_INTERP(2,1)))**2+(gama-(P_INTERP(1,1)))**2) 

!write(1,*)y/w0,rhoy*pi*gama,P_INTERP(1,1)/gama,P_INTERP(2,1)/gama

Do k=1,DATA_NUM

!YDmat(k)=-w1(k)*P_INTERP(1,1)/((w1(k))**2+y**2) 
YDmat(k)=-w1(k)*rhoy/((w1(k))**2+y**2) 

!YDmat((DATA_NUM)+k)=-y*P_INTERP(1,1)/((w1(k))**2+y**2) 
YDmat((DATA_NUM)+k)=-y*rhoy/((w1(k))**2+y**2) 

end do !k
!----------------------------------------------------------------------
end subroutine Fint
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!******* Transport Integrals: linear response ***********************
!-----------------------------------------------------------------------
subroutine FIn(NEQ, y, YIn, YDIn)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon

INTEGER::k,i,NEQ
DOUBLE PRECISION::x,xcommon,y,rhoy,drhoy
DOUBLE PRECISION,DIMENSION(800)::YFRGcommon
DOUBLE PRECISION,DIMENSION(NEQ)::YIn,YDIn
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,beta,fcoef,betaL,betaR&
                 ,gamaR,gamaL,muL,muR,fermi,fermiaux,betaph,betaaux
DOUBLE PRECISION::IntCaushy,ap,am,SHart,dimsigmaR,dresigmaR
DOUBLE PRECISION,DIMENSION(200)::w1
COMPLEX*16::ic
!------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
!INTEGER*4::DATA_NUM
INTEGER*4, parameter::DIM_P=2          !the spatial dimension.
INTEGER*4, parameter::DATA_NUM=200     !the number of data points
INTEGER*4, parameter::INTERP_NUM=2     !the number of points at which interpolation is to be done

REAL*8::T_DATA(DATA_NUM)               !the value of the independent variable at the sample points
REAL*8::P_DATA(DIM_P,DATA_NUM)         !the value of the dependent variables at the sample points
REAL*8::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
REAL*8:: P_INTERP(DIM_P,INTERP_NUM)    !the interpolated values of the dependent variables
!----------------------------------------------------------------------
!----------------------------------------------------------------------
ep=(lambda**2)/w0
SHart=-ep
ic=cmplx(0.0d+00,1.0d+00)
!DATA_NUM=200 
!----------------------------------------------------------------------
Do k=1,DATA_NUM
T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA(1,k)=YFRGcommon((DATA_NUM+1)-k)
P_DATA(2,k)=YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k)-(e0-ep)

! P_DATA(1,k)=(-1.0d+00/pi)*(-gama+(YFRGcommon((DATA_NUM+1)-k)))&
!              /((T_DATA(k)-YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k))**2+(gama-(YFRGcommon((DATA_NUM+1)-k)))**2) 


end do       
!----------------------------------------------------------------------
T_INTERP(1)=y
T_INTERP(2)=y+(w0/1000.0d+00)
call interp_linear(DIM_P,DATA_NUM,T_DATA,P_DATA,INTERP_NUM,T_INTERP,P_INTERP)
!---------------------------------------------------------------------- 

dimsigmaR=(P_INTERP(1,2)-P_INTERP(1,1))/(T_INTERP(2)-T_INTERP(1))
dresigmaR=(P_INTERP(2,2)-P_INTERP(2,1))/(T_INTERP(2)-T_INTERP(1))

rhoy=(-1.0d+00/pi)*(-gama+(P_INTERP(1,1)))&
              /((y-(e0-ep)-(P_INTERP(2,1)))**2+(gama-(P_INTERP(1,1)))**2) 

drhoy=(1.0d+00/pi)*((-dimsigmaR)*(((y-P_INTERP(2,1))**2)+((gama-P_INTERP(1,1))**2))&
                   -(gama-P_INTERP(1,1))*(2.0d+00*(y-P_INTERP(2,1))*(1.0d+00-dresigmaR)&
                                         +2.0d+00*(gama-P_INTERP(1,1))*(-dimsigmaR)))&
                 /((((y-P_INTERP(2,1))**2)+((gama-P_INTERP(1,1))**2))**2)

!write(1,*)y,rhoy*pi*gama

beta=(betaR+betaL)/2.0d+00

!----------------------------------------------------------------------
!I_0:
!----------------------------------------------------------------------
YDIn(1)=gama*pi*rhoy*beta/(4.0d+00*((cosh(beta*y/2.0d+00))**2))!
!----------------------------------------------------------------------
!I_1:
!----------------------------------------------------------------------
YDIn(2)=gama*pi*rhoy*y*beta/(4.0d+00*((cosh(beta*y/2.0d+00))**2))
!----------------------------------------------------------------------
!I_2:
!----------------------------------------------------------------------
YDIn(3)=gama*pi*rhoy*(y**2)*beta/(4.0d+00*((cosh(beta*y/2.0d+00))**2))
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end subroutine FIn
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***********************************************************************
DOUBLE PRECISION function funSK(w2)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon

DOUBLE PRECISION,INTENT(in)::w2

!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
INTEGER*4, parameter::DIM_P=4          !the spatial dimension.
integer*4, parameter::DATA_NUM=200    !the number of data points
INTEGER*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

DOUBLE PRECISION::T_DATA(DATA_NUM)               !the value of the independent variable at the sample points
DOUBLE PRECISION::P_DATA1(DIM_P,DATA_NUM)         !the value of the dependent variables at the sample points
DOUBLE PRECISION::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
DOUBLE PRECISION:: P_INTERP1(DIM_P,INTERP_NUM)      !the interpolated values of the dependent variables
DOUBLE PRECISION::ep,lambda,w0,gama,e0,pi,beta,x,xcommon,fcoef,betaL,betaR&
                  ,gamaR,gamaL,muL,muR,fermi,fermiaux,betaph,betaaux
DOUBLE PRECISION,DIMENSION(200)::w1

INTEGER::NEQ,k
DOUBLE PRECISION,DIMENSION(800)::YFRGcommon,YDFRG,SFRG,SDFRG
COMPLEX*16::ic,GLRFRG,GLAFRG,GLKFRG,RSelfinterp,ASelfinterp,KSelfinterp
COMPLEX*16::SLRinterp,SLAinterp,SLKinterp
!----------------------------------------------------------------------
!----------------------------------------------------------------------
ic=cmplx(0.0d+00,1.0d+00)
ep=(lambda**2)/w0

!print*,NEQ
!Do k=1,NEQ
Do k=1,DATA_NUM

T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA1(1,k)=YFRGcommon((DATA_NUM+1)-k)  
P_DATA1(2,k)=YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k)
P_DATA1(3,k)=YFRGcommon((2*DATA_NUM)+(DATA_NUM+1)-k)
P_DATA1(4,k)=YFRGcommon((3*DATA_NUM)+(DATA_NUM+1)-k)
!----------------------------------------------------------------------
end do
!----------------------------------------------------------------------
T_INTERP(INTERP_NUM)=w2
call interp_linear(DIM_P,DATA_NUM,T_DATA,P_DATA1,INTERP_NUM,T_INTERP,P_INTERP1)
!----------------------------------------------------------------------
RSelfinterp=ic*P_INTERP1(1,1)+P_INTERP1(2,1)
ASelfinterp=CONJG(RSelfinterp)
KSelfinterp=ic*P_INTERP1(3,1)+P_INTERP1(4,1)

GLRFRG=1.0d+00/(w2+ic*(xcommon+gama)-fcoef*RSelfinterp-(1.0d+00-fcoef)*e0)
GLAFRG=1.0d+00/(w2-ic*(xcommon+gama)-fcoef*ASelfinterp-(1.0d+00-fcoef)*e0)
GLKFRG=GLRFRG*GLAFRG*(fcoef*KSelfinterp-2.0d+00*ic*gama*(1.0d+00-2.0d+00*fermi(w2))&
                       -2.0d+00*ic*xcommon*((1.0d+00-2.0d+00*fermiaux(w2))))


SLRinterp=ic*((GLRFRG)**2)
SLAinterp=CONJG(SLRinterp)!-ic*((GLAFRG)**2)
SLKinterp=ic*GLRFRG*GLKFRG-ic*GLAFRG*GLKFRG+GLRFRG*GLAFRG*(2.0d+00*ic*(1.0d+00-2.0d+00*fermiaux(w2)))


funSK=aimag(SLKinterp)!*(w2-w0)

end function funSK
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
integer*8 function sgnn(w)
implicit none
double precision,intent(in)::w
if(w>0)then
sgnn=1
elseif(w==0)then
sgnn=1
else
sgnn=-1
endif
end function sgnn
!----------------------------------------------------------------------
!----------------------------------------------------------------------
double precision function fermi(w)
implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR

double precision,intent(in)::w
double precision::fermiL,fermiR,w0,lambda,gama,pi,e0,betaL,betaR&
                  ,muL,muR,gamaL,gamaR,betaph,betaaux
fermiL=(1.0d+00-tanh(betaL*(w-muL)/2.0d+00))/2.0d+00!1.0d+00/(1.0d+00+exp(betaL*(w-muL)))
fermiR=(1.0d+00-tanh(betaR*(w-muR)/2.0d+00))/2.0d+00!1.0d+00/(1.0d+00+exp(betaR*(w-muR)))
fermi=(gamaL/gama)*fermiL+(gamaR/gama)*fermiR
end function fermi
!----------------------------------------------------------------------
!----------------------------------------------------------------------
double precision function fermiaux(w)
implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR

double precision,intent(in)::w
double precision::fermiLaux,fermiRaux,w0,lambda,gama,pi,e0,betaL,betaR&
                  ,muL,muR,gamaL,gamaR,betaph,betaaux
fermiLaux=(1.0d+00-tanh(betaL*(w-muL)/2.0d+00))/2.0d+00!1.0d+00/(1.0d+00+exp(betaL*(w-muL)))
fermiRaux=(1.0d+00-tanh(betaR*(w-muR)/2.0d+00))/2.0d+00!1.0d+00/(1.0d+00+exp(betaR*(w-muR)))
fermiaux=(gamaL/gama)*fermiLaux+(gamaR/gama)*fermiRaux!(1.0d+00-tanh(betaL*(w)/2.0d+00))/2.0d+00!
end function fermiaux
!----------------------------------------------------------------------
!----------------------------------------------------------------------
double precision function bose(w)
implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR

double precision,intent(in)::w
double precision::fermiL,fermiR,w0,lambda,gama,pi,e0,betaL,betaR&
                  ,muL,muR,gamaL,gamaR,betaph,betaaux
bose=((1.0d+00/(tanh(betaph*w/2.0d+00)))-1.0d+00)/2.0d+00
end function bose
!-------------------------------------------------------------
!-----------------------------------------------------------------------
!******* currents integrals***********************
!-----------------------------------------------------------------------
subroutine Fcurrent(NEQ, y, Ycur, YDcur)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon

INTEGER::k,i,NEQ
DOUBLE PRECISION::x,xcommon,y,rhoy,drhoy,fermiL,fermiR,disy
DOUBLE PRECISION,DIMENSION(800)::YFRGcommon
DOUBLE PRECISION,DIMENSION(NEQ)::Ycur,YDcur
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,beta,fcoef,betaL,betaR&
                 ,gamaR,gamaL,muL,muR,fermi,fermiaux,betaph,betaaux
DOUBLE PRECISION::IntCaushy,ap,am,SHart
DOUBLE PRECISION,DIMENSION(200)::w1
COMPLEX*16::ic
!------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
!INTEGER*4::DATA_NUM
INTEGER*4, parameter::DIM_P=3          !the spatial dimension.
INTEGER*4, parameter::DATA_NUM=200     !the number of data points
INTEGER*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

REAL*8::T_DATA(DATA_NUM)               !the value of the independent variable at the sample points
REAL*8::P_DATA(DIM_P,DATA_NUM)         !the value of the dependent variables at the sample points
REAL*8::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
REAL*8:: P_INTERP(DIM_P,INTERP_NUM)    !the interpolated values of the dependent variables
!----------------------------------------------------------------------
!----------------------------------------------------------------------
ep=(lambda**2)/w0
SHart=-ep
ic=cmplx(0.0d+00,1.0d+00)
!DATA_NUM=200 
!----------------------------------------------------------------------
Do k=1,DATA_NUM
T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA(1,k)=YFRGcommon((DATA_NUM+1)-k)
P_DATA(2,k)=YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k)-(e0-ep)
P_DATA(3,k)=YFRGcommon((2*DATA_NUM)+(DATA_NUM+1)-k)

! P_DATA(1,k)=(-1.0d+00/pi)*(-gama+(YFRGcommon((DATA_NUM+1)-k)))&
!              /((T_DATA(k)-YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k))**2+(gama-(YFRGcommon((DATA_NUM+1)-k)))**2) 


end do       
!----------------------------------------------------------------------
T_INTERP(1)=y
call interp_linear(DIM_P,DATA_NUM,T_DATA,P_DATA,INTERP_NUM,T_INTERP,P_INTERP)
!---------------------------------------------------------------------- 
fermiL=(1.0d+00-tanh(betaL*(y-muL)/2.0d+00))/2.0d+00!1.0d+00/(1.0d+00+exp(betaL*(w-muL)))
fermiR=(1.0d+00-tanh(betaR*(y-muR)/2.0d+00))/2.0d+00!1.0d+00/(1.0d+00+exp(betaR*(w-muR)))

fermi=(gamaL/gama)*fermiL+(gamaR/gama)*fermiR
!---------------------------------------------------------------------- 

disy=(0.5d+00)-((P_INTERP(3,1)-(2.0d+00*gama*(1.0d+00-2.0d+00*fermi)))/(4.0d+00*(-gama+P_INTERP(1,1))))

rhoy=(-1.0d+00/pi)*(-gama+(P_INTERP(1,1)))&
              /((y-(e0-ep)-(P_INTERP(2,1)))**2+(gama-(P_INTERP(1,1)))**2) 

!write(1,*)y,disy              
!----------------------------------------------------------------------
!charge-current:
!----------------------------------------------------------------------
YDcur(1)=(4.0d+00*pi)*((gamaR*gamaL)/gama)*rhoy*(fermiL-fermiR)
!----------------------------------------------------------------------
!Heat-current-Right:
!----------------------------------------------------------------------
YDcur(2)=(-4.0d+00*pi)*(gamaR)*rhoy*(fermiR-disy)*(y-0.0d+00*muR)
!----------------------------------------------------------------------
!Hear-current-Left:
!----------------------------------------------------------------------
YDcur(3)=(-4.0d+00*pi)*(gamaL)*rhoy*(fermiL-disy)*(y-0.0d+00*muL)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
YDcur(4)=(4.0d+00*pi)*(gama)*rhoy*(fermi-disy)
!----------------------------------------------------------------------
end subroutine Fcurrent
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!******* Transport Integrals:***********************
!-----------------------------------------------------------------------
subroutine FInL(NEQ, y, YInL, YDInL)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon

INTEGER::k,i,NEQ
DOUBLE PRECISION::x,xcommon,y,rhoy,drhoy
DOUBLE PRECISION,DIMENSION(800)::YFRGcommon
DOUBLE PRECISION,DIMENSION(NEQ)::YInL,YDInL
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,beta,fcoef,betaph,betaaux&
                  ,betaL,betaR,gamaR,gamaL,muL,muR
DOUBLE PRECISION::IntCaushy,ap,am,SHart,dimsigmaR,dresigmaR
DOUBLE PRECISION,DIMENSION(200)::w1
COMPLEX*16::ic
!------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
!INTEGER*4::DATA_NUM
INTEGER*4, parameter::DIM_P=2          !the spatial dimension.
INTEGER*4, parameter::DATA_NUM=200     !the number of data points
INTEGER*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

REAL*8::T_DATA(DATA_NUM)               !the value of the independent variable at the sample points
REAL*8::P_DATA(DIM_P,DATA_NUM)         !the value of the dependent variables at the sample points
REAL*8::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
REAL*8:: P_INTERP(DIM_P,INTERP_NUM)    !the interpolated values of the dependent variables
!----------------------------------------------------------------------
!----------------------------------------------------------------------
ep=(lambda**2)/w0
SHart=-ep
ic=cmplx(0.0d+00,1.0d+00)
!DATA_NUM=200 
!----------------------------------------------------------------------
Do k=1,DATA_NUM
T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA(1,k)=YFRGcommon((DATA_NUM+1)-k)
P_DATA(2,k)=YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k)-(e0-ep)

! P_DATA(1,k)=(-1.0d+00/pi)*(-gama+(YFRGcommon((DATA_NUM+1)-k)))&
!              /((T_DATA(k)-YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k))**2+(gama-(YFRGcommon((DATA_NUM+1)-k)))**2) 


end do       
!----------------------------------------------------------------------
T_INTERP(1)=y
call interp_linear(DIM_P,DATA_NUM,T_DATA,P_DATA,INTERP_NUM,T_INTERP,P_INTERP)
!---------------------------------------------------------------------- 

rhoy=(-1.0d+00/pi)*(-gama+(P_INTERP(1,1)))&
              /((y-(e0-ep)-(P_INTERP(2,1)))**2+(gama-(P_INTERP(1,1)))**2) 

!----------------------------------------------------------------------
!I_0:
!----------------------------------------------------------------------
YDInL(1)=(4.0d+00*pi*gamaL*gamaR/gama)*rhoy*(betaL)/(4.0d+00*(cosh(betaL*(y-muL)/2.0d+00))**2)
!gama*pi*rhoy*beta/(4.0d+00*((cosh(beta*y/2.0d+00))**2))!
!----------------------------------------------------------------------
!I_1:
!----------------------------------------------------------------------
YDInL(2)=(4.0d+00*pi*gamaL*gamaR/gama)*rhoy*(y)*(betaL)/(4.0d+00*(cosh(betaL*(y-muL)/2.0d+00))**2)
!gama*pi*rhoy*y*beta/(4.0d+00*((cosh(beta*y/2.0d+00))**2))
!----------------------------------------------------------------------
!I_2:
!----------------------------------------------------------------------
YDInL(3)=(4.0d+00*pi*gamaL*gamaR/gama)*rhoy*((y)**2)*(betaL)/(4.0d+00*(cosh(betaL*(y-muL)/2.0d+00))**2)
!gama*pi*rhoy*(y**2)*beta/(4.0d+00*((cosh(beta*y/2.0d+00))**2))
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end subroutine FInL
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine FInR(NEQ, y, YInR, YDInR)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon

INTEGER::k,i,NEQ
DOUBLE PRECISION::x,xcommon,y,rhoy,drhoy
DOUBLE PRECISION,DIMENSION(800)::YFRGcommon
DOUBLE PRECISION,DIMENSION(NEQ)::YInR,YDInR
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,beta,fcoef,betaph,betaaux&
                  ,betaL,betaR,gamaR,gamaL,muL,muR
DOUBLE PRECISION::IntCaushy,ap,am,SHart,dimsigmaR,dresigmaR
DOUBLE PRECISION,DIMENSION(200)::w1
COMPLEX*16::ic
!------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
!INTEGER*4::DATA_NUM
INTEGER*4, parameter::DIM_P=2          !the spatial dimension.
INTEGER*4, parameter::DATA_NUM=200     !the number of data points
INTEGER*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

REAL*8::T_DATA(DATA_NUM)               !the value of the independent variable at the sample points
REAL*8::P_DATA(DIM_P,DATA_NUM)         !the value of the dependent variables at the sample points
REAL*8::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
REAL*8:: P_INTERP(DIM_P,INTERP_NUM)    !the interpolated values of the dependent variables
!----------------------------------------------------------------------
!----------------------------------------------------------------------
ep=(lambda**2)/w0
SHart=-ep
ic=cmplx(0.0d+00,1.0d+00)
!DATA_NUM=200 
!----------------------------------------------------------------------
Do k=1,DATA_NUM
T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA(1,k)=YFRGcommon((DATA_NUM+1)-k)
P_DATA(2,k)=YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k)-(e0-ep)

! P_DATA(1,k)=(-1.0d+00/pi)*(-gama+(YFRGcommon((DATA_NUM+1)-k)))&
!              /((T_DATA(k)-YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k))**2+(gama-(YFRGcommon((DATA_NUM+1)-k)))**2) 


end do       
!----------------------------------------------------------------------
T_INTERP(1)=y
call interp_linear(DIM_P,DATA_NUM,T_DATA,P_DATA,INTERP_NUM,T_INTERP,P_INTERP)
!---------------------------------------------------------------------- 

rhoy=(-1.0d+00/pi)*(-gama+(P_INTERP(1,1)))&
              /((y-(e0-ep)-(P_INTERP(2,1)))**2+(gama-(P_INTERP(1,1)))**2) 

!----------------------------------------------------------------------
!I_0:
!----------------------------------------------------------------------
YDInR(1)=(4.0d+00*pi*gamaL*gamaR/gama)*rhoy*(betaR)/(4.0d+00*(cosh(betaR*(y-muR)/2.0d+00))**2)
!gama*pi*rhoy*beta/(4.0d+00*((cosh(beta*y/2.0d+00))**2))!
!----------------------------------------------------------------------
!I_1:
!----------------------------------------------------------------------
YDInR(2)=(4.0d+00*pi*gamaL*gamaR/gama)*rhoy*(y)*(betaR)/(4.0d+00*(cosh(betaR*(y-muR)/2.0d+00))**2)
!gama*pi*rhoy*y*beta/(4.0d+00*((cosh(beta*y/2.0d+00))**2))
!----------------------------------------------------------------------
!I_2:
!----------------------------------------------------------------------
YDInR(3)=(4.0d+00*pi*gamaL*gamaR/gama)*rhoy*((y)**2)*(betaR)/(4.0d+00*(cosh(betaR*(y-muR)/2.0d+00))**2)
!gama*pi*rhoy*(y**2)*beta/(4.0d+00*((cosh(beta*y/2.0d+00))**2))
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end subroutine FInR
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine FIncor(NEQ, y, YIncor, YDIncor)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon,dnespec

INTEGER::k,i,NEQ
DOUBLE PRECISION::x,xcommon,y,drhoy
DOUBLE PRECISION,DIMENSION(800)::YFRGcommon
DOUBLE PRECISION,DIMENSION(NEQ)::YIncor,YDIncor
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,beta,fcoef,betaph,betaaux&
                  ,betaL,betaR,gamaR,gamaL,muL,muR
DOUBLE PRECISION::IntCaushy,ap,am,SHart,dimsigmaR,dresigmaR
DOUBLE PRECISION,DIMENSION(200)::w1
COMPLEX*16::ic
!------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
!INTEGER*4::DATA_NUM
INTEGER*4, parameter::DIM_P=1          !the spatial dimension.
INTEGER*4, parameter::DATA_NUM=200     !the number of data points
INTEGER*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

REAL*8::T_DATA(DATA_NUM)               !the value of the independent variable at the sample points
REAL*8::P_DATA(DIM_P,DATA_NUM)         !the value of the dependent variables at the sample points
REAL*8::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
REAL*8:: P_INTERP(DIM_P,INTERP_NUM)    !the interpolated values of the dependent variables
!----------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:DATA_NUM)::dnespec
DOUBLE PRECISION::fermiL,fermiR
!----------------------------------------------------------------------
ep=(lambda**2)/w0
SHart=-ep
ic=cmplx(0.0d+00,1.0d+00)
!DATA_NUM=200 
!----------------------------------------------------------------------
Do k=1,DATA_NUM
T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA(1,k)=dnespec((DATA_NUM+1)-k)
end do       
!----------------------------------------------------------------------
T_INTERP(1)=y
call interp_linear(DIM_P,DATA_NUM,T_DATA,P_DATA,INTERP_NUM,T_INTERP,P_INTERP)
!---------------------------------------------------------------------- 
drhoy=P_INTERP(1,1)
!---------------------------------------------------------------------- 
fermiL=(1.0d+00-tanh(betaL*(y-muL)/2.0d+00))/2.0d+00!1.0d+00/(1.0d+00+exp(betaL*(w-muL)))
fermiR=(1.0d+00-tanh(betaR*(y-muR)/2.0d+00))/2.0d+00!1.0d+00/(1.0d+00+exp(betaR*(w-muR)))
!---------------------------------------------------------------------- 
!----------------------------------------------------------------------
!I_0:
!----------------------------------------------------------------------
YDIncor(1)=(4.0d+00*pi*gamaL*gamaR/gama)*drhoy*(fermiL-fermiR)
!----------------------------------------------------------------------
!I_1:
!----------------------------------------------------------------------
YDIncor(2)=(4.0d+00*pi*gamaL*gamaR/gama)*drhoy*(y)*(fermiL)
!----------------------------------------------------------------------
YDIncor(3)=(4.0d+00*pi*gamaL*gamaR/gama)*drhoy*(y)*(fermiR)
!----------------------------------------------------------------------
end subroutine FIncor
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine FInnedis(NEQ, y, YInnedis, YDInnedis)

implicit none

common w0,lambda,gama,pi,e0,betaph,betaaux,betaL,betaR,gamaR,gamaL,muL,muR,xcommon,fcoef,w1,YFRGcommon&
       ,dnespec,dneSigmaRre,dneSigmaRim,dneSigmaK

INTEGER::k,i,NEQ
DOUBLE PRECISION::x,xcommon,y,rhoy,drhoy,dneSigmaRreint,dneSigmaRimint,dneSigmaKint,dnedisy,dfR,dfL,disy,drhoy2
DOUBLE PRECISION,DIMENSION(800)::YFRGcommon
DOUBLE PRECISION,DIMENSION(NEQ)::YInnedis, YDInnedis
DOUBLE PRECISION::lambda,gama,zero,pi,w0,ep,e0,beta,fcoef,betaph,betaaux&
                  ,betaL,betaR,gamaR,gamaL,muL,muR
DOUBLE PRECISION::IntCaushy,ap,am,SHart,dimsigmaR,dresigmaR
DOUBLE PRECISION,DIMENSION(200)::w1
COMPLEX*16::ic
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameters for INTERP library:
!----------------------------------------------------------------------
!INTEGER*4::DATA_NUM
INTEGER*4, parameter::DIM_P=7          !the spatial dimension.
INTEGER*4, parameter::DATA_NUM=200     !the number of data points
INTEGER*4, parameter::INTERP_NUM=1     !the number of points at which interpolation is to be done

REAL*8::T_DATA(DATA_NUM)               !the value of the independent variable at the sample points
REAL*8::P_DATA(DIM_P,DATA_NUM)         !the value of the dependent variables at the sample points
REAL*8::T_INTERP(INTERP_NUM)           !the value of the independent variable at the interpolation points
REAL*8:: P_INTERP(DIM_P,INTERP_NUM)    !the interpolated values of the dependent variables
!----------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:DATA_NUM)::dnespec
DOUBLE PRECISION::fermiL,fermiR,fermi
DOUBLE PRECISION,DIMENSION(1:DATA_NUM)::dneSigmaRre,dneSigmaRim,dneSigmaK
!----------------------------------------------------------------------
ep=(lambda**2)/w0
SHart=-ep
ic=cmplx(0.0d+00,1.0d+00)
!DATA_NUM=200 
!----------------------------------------------------------------------
!---------------------------------------------------------------------- 
!----------------------------------------------------------------------
Do k=1,DATA_NUM
T_DATA(k)=w1((DATA_NUM+1)-k)
P_DATA(1,k)=YFRGcommon((DATA_NUM+1)-k)
P_DATA(2,k)=YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k)-(e0-ep)
P_DATA(3,k)=YFRGcommon((2*DATA_NUM)+(DATA_NUM+1)-k)

! derivatives:

P_DATA(4,k)=dnespec((DATA_NUM+1)-k)
P_DATA(5,k)=dneSigmaRre((DATA_NUM+1)-k)
P_DATA(6,k)=dneSigmaRim((DATA_NUM+1)-k)
P_DATA(7,k)=dneSigmaK((DATA_NUM+1)-k)


! P_DATA(1,k)=(-1.0d+00/pi)*(-gama+(YFRGcommon((DATA_NUM+1)-k)))&
!              /((T_DATA(k)-YFRGcommon((1*DATA_NUM)+(DATA_NUM+1)-k))**2+(gama-(YFRGcommon((DATA_NUM+1)-k)))**2) 

end do       
!----------------------------------------------------------------------    
!----------------------------------------------------------------------
T_INTERP(1)=y
call interp_linear(DIM_P,DATA_NUM,T_DATA,P_DATA,INTERP_NUM,T_INTERP,P_INTERP)
!---------------------------------------------------------------------- 
fermiL=(1.0d+00-tanh(betaL*(y-muL)/2.0d+00))/2.0d+00!1.0d+00/(1.0d+00+exp(betaL*(w-muL)))
fermiR=(1.0d+00-tanh(betaR*(y-muR)/2.0d+00))/2.0d+00!1.0d+00/(1.0d+00+exp(betaR*(w-muR)))

fermi=(gamaL/gama)*fermiL+(gamaR/gama)*fermiR
!---------------------------------------------------------------------- 
if(betaL>betaR)then
dfR=((betaR**2)*(y-muR))/(4.0d+00*(cosh(betaR*(y-muR)/2.0d+00))**2)
dfL=0.0d+00
else
dfR=(-betaR)/(8.0d+00*(cosh(betaR*(y-muR)/2.0d+00))**2)
dfL=(betaL)/(8.0d+00*(cosh(betaL*(y-muL)/2.0d+00))**2)
endif

disy=(0.5d+00)-((P_INTERP(3,1)-(2.0d+00*gama*(1.0d+00-2.0d+00*fermi)))/(4.0d+00*(-gama+P_INTERP(1,1))))

rhoy=(-1.0d+00/pi)*(-gama+(P_INTERP(1,1)))&
              /((y-(e0-ep)-(P_INTERP(2,1)))**2+(gama-(P_INTERP(1,1)))**2) 
   
   
drhoy=P_INTERP(4,1)   
dneSigmaRreint=P_INTERP(5,1) 
dneSigmaRimint=P_INTERP(6,1) 
dneSigmaKint=P_INTERP(7,1) 
!----------------------------------------------------------------------
drhoy2=(-1.0d+00/pi)*(((dneSigmaRimint)/((y-(e0-ep)-(P_INTERP(2,1)))**2+(gama-(P_INTERP(1,1)))**2))&
                     +((gama-(P_INTERP(1,1)))*((2.0d+00*(y-(e0-ep)-(P_INTERP(2,1)))*(-dneSigmaRreint))&
                                              +(2.0d+00*(gama-(P_INTERP(1,1)))*(-dneSigmaRimint)))&
                                            /(((y-(e0-ep)-(P_INTERP(2,1)))**2+(gama-(P_INTERP(1,1)))**2)**2)))

!----------------------------------------------------------------------
dnedisy=(-(dneSigmaKint+(4.0d+00*(gamaR*dfR+gamaL*dfL)))/(4.0d+00*(-gama+P_INTERP(1,1))))&
       +(((P_INTERP(3,1)-(2.0d+00*gama*(1.0d+00-2.0d+00*fermi)))*(dneSigmaRimint))/(4.0d+00*((-gama+P_INTERP(1,1))**2)))
!----------------------------------------------------------------------
!I_0:
!----------------------------------------------------------------------
YDInnedis(1)=(4.0d+00*pi*gamaL*gamaR/gama)*drhoy*(disy)*y
!----------------------------------------------------------------------
!I_1:
!----------------------------------------------------------------------
YDInnedis(2)=(4.0d+00*pi*gamaL*gamaR/gama)*rhoy*(y)*(dnedisy)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end subroutine FInnedis
!----------------------------------------------------------------------

