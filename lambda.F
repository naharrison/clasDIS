       real function lambda(ip,beam,Mcut)
       implicit none
#include "ntupgkoko.inc"
#include "claslepto.inc"
       integer ip
       real Mcut
c
c      calculating azimuthal angle from lab variables
c      of hadron  
c
c      input variables: pie,pit,pif-momentum,theta,phi of hadron
c                       eleq2,elee,elef-Q^2,momentum,phi of eleitron                 
c
c     CALL CROSS(A,B,C) C=[AxB]
c     VDOT(A,B,N), CALL VMUL(A,B,X,N) X_i=A_i.B_i  N=3  reals
c     VDOTN(A,B,N)=ab/|a||b|
c     VMOD (A,N)  =|a|
       real pi,tetgam,anu,pien,beam
       real eleq,pitg1,cospff1
       real xm,amp,tetgam1,elety,ebeam,vdotm
       real pi4(4),qiu4(4),el04(4),elf4(4),tnorm(4)
       real pro4(4),tnorm2(4)
       real vmass,pt,phi,costh

       integer m
       real imass
       real theta,mom,pcut

       pi=acos(-1.0)
       amp=0.938
       ebeam=beam
c
c     define all 4momenta
c     
       if(ip.eq.0) then         !data
c     if(npim.le.0.or.npro.le.0) print *,'wrong set'
c
c pi-
          pien=sqrt(pime*pime-0.139*0.139)
          pi4(4)=pime
          pi4(1)=pien*cos(pimf)*sin(pimt) 
          pi4(2)=pien*sin(pimf)*sin(pimt)
          pi4(3)=pien*cos(pimt)
c
c pro
          pien=sqrt(proe*proe-0.88035)
          pro4(4)=proe
          pro4(1)=pien*cos(prof)*sin(prot) 
          pro4(2)=pien*sin(prof)*sin(prot)
          pro4(3)=pien*cos(prot)

          lambda=0.
          return
       else if(ip.eq.1) then    ! mc Lambda
c          if(gnpim.le.0.or.gnpro.le.0) print *,'wrong pi-',gnpim,gnpro
c pro
          pien=sqrt(gproe*gproe-0.88035)
          pro4(4)=gproe
          pro4(1)=pien*cos(gprof)*sin(gprot) 
          pro4(2)=pien*sin(gprof)*sin(gprot)
          pro4(3)=pien*cos(gprot)
          theta=gprot*180./pi
c          if ((pien.lt.0).or.(pien.gt.4.).or.(theta.lt.14)) then
c             lambda=0.
c             return 
c          endif

c pi-
          do m=1,N
             if (k(m,2).eq.-211) then
                pi4(1)=P(m,1)
                pi4(2)=P(m,2)
                pi4(3)=P(m,3)
                pi4(4)=P(m,4)
                mom= sqrt(pi4(1)*pi4(1)+pi4(2)*pi4(2)+pi4(3)*pi4(3))
                theta=acos(pi4(3)/mom)*180./pi
c	        pcut=0.63-0.009*theta;
c                if ((mom.le.1.5).and.(theta.ge.20).and.(mom.ge.pcut)) then

                   call vsumm(pro4,pi4,tnorm,4)
                   if (vmass(tnorm).le.Mcut) then
                      lambda=1.
                      return
                   endif

c                endif

             endif
          enddo

       else if(ip.eq.2) then    ! mc Lambda
c          if(gnpi0.le.0.or.gnneu.le.0) print *,'wrong pi0',gnpi0,gnpro
c     
c
c neutron
          pien=sqrt(gneue*gneue-0.88285)
          pro4(4)=gneue
          pro4(1)=pien*cos(gneuf)*sin(gneut) 
          pro4(2)=pien*sin(gneuf)*sin(gneut)
          pro4(3)=pien*cos(gneut)

c pi0
          do m=1,N
             if (k(m,2).eq.111) then
                pi4(1)=P(m,1)
                pi4(2)=P(m,2)
                pi4(3)=P(m,3)
                pi4(4)=P(m,4)

                call vsumm(pro4,pi4,tnorm,4)
                if (vmass(tnorm).le.Mcut) then
                   lambda=1.
                   return
                endif

             endif
          enddo

       else
          lambda=0
          print *,'Unknown Lambda',ip
          stop
       endif
c
c     
c       lambda=vmass(tnorm)
       lambda=0.

       return
       end

