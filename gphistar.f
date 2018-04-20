       real function gphistar(ip,beam)
       implicit none
       include 'nt22.inc'
       integer ip
c
c      calculating azimuthal angle from lab variables
c      of hadron  
c
c      input variables: pie,pit,pif-momentum,theta,phi of hadron
c                       posq2,pose,posf-Q^2,momentum,phi of positron                 
c
c     CALL CROSS(A,B,C) C=[AxB]
c     VDOT(A,B,N), CALL VMUL(A,B,X,N) X_i=A_i.B_i  N=3  reals
c     VDOTN(A,B,N)=ab/|a||b|
c     VMOD (A,N)  =|a|
       real pi,tetgam,anu,pien,beam
       real eleq,pitg1,cospff1
       real xm,amp,tetgam1,elety,ebeam
       real pi4(4),qiu4(4),el04(4),elf4(4),tnorm(4)
       real pro4(4),tnorm2(4)
       real vmass,vangle,vdotm,phigstar
       pi=acos(-1.0)

        amp=0.938
        ebeam=beam

c
c     define all 4momenta
c
       if(ip.eq.1) then
c
c      pi+
c
        if(gnpip.eq.0) print *,'Bad pip'
         pien=gpipe
         pi4(4)=sqrt(pien*pien+0.135*0.135)
         pi4(1)=pien*cos(gpipf)*sin(gpipt) 
         pi4(2)=pien*sin(gpipf)*sin(gpipt)
         pi4(3)=pien*cos(gpipt)
c

       elseif(ip.eq.-1) then
c      pi-
c
        if(gnpim.eq.0) print *,'Bad pim'
        pien=gpime
        pi4(4)=sqrt(pien*pien+0.135*0.135)
        pi4(1)=pien*cos(gpimf)*sin(gpimt) 
        pi4(2)=pien*sin(gpimf)*sin(gpimt)
        pi4(3)=pien*cos(gpimt)
c
       elseif(ip.eq.0) then
c      pi0
c
        if(gnpi0.eq.0) print *,'Bad pi0'
        pien=gpi0e
        pi4(4)=sqrt(pien*pien+0.135*0.135)
        pi4(1)=pien*cos(gpi0f)*sin(gpi0t) 
        pi4(2)=pien*sin(gpi0f)*sin(gpi0t)
        pi4(3)=pien*cos(gpi0t)
c
       elseif(ip.eq.2) then
c      real gamma
c
        if(gngam.eq.0) print *,'Bad gamma'
        pien=ggame
        pi4(4)=pien
        pi4(1)=pien*cos(ggamf)*sin(ggamt) 
        pi4(2)=pien*sin(ggamf)*sin(ggamt)
        pi4(3)=pien*cos(ggamt)
c

       elseif(ip.eq.4) then
c      proton
c
        if(gnpro.eq.0) print *,'Bad pro'
        pien=gproe
        pi4(4)=sqrt(pien*pien+amp*amp)
        pi4(1)=pien*cos(gprof)*sin(gprot) 
        pi4(2)=pien*sin(gprof)*sin(gprot)
        pi4(3)=pien*cos(gprot)
c
c
       elseif(ip.eq.6) then
c      rho 0 
c
        if(gnpip*gnpim.eq.0) print *,'mismas:Bad rho',gnpip,gnpim
        pien=sqrt(gpipe*gpipe+0.0193)+sqrt(gpime*gpime+0.0193)
        pi4(4)=pien
        pi4(1)=gpipe*cos(gpipf)*sin(gpipt)+gpime*cos(gpimf)*sin(gpimt) 
        pi4(2)=gpipe*sin(gpipf)*sin(gpipt)+gpime*sin(gpimf)*sin(gpimt)
        pi4(3)=gpipe*cos(gpipt)+gpime*cos(gpimt)
c
       elseif(ip.eq.7) then
c      rho + 
c
        if(gnpip*gnpi0.eq.0) print *,'Bad rho+'
        pien=sqrt(gpipe*gpipe+0.0193)+sqrt(gpi0e*gpi0e+0.0182)
        pi4(4)=pien
        pi4(1)=gpipe*cos(gpipf)*sin(gpipt)+gpi0e*cos(gpi0f)*sin(gpi0t) 
        pi4(2)=gpipe*sin(gpipf)*sin(gpipt)+gpi0e*sin(gpi0f)*sin(gpi0t)
        pi4(3)=gpipe*cos(gpipt)+gpi0e*cos(gpi0t)
c
       elseif(ip.eq.-7) then
c      rho - 
c
        if(gnpim*gnpi0.eq.0) print *,'Bad rho-'
        pien=sqrt(gpime*gpime+0.0193)+sqrt(gpi0e*gpi0e+0.0182)
        pi4(4)=pien
        pi4(1)=gpime*cos(gpimf)*sin(gpimt)+gpi0e*cos(gpi0f)*sin(gpi0t) 
        pi4(2)=gpime*sin(gpimf)*sin(gpimt)+gpi0e*sin(gpi0f)*sin(gpi0t)
        pi4(3)=gpime*cos(gpimt)+gpi0e*cos(gpi0t)
       endif

c       Now el and gamma
c
c
c       el0
c
        el04(1)=0
        el04(2)=0
        el04(3)=ebeam
        el04(4)=ebeam
        elf4(1)=gelee*sin(gelet)*cos(gelef)
        elf4(2)=gelee*sin(gelet)*sin(gelef)
        elf4(3)=gelee*cos(gelet)
        elf4(4)=gelee
c
c     gamma*
c
c proton
c
       pro4(4)=amp
       pro4(3)=0
       pro4(2)=0
       pro4(1)=0
c
c
        anu=ebeam-abs(gelee)
        eleq=sqrt(gq2+anu*anu)
        elety=anu/ebeam
        xm=(amp*gxb)**2
        tetgam1=sqrt(4*xm*(1.0-elety-xm*elety*elety/gq2)/(gq2+4*xm))
        tetgam=asin(tetgam1)
         pitg1=tetgam

         call vdifm(el04,elf4,qiu4,4)
        
c$$$        qiu4(4)=anu
c$$$        qiu4(1)=eleq*cos(pi+elef)*sin(pitg1) ! photon momentum
c$$$        qiu4(2)=eleq*sin(pi+elef)*sin(pitg1)
c$$$        qiu4(3)=eleq*cos(pitg1)
c$$$c
         call crossm(qiu4,el04,tnorm)
         if (vdotm(tnorm,pi4,3).ge.0) then    ! for pi+
          gphistar=vangle(qiu4,el04,qiu4,pi4)
         else
          gphistar=2*pi-vangle(qiu4,el04,qiu4,pi4)
         endif

         return
         end



c/HISTOGRAM/CREATE/PROFILE '103' ' ' '6' '0.2' '2' '-0.1' '0.1'
c        A set of vector utilities
c
c
c 1.    call crossm(a,b,c)
c       c=a x b for  3 vectors a,b
c
c 2.    real function vangle(a,b,c,d)
c      angle between 2 planes a x b  and c x d
c
c 3.    call vmult(a,r,n)
c      a=r*a n-dimension vectors a multiply with real number r
c
c 4.    call vsumm(a,b,c,n)
c      c=a-b n-dimension vectors a,b
c
c 5.   call vdifm(a,b,c,n)
c      c=a-b n-dimension vectors a,b
c
c 6.   real function vdotm(a,b,n)
c      return a.b for n a,b vectors of n -dimension

       subroutine crossm(a,b,c)
       real a(3),b(3),c(3)
       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       return
       end
c
c
       subroutine vmult(a,r,n)
       real a(n),r
       integer i,n
        do i=1,n
          a(i)=r*a(i)
        enddo
       return
       end
c
c
       subroutine vsumm(a,b,c,n)
       real a(n),b(n),c(n)
       integer i,n
        do i=1,n
          c(i)=a(i)+b(i)
        enddo
       return
       end
c
       subroutine vdifm(a,b,c,n)
       real a(n),b(n),c(n)
       integer i,n
        do i=1,n
          c(i)=a(i)-b(i)
        enddo
       return
       end
c
c
       real function vdotm(a,b,n)
       real a(n),b(n),s
       integer i,n
       s=0.0
       do i=1,3
         s=s+a(i)*b(i)
       enddo
       if(n.eq.4) s=s-a(n)*b(n)
       vdotm=s
       return
       end
c
c
       real function vmass(a)
       real vm
        vm= vdotm(a,a,4)
        if (vm.lt.0.0) then
          vmass=sqrt(-vm)
        else
          vmass=-1.0
        endif 
       return
       end
c   
       real function vangle(a,b,c,d)
       real a(3),b(3),c(3),d(3),xm,ym,vcos
       real x(3),y(3),pi
       pi=acos(-1.0)
       call crossm(a,b,x)
       call crossm(c,d,y)
       xm=vdotm(x,x,3)
       ym=vdotm(y,y,3)
       if(xm.gt.0.0 .and. ym.gt.0.0) then
         vcos=vdotm(x,y,3)/sqrt(xm)/sqrt(ym)
         if(abs(vcos).lt.1.0) then
            vangle=acos(vcos)
         else
            if(vcos.ge.1.0)  vangle=0
            if(vcos.le.-1.0)  vangle=pi
         endif 
       else
         vangle=0
       endif
       return
       end
