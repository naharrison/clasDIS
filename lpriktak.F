
********************************************************************

      subroutine lpriktak(s,pt,ptmaxpr,phi)
c
      implicit none
c
cak modified to include cahn cos(phi) correction and sivers effect
c...size (pt) and azimuthal angle (phi) of primordial kt according
c...to a gaussian distribution.

      real s,pt,ptmaxpr,phi
c
      integer  nlupdm,nplbuf
      parameter  (nlupdm=4000,nplbuf=5)
      common/lujets/n,k(nlupdm,5),p(nlupdm,nplbuf),v(nlupdm,5)
      integer  n,k
      real  p,v
c
      common /leptou/ cut(14),lst(40),parl(30),
     &                x,y,w2,q2,u
      real cut,parl,x,y,w2,q2,u
      integer lst
c
      real ptq,phiq,phis,sivq,st
      common/quark/ptq,phiq,phis,sivq,st
      integer iazuser,iaz,kfq,i
      real pi,twopi,phimu,rlu
      common/twopi/pi,twopi
      common/azopt/iazuser,phimu
      real pbeam,pvg
      common/labmom/pbeam(3),pvg(3)
      real cahn,cahnmax,cahn1,cahn2
      real dstr,q2pt,q2pt2,y1,y11,y12,y2
      
      iaz=0
      if(iazuser .ne. 0)iaz=iazuser
 111  pt=s*sqrt(-alog(rlu(0)))
      if(pt.gt.ptmaxpr)goto 111
      ptq=pt

      if(iaz .eq. 0) then
         phi=rlu(0)*twopi
         phiq=phi
         return

c  cahn cos\phi

      elseif(iaz .eq. 1) then
         cahn=4.*sqrt(1.-y)*(2.-y)/(1.+(1-y)**2)*pt/sqrt(q2)
         cahnmax=1.+cahn
 1       continue
         phi=rlu(0)*twopi
         dstr=1.-cahn*cos(phi)
         if(dstr/cahnmax .lt. rlu(0)) goto 1
         phiq=phi

c  cahn cos\phi + cos2\phi

      elseif(iaz .eq. 2) then
         q2pt  =  pt/sqrt(q2)
         q2pt2 =  pt**2/q2
         y1  =  (1.-y)
         y11 =  1.+y1**2+4.*y1*q2pt2
         y12 =  sqrt(y1)
         y2  =  2.-y
         cahn1 = 4.*y12*y2*q2pt
         cahn2 = 4.*y1*q2pt2
	 cahnmax=y11+cahn1+2.*cahn2

 2       continue
         phi=rlu(0)*twopi
         dstr=y11-cahn1*cos(phi)+cahn2*cos(2.*phi)
         if(dstr/cahnmax .lt. rlu(0)) goto 2
         phiq=phi
         return

c  sivers

      elseif(iaz .eq. 3) then
         kfq=k(lst(26),2)       !  outgoing quark type
         lst(29)=0
         call lframe(3,0)
         phimu=rlu(0)*twopi
         call ludbrb(0,0,0.,phimu,0.d0,0.d0,0.d0)
         do i=1,3
            pbeam(i)=p(1,i)
            pvg(i)=p(3,i)
         enddo
         call spinaz(phis,st)
         call ludbrb(0,0,0.,-phimu,0.d0,0.d0,0.d0)
         call lframe(1,0)
         lst(29)=1
c         call siv_ansel(pt,kfq,phis,st,phi,sivq)
         call siv_examp(pt,kfq,phis,st,phi,sivq)
c         call siv_examp1(pt,kfq,phis,st,phi,sivq)
         phiq=phi
         return

c  sivers+cahn

      elseif(iaz .eq. 4) then
         kfq=k(lst(26),2)       !  outgoing quark type
         lst(29)=0
         call lframe(3,0)
         phimu=rlu(0)*twopi
         call ludbrb(0,0,0.,phimu,0.d0,0.d0,0.d0)
         do i=1,3
            pbeam(i)=p(1,i)
            pvg(i)=p(3,i)
         enddo
         call spinaz(phis,st)
         call ludbrb(0,0,0.,-phimu,0.d0,0.d0,0.d0)
         call lframe(1,0)
         lst(29)=1
         call siv_cahn(pt,kfq,phis,st,phi,sivq)
         phiq=phi
         return
      endif
c     
      end

***********************************************************

      subroutine spinaz(phis,st)
c  calculates azimuth and transverse component of spin     
      common/twopi/pi,twopi
      common/labmom/pbeam(3),pvg(3)
      dimension vng(3),vns(3)
      common/spin/spin(3)
      common/lujets/n,k(4000,5),p(4000,5),v(4000,5)
      call cross(pbeam,pvg,vng)
      call cross(spin,pvg,vns)
      cosphi=vdotn(vng,vns,3)
      phis=acos(cosphi)
      soutpl=vdot(spin,vng,3)
      if(soutpl.gt.0.) phis=twopi-phis
      costetsgam=vdotn(spin,pvg,3)
      st=sqrt(1.-costetsgam**2)
      return
      end

**********************************************************************

      subroutine siv_ansel(pt,kfq,phis,st,phi,sivq)
      common/twopi/pi,twopi
      common /leptou/ cut(14),lst(40),parl(30),x,y,w2,q2,u
      if(kfq.ne.1.and.kfq.ne.2)then
         phi=rlu(0)*twopi
         sivq=mod(phi-phis+twopi,twopi)
      elseif(kfq.eq.2)then
         a=2.
         b=0.3
         ab=a+b
         c=0.5
         siv=c*x**a*(1.-x)**b*ab**ab/a**a/b**b
         bet=1.25
         r=0.7
         rr=(1.-r)/r
         siv=siv*bet*sqrt(exp(1.)*rr)*pt*exp(-bet**2*rr*pt**2)
         siv=siv*st
         sivmax=1.+abs(siv)
 2       continue
         phi=rlu(0)*twopi
         sivq=mod(phi-phis+twopi,twopi)
         dstr=1.+siv*sin(sivq)
         if(dstr/sivmax .lt. rlu(0)) goto 2
      elseif(kfq.eq.1)then
         a=1.5
         b=0.2
         ab=a+b
         c=-1.
         siv=c*x**a*(1.-x)**b*ab**ab/a**a/b**b
         bet=1.25
         r=0.7
         rr=(1.-r)/r
         siv=siv*bet*sqrt(2.*exp(1.)*rr)*pt*exp(-bet**2*rr*pt**2)
         siv=siv*st
         sivmax=1.+abs(siv)
 3       continue
         phi=rlu(0)*twopi
         sivq=mod(phi-phis+twopi,twopi)
         dstr=1.+siv*sin(sivq)
         if(dstr/sivmax .lt. rlu(0)) goto 3
      endif
      return
      end

**********************************************************************

      subroutine siv_examp(pt,kfq,phis,st,phi,sivq)
      common/twopi/pi,twopi
      common /leptou/ cut(14),lst(40),parl(30),x,y,w2,q2,u
      common/sfun/cn(4),a(4),b(4),cm1,iqs
      if(abs(kfq).eq.3)then   !strange sea
         phi=rlu(0)*twopi
         sivq=mod(phi-phis+twopi,twopi)
      elseif(kfq.eq.2)then
         if(lst(27).eq.0)then
            iqs=1
            siv=cn(iqs)
            siv=siv*st
            siv=siv*sivx(x)
            siv=siv*sivpt(pt)
            sivmax=1.+abs(siv)
 1          continue
            phi=rlu(0)*twopi
            sivq=mod(phi-phis+twopi,twopi)
            dstr=1.+siv*sin(sivq)
            if(dstr/sivmax .lt. rlu(0)) goto 1
         elseif(lst(27).eq.1)then
            iqs=3
            siv=cn(iqs)
            siv=siv*st
            siv=siv*sivx(x)
            siv=siv*sivpt(pt)
            sivmax=1.+abs(siv)
 3          continue
            phi=rlu(0)*twopi
            sivq=mod(phi-phis+twopi,twopi)
            dstr=1.+siv*sin(sivq)
            if(dstr/sivmax .lt. rlu(0)) goto 3
         endif
      elseif(kfq.eq.1)then
         if(lst(27).eq.0)then
            iqs=2
            siv=cn(iqs)
            siv=siv*st
            siv=siv*sivx(x)
            siv=siv*sivpt(pt)
            sivmax=1.+abs(siv)
 2          continue
            phi=rlu(0)*twopi
            sivq=mod(phi-phis+twopi,twopi)
            dstr=1.+siv*sin(sivq)
            if(dstr/sivmax .lt. rlu(0)) goto 2
         elseif(lst(27).eq.1)then
            iqs=4
            siv=cn(iqs)
            siv=siv*st
            siv=siv*sivx(x)
            siv=siv*sivpt(pt)
            sivmax=1.+abs(siv)
 4          continue
            phi=rlu(0)*twopi
            sivq=mod(phi-phis+twopi,twopi)
            dstr=1.+siv*sin(sivq)
            if(dstr/sivmax .lt. rlu(0)) goto 4
         endif
      endif
      return
      end

**********************************************************************

      function sivx(x)
      common/sfun/cn(4),ai(4),bi(4),cm1,iqs
      a=ai(iqs)
      b=bi(iqs)
      c=(a+b)**(a+b)/a**a/b**b
      sivx=c*x**a*(1.-x)**b
      return
      end

**********************************************************************

      function sivpt(pt)
      common/sfun/cn(4),a(4),b(4),cm1,iqs
cold      sivpt=1.05*pt/(pt**2+0.3)
      ppt=pt/cm1
      sivpt=2.331643982*ppt*exp(-ppt**2)
      return
      end

**********************************************************************

      subroutine siv_examp1(pt,kfq,phis,st,phi,sivq)
      common/twopi/pi,twopi
      common /leptou/ cut(14),lst(40),parl(30),x,y,w2,q2,u
      siv=200.*x*pt
      siv=siv*st
      sivmax=1.+abs(siv)
 1    continue
      phi=rlu(0)*twopi
      sivq=mod(phi-phis+twopi,twopi)
      dstr=1.+siv*sin(sivq)
      if(dstr/sivmax .lt. rlu(0)) goto 1
      return
      end

**********************************************************************

      subroutine siv_cahn(pt,kfq,phis,st,phi,sivq)
c   just an example, not well tested
      common/twopi/pi,twopi
      common /leptou/ cut(14),lst(40),parl(30),x,y,w2,q2,u
      siv=0.5
      cahn=-0.5
      siv=siv*st
      phi0=atan(siv*cos(phis)/(cahn-siv*sin(phis)))
      sivmax=1.+abs(siv*sin(phi0-phis)+cahn*cos(phi0))
 1    continue
      phi=rlu(0)*twopi
      sivq=mod(phi-phis+twopi,twopi)
      dstr=1.+siv*sin(sivq)+cahn*cos(phi)
      if(dstr/sivmax .lt. rlu(0)) goto 1
      return
      end

**********************************************************************
