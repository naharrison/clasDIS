      integer function accep12(torus)
      include 'mc22.inc'
       real newphi
       pi=acos(-1.0)
       torcur=torus
       phis=newphi(gelef*180.0/pi)
       thetad=gelet*180.0/pi
       p=gelee
       call clas_at12g(11,p,thetad,phis,torcur,d_phi,weight)
c       print *,thetad,gelef,phis,weight
       accep12=weight
       end
*     
c==================================================
      subroutine clas_at12g(id,p,thetad,phis,torcur,d_phi,weight)
*     
*     Version 1.0 / November 30, 1999
* Inputs -
*     id:     Particle ID according to PDG (11-electron, 2212-proton, 
*     211-pi+, -211-pi-, 22-photon ...)
*     p:      momentum, in GeV, on output will be smeared
*     phis:   phi in the sector in degrees (at mid plane phi=0.), 
*     thetad: scattering angle in degrees at dc layer the DC layer 1 
*     torcur: torus current (e.g. -1500.)
* Outputs -
*     d_phi:  part of the phi acceptance, full acceptance should be 
*     delta_phi=2*pi*d_phi
*     weight: =1. if accepted,  0. if not
c     
      implicit none
      real p,thetad,phis,torcur,d_phi,weight
      integer id
c     
      real phie_mod,thcut,pnorm,expon
      logical pcutl
      real ec_r,Ec_th,rl,rh
      data Ec_r,ec_th,rl,rh/510.32,0.436332,562.22,215.67/
      real ylow,yhi,tgrho,dl,rc0
      data ylow,yhi,tgrho,dl,rc0/-182.974,189.956,1.95325,100.,120./
      real r,xh,yh,zh,xcut
      real uh,vh,wh,xi,yi
      real sinth,costh,sinph,cosph,dth
      real rcoil,hcoil,dcoil
      real thetamin, thetaf_clas, dgap
      data thetamin, thetaf_clas,dgap/5.,45.,1./
c     
      real d2r,pi,r2d,rran(3)
      data d2r,pi,r2d/0.01745,3.1415926,57.299578/
      integer lrun
      data lrun/3/
c     
      weight=0.
      d_phi=0.
      if(thetad.lt.thetamin)return
      if(id.eq.11.and.p.lt.0.5)return
      if(abs(id).eq.211.and.p.lt.0.2)return
      call rnorml(rran,lrun)
c     
      IF(id.eq.11.or.id.eq.-211)THEN
*     Define acceptance cuts for electrons
         if(thetad.gt.thetamin.and.thetad.lt.thetaf_clas)then
            If(torcur.lt.0.)then
               thcut = 5.0  - 25.*(1.-p*1500./abs(torcur)/5.)**24
               pnorm = p*3375./abs(torcur)
               pcutl=pnorm.gt.(.6+4.*(thetad/45.)**6)*abs(torcur)/3375.
               if(thetad.gt.thcut.and.pcutl)then
                  d_phi=25.*sin((thetad-5.)*d2r)**0.15* 
     *                 (p*1500./abs(torcur)/5.)**(1./24.)
                  if(abs(phis).lt.d_phi)then 
                     weight=1.
                  else
                     weight=0.
                  endif
               else
                  weight=0.
                  d_phi=0.
               endif
            ElseIf(torcur.gt.0.)then
               thcut =5.5 + 17.5/((p+0.050)*3375./abs(torcur))
               if(thetad.gt.thcut)then
                  expon = 0.25*(p*3375./torcur)**0.33
                  d_phi = 27.*sin((thetad-thcut)*d2r)**expon 
                  if(abs(phis).lt.d_phi)then
                     weight=1.
                  else
                     weight=0.
                  endif
               else
                  d_phi= 0.
                  weight=0.
               endif
            Endif
            p=p+rran(1)*p*0.005
            thetad=thetad+rran(2)/d2r/1000.
            phis=phis+rran(3)/d2r/500.
         elseif(thetad.gt.thetaf_clas+dgap.and.id.eq.-211)then
            weight=1.
            d_phi=30.
            p=p+rran(1)*p*0.02
         endif
c     
         d_phi=d_phi/30.
c     End of simulation for negatives
      ELSEIF(id.eq.2212.or.id.eq.211.or.id.eq.45)THEN
c     Positives come here
         if(thetad.gt.thetaf_clas.and.thetad.lt.thetaf_clas+dgap)return
*     minimum momentum cut
         if(p.gt.0.4.and.thetad.lt.thetaf_clas)then
            If(torcur.lt.0.)then
               thcut = 5.0  - 25.*(1.-p*1500./abs(torcur)/5.)**24
               pnorm = p*3375./abs(torcur)
               pcutl=pnorm.gt.(.6+4.*(thetad/45.)**6)*abs(torcur)/3375.
               if(thetad.gt.thcut.and.pcutl)then
                  d_phi=25.*sin((thetad-5.)*d2r)**0.15* 
     *                 (p*1500./abs(torcur)/5.)**(1./24.)
               else
                  d_phi=0.
               endif
            ElseIf(torcur.gt.0.)then
               d_phi = 28.*cos((thetad-37.)*d2r)*
     *              (p*3375./abs(torcur)/5.)**(1./64.)
            EndIf
         elseif(thetad.gt.thetaf_clas+dgap.and.thetad.lt.100.
     &           .and.p.gt.0.2.and.p.lt.2.)then
            d_phi=30.
         endif 
         if(abs(phis).lt.d_phi)then 
            weight=1.
         else
            weight=0.
         endif
         d_phi=d_phi/30.
         if(thetad.lt.thetaf_clas)then
            p=p+rran(1)*p*0.007
         else
            p=p+rran(1)*p*0.02
         endif
         thetad=thetad+rran(2)/d2r/500.
         phis=phis+rran(3)/d2r/250.
c     End of simulation for positives
      ELSEIF(id.eq.22)THEN
c     Photon
         if(p.lt.0.1)return
         dth=atan(2./r/sqrt(p))*r2d
         if(thetad.gt.thetaf_clas+dgap.and.thetad.lt.100.)then
            d_phi=30.
            weight=1.
            p=p+0.05*rran(1)*sqrt(p)
            thetad=thetad+rran(2)*dth
            phis=phis+rran(3)*dth
         elseif(thetad.lt.thetaf_clas)then
            if(p.lt.0.2)return
            sinth=sin(thetad*d2r)
            costh=cos(thetad*d2r)
            sinph=sin(phis*d2r)
            cosph=cos(phis*d2r)
c     
            r=(ec_r*(sin(ec_th)**2+cos(ec_th))+dl)/
     /           (costh+sinth*cosph*sin(ec_th))
            xh=r*cosph*sinth
            yh=r*sinph*sinth
            zh=r*costh-dl
c     
            call ec_xyz_duvw(xh,yh,zh,uh,vh,wh,xi,yi)
c     
            rcoil=dl*costh+sqrt(rc0**2-(dl*sinth)**2)
            hcoil=rcoil*sinth*cosph
            dcoil=hcoil*sin(0.5236)/cos(0.5236)-abs(rcoil*sinth*sinph)
c     
            if(yi.gt.ylow+10..and.yi.lt.yhi-5.)then
               xcut=(yi-ylow)/tgrho
               if(abs(xi).lt.xcut-5.)then
                  weight=1
                  p=p+0.1*rran(1)*sqrt(p)
               elseif(abs(xi).gt.xcut.and.dcoil.lt.8.)then
                  weight=1
                  p=p+0.05*rran(1)*sqrt(p)
               endif
               thetad=thetad+rran(2)*dth
               phis=phis+rran(3)*dth
               d_phi=1.-5./(abs(r*sinph*sinth)+8.)
            endif
         endif
      ELSE
c     No acceptance and smearing for this particle
         weight=0.
         d_phi=0.
      ENDIF
*     
      return
      end
C=======================================================================
      Subroutine ec_xyz_duvw(x,y,z,u,v,w,xi,yi)
      REAL x,y,z,u,v,w,xi,yi,zi
      REAL EC_the,ec_phi,phi,ylow,yhi,tgrho,sinrho,cosrho
      data EC_the/0.4363323/
      data ylow,yhi/-182.974,189.956/
      data tgrho,sinrho,cosrho/1.95325,0.8901256,0.455715/
      real rot(3,3)
c
      phi=atan2(y,x)*57.29578
      if(phi.lt.0.)phi=phi+360.
      phi=phi+30.
      if(phi.ge.360.)phi=phi-360. 	
      Ec_phi=int(phi/60.)*1.0471975
c     
      rot(1,1)=cos(Ec_the)*cos(Ec_phi)
      rot(1,2)=-sin(Ec_phi)
      rot(1,3)=sin(Ec_the)*cos(Ec_phi)
      rot(2,1)=cos(Ec_the)*sin(Ec_phi)
      rot(2,2)=cos(Ec_phi)
      rot(2,3)=sin(Ec_the)*sin(Ec_phi)
      rot(3,1)=-sin(Ec_the)
      rot(3,2)=0.
      rot(3,3)=cos(Ec_the)
c     
      yi=x*rot(1,1)+y*rot(2,1)+z*rot(3,1)
      xi=x*rot(1,2)+y*rot(2,2)+z*rot(3,2)
      zi=x*rot(1,3)+y*rot(2,3)+z*rot(3,3)
      zi=zi-510.32
      u=(yi-ylow)/sinrho
      v=(yhi-ylow)/tgrho-xi+(yhi-yi)/tgrho
      w=((yhi-ylow)/tgrho+xi+(yhi-yi)/tgrho)/2./cosrho
      end

C
      real Function newphi(phi)
      real phi,phinew
      if (phi.gt.330.) then
        phinew = phi-360.
      elseif (phi.ge.0.0.and.phi.le.30.) then
        phinew = phi
      elseif (phi.gt.30.0.and.phi.le.90.) then
        phinew = phi-60.
      elseif (phi.gt.90.0.and.phi.le.150.) then
        phinew = phi-120.
      elseif (phi.gt.150.0.and.phi.le.210.) then
        phinew = phi-180.
      elseif (phi.gt.210.0.and.phi.le.270.) then
        phinew = phi-240.
      elseif (phi.gt.270.0.and.phi.le.330.) then
        phinew = phi-300.
      endif
c       print *,phi,phinew
      newphi=phinew
      end  

