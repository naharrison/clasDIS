c **********************************************************************

      subroutine lepto

      implicit none

c...administer the generation of an event.
c...note: if error flag lst(21) is non-zero, no proper event generated.

      common /lintrl/ psave(3,4,5),ksave(4),xmin,xmax,ymin,ymax,
     &q2min,q2max,w2min,w2max,ilep,inu,ig,iz
      real psave,xmin,xmax,ymin,ymax,q2min,q2max,w2min,w2max
      integer ksave,ilep,inu,ig,iz
      save /lintrl/

*
* to avoid variable conflictions, a second keep element is necessary
* with the same common block name (see lptou2)
*
      common /leptou/ cut(14),lst(40),parl(30),
     &                x,y,w2,q2,u
      real cut,parl,x,y,w2,q2,u
      integer lst
      save /leptou/

      common /linter/ pari(50),ewqc(2,2,8),qc(8),zl(2,4),zq(2,8),pq(17)
      real  pari,ewqc,qc,zl,zq,pq
      save /linter/

      integer  nlupdm,nplbuf
      parameter  (nlupdm=4000,nplbuf=5)
      common/lujets/n,k(nlupdm,5),p(nlupdm,nplbuf),v(nlupdm,5)
      integer  n,k
      real  p,v
      save /lujets/

      common/ludat1/mstu(200),paru(200),mstj(200),parj(200)
      integer  mstu,mstj
      real  paru,parj
      save /ludat1/

      common /lboost/ dbeta(2,3),stheta(2),sphi(2),pb(5),phir
      double precision dbeta
      real stheta,sphi,pb,phir
      save /lboost/


* if ariadne is used to simulate the parton shower evolution, the
* ardat1 common block is neccessary for a proper interface.
* if ariadne is used to simulate the parton shower evolution, the
* ardat1 common block is neccessary for a proper interface.
      common /ardat1/ para(40),msta(40)
      real para
      integer msta
      save /ardat1/ 


      integer nummis,nwarn,i,j,ns,l17
      real ulmass,ulangl,rlu,ulalps
      real qg,qqb,spq,srlu,plu,pt,phi,energy,p2,atan
      double precision detot,dari29,dari30
      dimension spq(17)
      data nummis,nwarn/0,10/,dari29,dari30/2*0.d0/

      l17=0
    1 lst(21)=0
      do 10 i=1,10
      do 10 j=1,5
      k(i,j)=0
   10 v(i,j)=0.
      do 15 i=1,4
      k(i,1)=21
   15 k(i,2)=ksave(i)
      k(4,1)=1
      n=2

      if(lst(17).ne.0.and.lst(2).gt.0) then
c...lepton and/or nucleon energy may vary from event to event,
        if(l17.eq.0) then
c...momentum vectors from p(i,j) i=1,2 j=1,2,3 on entry in lepto
          do 20 i=1,2
          p(i,5)=ulmass(k(i,2))
          p(i,4)=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
          do 20 j=1,5
   20     psave(3,i,j)=p(i,j)
        else
c...momentum vectors from psave if new try, i.e. jump back to 1
          do 25 i=1,2
          do 25 j=1,5
   25     p(i,j)=psave(3,i,j)
        endif
        l17=1
c...transform to cms of incoming particles, lepton along +z axis.
        do 30 j=1,3
   30   dbeta(1,j)=(dble(p(1,j))+dble(p(2,j)))/
     &             (dble(p(1,4))+dble(p(2,4)))
        call ludbrb(0,0,0.,0.,-dbeta(1,1),-dbeta(1,2),-dbeta(1,3))
        sphi(1)=ulangl(p(1,1),p(1,2))
        call ludbrb(0,0,0.,-sphi(1),0.d0,0.d0,0.d0)
        stheta(1)=ulangl(p(1,3),p(1,1))
        call ludbrb(0,0,-stheta(1),0.,0.d0,0.d0,0.d0)
        lst(28)=2
        parl(21)=2.*(p(1,4)*p(2,4)-p(1,3)*p(2,3))
      else
c...initial state momenta fixed from linit call.
        do 42 i=1,2
        do 40 j=1,5
   40   p(i,j)=psave(3,i,j)
   42   if(psave(3,1,3).lt.0.) p(i,3)=-psave(3,i,3)
        lst(28)=3
      endif
      call leptox
c...return if error or if no event to be generated.
      if(lst(21).ne.0.or.lst(2).le.0.or.lst(7).eq.-1) return

      if(pari(29).lt.0.5) then
c...for first call, reset double precision counters.
        dari29=0.d0
        dari30=0.d0
      endif
      dari29=dari29+1.d0
      pari(29)=dari29
c     call gulist(-3,2)
c...scattered lepton and exchanged boson added to event record in lkinem
c...transform to lepton-nucleon cms if not made earlier
      if(lst(17).eq.0) then
        do 46 i=3,4
        do 45 j=1,5
   45   psave(3,i,j)=p(i,j)
   46   if(psave(3,1,3).lt.0.) psave(3,i,3)=-p(i,3)
        call ludbrb(0,0,0.,0.,0.d0,0.d0,-dbeta(1,3))
        lst(28)=2
      endif
      do 50 i=1,4
      do 50 j=1,5
   50 psave(2,i,j)=p(i,j)
c     call gulist(-2,2)

c...prepare for parton cascade.
      if(lst(8).ge.2.and.mod(lst(8),10).ne.9) call lshowr(0)

c...transform to hadronic cms, boost parameters in double precision.
      detot=dble(p(1,4))-dble(p(4,4))+dble(p(2,4))
      dbeta(2,1)=-dble(p(4,1))/detot
      dbeta(2,2)=-dble(p(4,2))/detot
      dbeta(2,3)=(dble(p(1,3))-dble(p(4,3))+dble(p(2,3)))/detot
      call ludbrb(0,0,0.,0.,-dbeta(2,1),-dbeta(2,2),-dbeta(2,3))
      sphi(2)=0.
      stheta(2)=ulangl(p(3,3),p(3,1))
      call ludbrb(0,0,-stheta(2),0.,0.d0,0.d0,0.d0)
      lst(28)=1
      do 60 i=1,4
      do 60 j=1,5
   60 psave(1,i,j)=p(i,j)
c...save momentum of exchanged boson (used in subroutine lframe).
      do 70 j=1,5
   70 pb(j)=p(3,j)
c     call gulist(-1,2)

   90 n=4
      mstu(1)=n+1
      lst(26)=n+1
      lst(27)=0
      parl(25)=ulalps(q2)
      if(lst(8).eq.1.or.lst(8)/10.eq.1.or.mod(lst(8),10).eq.9) then
c...probabilities for hard, first order qcd events.
cae...corrected what to do when lqgev or lqqbev fail. now make lqev.
         call lqcdpr(qg,qqb)

        do 100 i=1,17
  100   spq(i)=pq(i)
  200   srlu=rlu(0)

        if(srlu.gt.qqb+qg) then
          call lqev
        elseif(srlu.gt.qqb) then
          if(lst(8).eq.9) then
            do 211 i=1,17
  211         pq(i)=spq(i)
            call lqev
          else
            call lqgev
            do 212 i=1,17
              pq(i)=spq(i)
  212       continue
          endif
        else
          call lqqbev
          do 213 i=1,17
            pq(i)=spq(i)
  213     continue
          if(lst(8).eq.9.and.lst(21).eq.0) then
            if(plu(5,11).lt.q2*para(20)) then
              do 220 i=1,17
  220         pq(i)=spq(i)
              call lqevar(k(5,2),k(7,2))
            endif
          endif
        endif
        if(lst(21).ne.0) then
  230     call lqev
          if(lst(21).ne.0) goto 230
        endif   
      else
c...qpm model without qcd corrections (cascade applied later).
  300   call lqev
        if(lst(21).ne.0) goto 300
      endif

      ns=mstu(1)
      mstu(1)=0
c     call gulist(-3,2)
c     write(6,*) ' lst(24)=',lst(24)
cjr--       no preclustering of small systems
          mstj(14)=-1
cjr--            
      if(lst(8).le.1.or.mod(lst(8),10).eq.9) then
c...no parton cascade, introduce primordial kt.
        if(parl(3).gt.1.e-03) then
          call lprikt(parl(3),pt,phi)
          call ludbrb(ns,n,0.,-phi,0.d0,0.d0,0.d0)
          call ludbrb(ns,n,atan(2.*pt/sqrt(w2)),phi,0.d0,0.d0,0.d0)
        endif
        if(mod(lst(8),10).ne.9) then
c...check system against fragmentation cuts.
          mstu(24)=0
          call luprep(0)
          if(mstu(24).ne.0) then
            if(lst(3).ge.1) write(6,*)'luprep error mstu(24)=',mstu(24),
     &                                ', new event generated'
            lst(21)=1
            goto 1
          endif
        endif
      elseif(lst(24).eq.1) then
c...include parton cascades (+ remnant & kt) on q-event
        call lshowr(1)
      else
c...include parton cascades (+ remnant & kt) on qg- or qqbar-event
        call lmeps
      endif
      if(lst(21).ne.0) then
c        if(lst(3).ge.1)
c     &     write(6,*)'cascade error lst(21)= ',lst(21),
c     &               ', new event generated'
        goto 1
      endif

cjr--       soft colour interactions
      if(lst(34).eq.1) call lsci(parl(7))
      if(lst(21).ne.0) goto 1
cjr--       take care of small systems
      call lsmall
      if(lst(21).ne.0) then
        if(lst(3).ge.1) write(6,*)' lsmall error lst(21)= ',lst(21),
     &                            ', new event generated'
        goto 1
      endif
      mstj(14)=1
      call luprep(0)
      if(mstu(24).ne.0) then
         if(lst(3).ge.1) write(6,*)' luprep error mstu(24)= ',mstu(24),
     &                             ', new event generated'
         lst(21)=1
      endif
cjr--            
      if(lst(21).ne.0) goto 1

      do 400 i=1,n
c...correct energy-momentum-mass mismatch for real particle
      if(p(i,5).lt.0.) goto 400
      energy=sqrt(dble(p(i,5))**2+dble(p(i,1))**2+dble(p(i,2))**2+
     &dble(p(i,3))**2)
      p2=dble(p(i,4))**2-dble(p(i,1))**2-dble(p(i,2))**2-dble(p(i,3))**2
      if(abs(energy-p(i,4))/(psave(3,1,4)+psave(3,2,4)).gt.paru(11))then
        nummis=nummis+1
c...for testing purposes
c       if(lst(3).ge.1.and.nummis.le.nwarn) then
c         write(6,1000) i,(k(i,j),j=1,2),(p(i,j),j=1,5),
c    &    sign(sqrt(abs(p2)),p2),energy,int(dari29),nwarn
c         if(abs(p2-p(i,5)**2).gt.400.) call lulist(2)
c       endif
cae        write(6,*) 'energy mismatch',lst(24),parl(28),parl(29),nummis
        goto 90
      endif
      p(i,4)=energy
  400 continue

      dari30=dari30+1.d0
      pari(30)=dari30
ctest if(lst(23).eq.2) parl(24)=parl(24)*dari30/dari29

      do 500 i=1,n
      do 500 j=1,5
  500 v(i,j)=0.
      if(lst(7).eq.1) then
        call luexec
        if(mstu(24).ne.0) then
          write(6,*) ' error from jetset, new event made'
          goto 90
        endif
      endif

c     call gulist(-1,2)
c...transform to desired frame
c     lst(28)=1
      lst(29)=0
      phir=6.2832*rlu(0)
      if(lst(17).eq.0) then
        if(lst(5).ge.2) call lframe(lst(5),0)
c...restore momenta (e,p,boson,l) due to numerical errors from boosts
        do 600 i=1,4
        do 600 j=1,5
  600   p(i,j)=psave(lst(28),i,j)
        if(lst(6).eq.1.and.lst(28).ge.2) then
c...random rotation in azimuthal angle
          call ludbrb(0,0,0.,phir,0.d0,0.d0,0.d0)
          lst(29)=1
        endif
      else
        if(lst(5).ge.2) call lframe(lst(5),lst(6))
      endif
c...deactivate scattered lepton
      if(mod(lst(4),10).eq.0) k(4,1)=21
c     call gulist(0,2)

      return
 1000 format(' warning: too large numerical mismatch in ',
     &'particle energy-momentum-mass',
     &/,3x,'i k(i,1) ..2)  p(i,1)  p(i,2)  p(i,3)',
     &'  p(i,4)  p(i,5)    mass  energy',/,i4,2i6,7f8.3,/,
     &' event no.',i8,' regenerated. only first',i5,' warnings printed')
      end

