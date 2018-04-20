
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
      return
      end
