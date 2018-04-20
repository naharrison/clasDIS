
C***********************************************************************
          subroutine write2filep(jj,il)
C***********************************************************************
c write the LUND event with limited to cl_printid-degrees events:
c
c
          implicit none
#include "names.inc"
#include "claslepto.inc"
#include "ntupgkoko.inc"
          integer ifile,itotal,iparl23
          DATA ifile /0/
          DATA itotal /0/
c          character*100  fname
          character*200  fname
c
          real plu,z_pos,x_pos, y_pos, random_num,random_norm, trad
          real thetap
          integer jj,il,i,j,istat,elpol
          logical lFirst
          DATA lFirst /.TRUE./
          DATA elpol /1/

          integer jc,m,lund2geantid,gamp_id,icharge

c
          iparl23=PARL(23)
          if(itotal.ge.IPARL23) then   ! each file correspond to number ~ xsection
           close (41)
           lFirst=.TRUE.
           ifile=ifile+1
           itotal=0
          endif



          if (lFirst) then   ! open the file
             call get_datfilenamep()

             jc = 0
             do m=1,len(datfilename)
                if ((datfilename(m:m).eq.' ').and.(jc.eq.0)) then
                   jc=m-1
                endif
             enddo
             if (jc.eq.0) jc=len(datfilename)

           if(ifile.le.9) then
             write(fname,'(A,i1,A)') datfilename(1:jc),ifile,'.dat'
             write(6,'(A,i1,A)') datfilename(1:jc),ifile,'.dat'
           else if(ifile.ge.10.and.ifile.le.99) then
            write(fname,'(A,i2,A)') datfilename(1:jc),ifile,'.dat'
            write(6,'(A,i2,A)') datfilename(1:jc),ifile,'.dat'
           else if(ifile.ge.100.and.ifile.le.999) then
            write(fname,'(A,i3,A)') datfilename(1:jc),ifile,'.dat'
            write(6,'(A,i3,A)') datfilename(1:jc),ifile,'.dat'
           else
            write(fname,'(A,i4,A)') datfilename(1:jc),ifile,'.dat'
            write(6,'(A,i4,A)') datfilename(1:jc),ifile,'.dat'
           endif 


             OPEN(unit=41,file=fname,status='unknown')
             lFirst=.FALSE.
          endif
c
                itotal=itotal+1


            if (lFirst) then   ! open the file
             OPEN(unit=41,file=datfilename,status='unknown')
             lFirst=.FALSE.
            endif
c
c           only exclusive pi+
c
c            if(.not.(gpipe+gneue+gelee.le.cl_beam_energy))return
c
c
c
                 j=0
                 do i=1,N
        thetap=acos(p(i,3)/sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2))/0.01745
                  if(k(i,1).lt.11.and.thetap.gt.cl_tmin.and.
     6             P(i,5).gt.-100.and.thetap.lt.1.0*cl_printid) j=j+1   ! save space
                 enddo
c23456789012345678901234567890123456789012345678901234567890123456789012  
               if(il.eq.2) then
                write(41,111) j,PARL(1),radQ2,radnu,1.0*cl_pol,1.0*elpol
               else
                write(41,101) j,PARL(1),PARL(2),cl_pol,elpol,
     6                   LEX,LEYY,LEW2,radq2,radnu
               endif
c
 111           FORMAT(2x,I5,10F10.4)
 110           FORMAT(2x,10I10)
 112           FORMAT(2x,10F10.4)
 113           FORMAT(2x,2I4,10F10.5)
c 101      FORMAT(2x,I10,2F4.0,2I3,5F6.3)
 101      FORMAT(2x,I10,2F4.0,2I3,5(F6.3,2x))
 102      FORMAT(2x,I3,F4.0,2I5,2I3,5F10.4,2X,5F10.4)
c
c        what is this?
c
c          z_pos=cl_zpos+(2.*random_num()-1)*cl_zwidth
c
c 120                x_pos=0.1*random_norm()
c                    if(x_pos.GT.cl_trad.OR.x_pos.LT.-cl_trad)goto 120 
c                 
c 121                y_pos=0.1*random_norm()
c                    if(y_pos.GT.cl_trad.OR.y_pos.LT.-cl_trad)goto 121 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

              do i=1,N
               if(abs(k(i,2)).gt.9999) k(i,2)=k(i,2)/100.0
 
              if(il.eq.2.or.il.eq.3) then
                 if(k(i,1).lt.11) then   ! save space
                  if(il.eq.2) then
                   write(41,110) (k(i,j),j=1,2),k(k(i,3),2),i,k(i,3)
                   write(41,112) (p(i,j),j=1,5) 

c23456789012345678901234567890123456789012345678901234567890123456789012  

                  write(41,112) (gelex+V(i,1)*0.1),(geley+V(i,2)
     6                           *0.1),(V(i,3)*0.1+gelez), 0.0,0.0

                 else    ! gamp
              gamp_id=lund2geantid(k(i,2))
              icharge=plu(i,6)
              write(41,113) gamp_id,icharge,(p(i,j),j=1,4)    
                 endif
                 endif 
               else     !(il.eq.2.or.il.eq.3
c23456789012345678901234567890123456789012345678901234567890123456789012 
         thetap=acos(p(i,3)/sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2))/0.01745
                  if(k(i,1).lt.11.and.thetap.gt.cl_tmin.and.
     6             P(i,5).gt.-100.and.thetap.lt.1.0*cl_printid) then
                 write (41,102) i,plu(i,6),(k(i,j),j=1,4),(P(i,j),j=1,5)
     6                                              ,(V(i,j)*0.1,j=1,3)
                 endif
               endif
              enddo
          return
          end
c
c
c
      subroutine get_datfilenamep()
      implicit none
#include "names.inc"
#include "claslepto.inc"
c      character*4  cpol,cdis
      character*4  cpol
      character*5  cdis
      character*10  cxsec
      character*21  ccuts
c
c
      integer jc,m,ncdis
ccc
      print *,'make the data file name',cl_beam_energy
      print *,'************************'
c
      if(cl_beam_energy.lt.10.0) then
      write(ccuts,'(A,F5.3,A,F3.2,A,F3.2)') 'e',cl_beam_energy,
     6'.emin',cl_emin,'tmin',cl_tmin 
      else
      write(ccuts,'(A,F6.3,A,F4.2,A,F3.2)') 'e',cl_beam_energy,
     6'.emn',cl_emin,'tmn',cl_tmin 
      endif
c
      print *,'get the x-section for ',PARL(23)
c
      if(cl_beam_energy.lt.15.0) then
      write(cxsec,'(A,F5.2,A)') '.xs',PARL(23)/1000.0,'nb'
      else
      write(cxsec,'(A,F5.1,A)') '.xs',PARL(23)/1000.0,'nb'
      endif
      print *,'**',ccuts,'++',cxsec
c
      if(cl_pol.eq.1) then
       cpol='.p1.'
      else if (cl_pol.eq.-1) then
       cpol='.m1.'
      else
       cpol='.00.'
      endif
c
      if (cl_pid.eq.0) then 
       cdis='.dis'
       ncdis=4
      else
         if (cl_pid.lt.1000) then
            write (cdis,'(A,I3,A)') '.',cl_pid
            ncdis=4
         else
            write (cdis,'(A,I4,A)') '.',cl_pid
            ncdis=5
         endif
c       write (cdis,'(A,F3.1,I4)') '.',cl_zmin,cl_pid
      endif

      jc = 0
      do m=1,len(cl_path)
         if ((cl_path(m:m).eq.' ').and.(jc.eq.0)) then
            jc=m-1
         endif
      enddo
      if (jc.eq.0) jc=len(cl_path)


      print *,'PARL23,24',PARL(23),PARL(24)
      write(datfilename,'(A,I2,A)') cl_path(1:jc)//'clasdis'//
     6cl_target(1:2)//cpol//ccuts//cxsec//cdis(1:ncdis)//'.'
     6,cl_printid,'.'
      print *,'filename ',datfilename
      end
c
