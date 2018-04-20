
C***********************************************************************
          subroutine write2filed()
C***********************************************************************
c write a txt file for studies
c
c
          implicit none
#include "names.inc"
#include "claslepto.inc"
#include "ntupgkoko.inc"
       integer j
       data j/0/
c
*
*--   Enter user code here
        if (j.eq.0) then
        OPEN(unit=41,file='kaon.dat',status='unknown')
        print *,j
        j=1
        endif
        if(gnkap.gt.0) then
c        write(6,'(I1,6F10.4)') gnkap,gelex,geley,gkapt,gkapf,
c     6   sqrt(gkape*gkape-0.24371),gkape/(cl_beam_energy-gelee)
        write(41,'(6F10.4)') gelex,geley,gkapt,gkapf,
     6   sqrt(gkape*gkape-0.24371),gkape/(cl_beam_energy-gelee)
        endif
         return
         end
