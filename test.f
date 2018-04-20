      subroutine get_datfilename()
       implicit none
#include "names.inc"
      character*4  cpol
      character*21  ccuts
      integer      ibeam
      character*10 cl_target
      character*60 cl_path,fname
      REAL cl_beam_energy,cl_zpos,cl_zwidth,cl_zmin
      REAL  cl_emin,cl_emax,cl_tmin,cl_tmax
      integer cl_triggers,cl_pol,cl_mstp,cl_vmstp
      cl_target='proton'
      cl_path='/home/avakian/w9/'
      cl_pol=-1
      cl_beam_energy=5.735
      cl_emin=0.75           ! def e-min
      cl_tmin=0.24           ! e16->0.3           ! def e'tmin
c
c
c
      write(ccuts,'(A,F4.2,A,F3.2,A,F3.2)') '.e',cl_beam_energy,
     6'.emin',cl_emin,'tmin',cl_tmin 
      print *,'**',ccuts
c
      if(cl_pol.eq.1) then
       cpol='.p1.'
      else if (cl_pol.eq.-1) then
       cpol='.m1.'
      else
       cpol='.00.'
      endif
      write(fname,'(A)') cl_path(1:17)//'clasdis'//
     6cl_target(1:2)//cpol//ccuts//'.dat'
      print *,fname           
      end
