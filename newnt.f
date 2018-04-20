       Subroutine Newnt 
* * Example of a COMIS subroutine to create a Ntuple interactively. 
* Data is read from a text input file * 
       character*8 mother,in1,in2 
       common/ntupc/mother,in1,in2 
       common/ntupr/xover 
* 
       lin=41 
       lout=42 
       id=1 
       open(unit=lin,file='datafile.dat',status='old') 
       call hropen(lout,'NTUPLE','New_Ntuple.hbook','N',1024,istat) 
* 
       call hbnt(id,'New Ntuple',' ') 
       call hbname(id,'ntupr',xover,'XOVER') 
       call hbnamc(id,'ntupc',mother,'MOTHER:c*8,in1:c*8,in2:c*8') 
* 
10     read(lin,1000,end=20,err=20)xover,mother,in1,in2 1000 
       format(e15.7,2x,a,7x,a,7x,a) 
       call hfnt(1) 
       go to 10 
* 
20     call hrout(id,icycle,' ') 
       call hrend('NTUPLE') 
       close (lin) 
       close (lout) 
       end 


