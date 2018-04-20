      character*10 cl_target
      REAL cl_beam_energy,cl_zpos,cl_zwidth,cl_trad,cl_zmin,cl_rast
      REAL  cl_emin,cl_emax,cl_tmin,cl_tmax
      REAL  cl_xmin,cl_xmax,cl_qmin,cl_qmax,cl_wmin,cl_wmax
      REAL  rdtiv,tiv,tetacm,alf05
      integer cl_verblev,cl_rand,cl_nprint,cl_mstu
      integer cl_pid,cl_triggers,cl_pol,cl_mstp,cl_vmstp
      integer  cl_smear
      integer  cl_ht,cl_outform,cl_printid,cl_lst15,cl_lst37,cl_lst8
      real  cl_sma,cl_smb,cl_smc,cl_smd,cl_parl3
      real  cl_parj12,cl_parj2,cl_parj3,cl_parj21,cl_parj23,cl_parj24
     6,cl_parj33,cl_parj41,cl_parj42,cl_parj4,cl_parj5,cl_parj6,cl_parj7
      logical clasdisOK,ntOK,cl_radgen,cl_radtab,datfileOK
      logical cl_her,acce16,acceg1
      character*70 bosout,cl_path
      character*180 datfilename
      common /OUT_NAMES/ bosout
      common /claspepsi1/ cl_verblev,cl_nprint,cl_rand,cl_mstu
     6   ,cl_mstp,cl_vmstp,cl_ht,cl_printid,cl_lst15,cl_lst37,cl_lst8
      common /claspepsi2/ cl_target
      common /claspepsi3/ cl_triggers,cl_pid,cl_beam_energy,cl_pol
      common /claspepsi4/ cl_zpos,cl_zwidth,cl_rast
      common /claspepsi44/ cl_trad,cl_zmin,cl_parl3
      common /claspepsi5/ cl_emin,cl_emax,cl_tmin,cl_tmax
      common /claspepsi55/ cl_parj2,cl_parj3,cl_parj4,cl_parj5,cl_parj6
     6     ,cl_parj7,cl_parj12,cl_parj21,cl_parj23,cl_parj24,cl_parj33
     6     ,cl_parj41,cl_parj42
      common /claspepsi6/ cl_smear,cl_sma,cl_smb,cl_smc,cl_smd
      common /claspepsi6/ cl_xmin,cl_xmax,cl_qmin,cl_qmax,
     6cl_wmin,cl_wmax
      common /clasdisrad/ cl_radgen,cl_radtab
      common /clasdiscont1/ clasdisOK,ntOK,cl_her,acce16,acceg1
      common /clasdisfile/ datfileOK,datfilename,cl_path,cl_outform
