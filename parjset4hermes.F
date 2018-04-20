       Subroutine parjset4hermes()
#include "claslepto.inc"
C--------------------------------------
C  set some HERMES jetset values -------
C--------------------------------------
C
C   start MSTJ
       MSTJ(12)=1 ! (D=2) choice of baryon production model. 2->with possib. for diquark to be split
C-MSTJ(40)=0  ! (D=0) regulate time-like parton showers
       MSTJ(45)=4  ! (D=5) maximum flavour that can be produced in shower by $\mathrm{g}\to 
C \mathrm{q}\overline{\mathrm{q}}$; also used to determine  the maximum number of active flavours 
Cin the $\alpha_{\mathrm{s}}$ factor in parton showers (here with a minimum of 3). 
C------------
C   start MSTU            
C-MSTU(112)=4  !switches related to couplings
C-MSTU(113)=4
C-MSTU(114)=4
C           
C-------------------------------------------
C   start PARJ
       PARJ(1)=0.030  !(D=0.10) is ${\cal P}(\mathrm{q}\mathrm{q})/{\cal P}(\mathrm{q})$, the suppression of diquark-antidiquark pair production in the colour field, compared with quark-antiquark production
       PARJ(2)=0.15   !(D=0.30) is ${\cal P}(\mathrm{s})/{\cal P}(u)$, the suppression of $\mathrm{s}$ quark pair production in the field compared with $\u $ or $\d $ pair production. 
       PARJ(3)=0.25   !(D=0.4) is $({\cal P}(\u\mathrm{s})/{\cal P}(\u\d ))/({\cal P}(\mathrm{s})/
C {\cal P}(\d ))$, the extra suppression of strange diquark production compared with the 
C normal suppression of strange quarks.
       PARJ(4)=0.05   !(D=0.05) is $(1/3){\cal P}(\u\d _1)/{\cal P}(\u\d _0)$, 
C the suppression of spin 1 diquarks compared with spin 0 ones (excluding the factor 3 coming from spin counting).
       PARJ(11)=0.3   ! (D=0.5) is the probability that a light meson 
C (containing $\u $ and $\d $ quarks only) has spin 1 (with 1-PARJ(11) the probability for spin 0). 
       PARJ(14)=0.0   ! (D=0.) is the probability that a spin = 0 meson 
C is produced with an orbital angular momentum 1, for a total spin = 1. 
       PARJ(15)=0.0   ! (D=0.) is the probability that a spin = 1 meson 
C is produced with an orbital angular momentum 1, for a total spin = 0. 
       PARJ(16)=0.0   ! (D=0.) is the probability that a spin = 1 meson 
C is produced with an orbital angular momentum 1, for a total spin = 1. 
C-----------
       PARJ(23)=0.03  !(D=0.01, 2.) a fraction PARJ(23) of the Gaussian transverse momentum 
C distribution is taken to be a factor PARJ(24) larger than input in PARJ(21). 
C This gives a simple parameterization of non-Gaussian tails to the Gaussian shape assumed above. 
       PARJ(24)=2.5
       PARJ(32)=1.0  !    (D=1. GeV) is, with quark masses added, used to define the minimum 
C allowable energy of a colour-singlet parton system.
C-------------
       PARJ(41)= 1.1265625  !    (D=0.3, 0.58 GeV$^{-2}$) give the $a$ and $b$ 
C parameters of the symmetric Lund fragmentation function for MSTJ(11)=1, 4 and 5 (and MSTJ(11)=3 for ordinary hadrons).
       PARJ(42)= 0.3671875
       PARJ(45 )= 0.8   !    (D=0.5) the amount by which the effective $a$ parameter 
C in the Lund flavour dependent symmetric fragmentation function is assumed to be larger than 
C the normal $a$ when diquarks are produced. More specifically, referring to eq. ([*]), 
C C $a_{\alpha} = $PARJ(41) when considering the fragmentation of a quark and = PARJ(41) + PARJ(45) 
C for the fragmentation of a diquark, with corresponding expression for $a_{\beta}$ depending 
C on whether the newly created object is a quark or diquark (for an independent gluon jet 
C generated with MSTJ(2)= 2 or 4, replace PARJ(41) by PARJ(43)). In the popcorn model, 
C a meson created in between the baryon and antibaryon has $a_{\alpha} = a_{\beta} = $PARJ(41) + PARJ(45).
       return
       end
