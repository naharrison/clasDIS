# clasDIS
Pythia based Deep Inelastic Scattering event generator for CLAS

Compiles but doesn't run.
Most likely a problem with clas-trunk;
Probably need to compile it, not just download a copy (see https://github.com/naharrison/docker-clas6).

## Requirements
python2-scons.noarch <br>
gcc-gfortran.x86_64 <br>
libX11-devel.x86_64 <br>
cernlib (https://github.com/naharrison/cernlib-repo) <br>
clas-trunk (https://github.com/naharrison/clas-trunk) <br>
<br>
setenv PYTHONPATH /.../clas-trunk/ <br>
setenv CERN /.../cernlib-repo/cernlib/x86_64_rhel7/2005/

## Example
`clasdis --trig 20000 --datf --outform 2 --beam 5498 --zpos -250 --zwidth 25 --t 5 60 --parl3 0.7 --lst15 145 --path ./CDOUT --parj33 0.3`
