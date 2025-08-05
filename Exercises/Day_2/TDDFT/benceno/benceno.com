%nprocshared=4
%chk=benceno_GS.chk
%mem=2GB
#P B3LYP/Def2TZVP opt freq 

benceno_GS

0 1
  6                  -0.000000    1.391233    0.000000
  6                  -1.204843    0.695616    0.000000
  6                   1.204843    0.695616    0.000000
  1                  -2.142801    1.237147    0.000000
  1                   2.142801    1.237147    0.000000
  6                  -1.204843   -0.695616    0.000000
  6                   1.204843   -0.695616    0.000000
  1                  -2.142801   -1.237147    0.000000
  1                   2.142801   -1.237147    0.000000
  6                   0.000000   -1.391233    0.000000
  1                   0.000000   -2.474293    0.000000
  1                  -0.000000    2.474293    0.000000

--Link1--
%oldchk=benceno_GS.chk
%chk=benceno_FC.chk
#P geom=check guess=read B3LYP/Def2TZVP td(singlet,root=1,NStates=20) density=current 

benceno_FC

0 1

