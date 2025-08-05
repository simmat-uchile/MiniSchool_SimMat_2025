%oldchk=benceno_FC.chk
%chk=benceno_FC_NO.chk
# Guess=(read,Save,Only,NaturalOrbitals) Geom=AllCheck ChkBasis

--link1--
%OldChk=benceno_FC.chk
%Chk=benceno_state1.chk
# Geom=AllCheck Guess=(Read,Only) Density=(Check,Transition=1) Pop=(Minimal,NTO,SaveNTO)

