In this folder there is a comparison between the kpoints mesh
generated with distance 0.06 1/ang and with distance 0.045 1/ang
performed with the WIEN2K code. 

The sets:
results-oxides-0045.json  has original name results-wien2k-finite-T_delta-k_0.045_prec3.json
results-oxides-006.json   has original name results-wien2k-finite-T_delta-k_0.06_prec3.json
results-unaries-0045.json was created in two stages by WIEN2K team. File here the only available copy.
results-unaries-006.json  has original name results-unaries-verification-PBE-v1-wien2k-dk_0.06.json 

IMPORTANT NOTE: the oxides results are not the final results used to create the referencee AE
set, they are an iteration before that. This means the central volumes are not the same as the
final reference, but also parameters for few materials have changed.

The script also prints the worst cases. Worst cases means % relavive difference more than
0.15 % for "V0_rel_diff", more than 1 % for "B0_rel_diff and more than 5%
for "B1_rel_diff". Here the result:

unaries V0_rel_diff

oxides V0_rel_diff
He-XO -0.16452442159383504
Rb-XO3 3.7187973022597616

unaries B0_rel_diff
C-X/SC 1.0694927226764814
Hg-X/SC 1.5698113400457128
Pa-X/SC -1.6633317272992867
Rn-X/Diamond -1.9233115968839603

oxides B0_rel_diff
Ba-XO3 -1.0488452980189695
Be-XO 1.0808517945200096
He-XO 3.332446707451444
Ra-XO3 1.1115049419760876
Rb-XO3 -142.29005661636492

unaries B1_rel_diff
Ce-X/SC -8.656625878593145
Er-X/Diamond 6.81933569895974
Fr-X/SC -5.558639915382058
Hg-X/BCC 8.159067558522283
Hg-X/FCC -14.737736507385087
N-X/SC -12.025286772438038
Rn-X/BCC 5.664839549341651
Rn-X/Diamond -23.52180516939958
Rn-X/FCC -14.464142350519863
Xe-X/Diamond -16.220704351427464
Xe-X/SC -8.07486322815154

oxides B1_rel_diff
Be-XO -15.322751851343087
F-XO3 7.268522342446507
He-XO -7.898116438356166
La-XO3 -7.300487220239187
Rb-XO3 149.34950503336557
Th-XO2 20.32113304453841

