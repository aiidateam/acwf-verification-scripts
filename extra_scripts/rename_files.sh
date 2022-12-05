#!/bin/bash

#plugin_names_start = [
#        'wien2k',
#        'vasp',
#        'siesta',
#        'quantum_espresso',
#        'gpaw',
#        'fleur',
#        'cp2k',
#        'castep',
#        'bigdft',
#        'abinit'
#        ]

echo "BE AWARE OF WHAT YOU ARE DOING!"
echo "This scripts simply takes the files from an input folder and copy them into the 'data' folder, changing their names."
echo "The input folder should be a folder coming from our internal google drive, that has a fixed structure of subfolders."
echo "Moreover this script assumes you already have correctly calculated the average of the all electron set, in 'average_ae_set' folder."


mkdir data

GDRIVE_FOLDER=$1
OXIDES_FOLD=$GDRIVE_FOLDER/outputs-oxides-verification-PBE-v1
UNARIES_FOLD=$GDRIVE_FOLDER/outputs-unaries-verification-PBE-v1

cp $OXIDES_FOLD/results-oxides-verification-PBE-v1-abinit-PseudoDojo-0.5b1-PBE-SR-standard-psp8.json data/results-oxides-verification-PBE-v1-abinit.json 
cp $OXIDES_FOLD/results-oxides-verification-PBE-v1-bigdft.json data/results-oxides-verification-PBE-v1-bigdft.json
cp $OXIDES_FOLD/results-oxides-verification-PBE-v1-cp2k_TZV2P.json data/results-oxides-verification-PBE-v1-cp2k.json
cp $OXIDES_FOLD/results-oxides-verification-PBE-v1-fleur.json data/results-oxides-verification-PBE-v1-fleur.json
cp $OXIDES_FOLD/results-oxides-verification-PBE-v1-gpaw.json data/results-oxides-verification-PBE-v1-gpaw.json
cp $OXIDES_FOLD/results-oxides-verification-PBE-v1-mk2-castep.json data/results-oxides-verification-PBE-v1-castep.json
cp $OXIDES_FOLD/results-oxides-verification-PBE-v1-quantum_espresso.json data/results-oxides-verification-PBE-v1-quantum_espresso.json
cp $OXIDES_FOLD/results-oxides-verification-PBE-v1-siesta.json data/results-oxides-verification-PBE-v1-siesta.json
cp $OXIDES_FOLD/results-oxides-verification-PBE-v1-vasp.json data/results-oxides-verification-PBE-v1-vasp.json
cp $OXIDES_FOLD/results-oxides-verification-PBE-v1-wien2k_dk_0.06.json data/results-oxides-verification-PBE-v1-wien2k.json

cp $UNARIES_FOLD/results-unaries-verification-PBE-v1-abinit-PseudoDojo-0.5b1-PBE-SR-standard-psp8.json data/results-unaries-verification-PBE-v1-abinit.json
cp $UNARIES_FOLD/results-unaries-verification-PBE-v1-bigdft.json data/results-unaries-verification-PBE-v1-bigdft.json
cp $UNARIES_FOLD/results-unaries-verification-PBE-v1-cp2k_TZV2P.json data/results-unaries-verification-PBE-v1-cp2k.json
cp $UNARIES_FOLD/results-unaries-verification-PBE-v1-fleur.json data/results-unaries-verification-PBE-v1-fleur.json
cp $UNARIES_FOLD/results-unaries-verification-PBE-v1-gpaw.json data/results-unaries-verification-PBE-v1-gpaw.json
cp $UNARIES_FOLD/results-unaries-verification-PBE-v1-mk2-castep.json data/results-unaries-verification-PBE-v1-castep.json
cp $UNARIES_FOLD/results-unaries-verification-PBE-v1-quantum_espresso.json data/results-unaries-verification-PBE-v1-quantum_espresso.json
cp $UNARIES_FOLD/results-unaries-verification-PBE-v1-siesta.json data/results-unaries-verification-PBE-v1-siesta.json
cp $UNARIES_FOLD/results-unaries-verification-PBE-v1-vasp.json data/results-unaries-verification-PBE-v1-vasp.json
cp $UNARIES_FOLD/results-unaries-verification-PBE-v1-wien2k-dk_0.06.json data/results-unaries-verification-PBE-v1-wien2k.json


cp average_ae_set/results-oxides-verification-PBE-v1-ae.json data/results-oxides-verification-PBE-v1-ae.json
cp average_ae_set/results-unaries-verification-PBE-v1-ae.json data/results-unaries-verification-PBE-v1-ae.json
