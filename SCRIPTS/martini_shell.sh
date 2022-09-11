#!/bin/bash

source /home/ubuntu/DE_NOVO_PROTEINS/martini_3.0.b.3.2/PythonBollocks/bin/activate
/home/ubuntu/DE_NOVO_PROTEINS/martini_3.0.b.3.2/martinize -f $1 -o $2 -x $3 -ss /home/ubuntu/DE_NOVO_PROTEINS/SCRIPTS/helix.ssd -ff /home/ubuntu/DE_NOVO_PROTEINS/martini_3.0.b.3.2/martini303v.partition -elastic
cp Protein_A.itp $4
/home/ubuntu/DE_NOVO_PROTEINS/martini_3.0.b.3.2/bbsc.sh $4 $1
deactivate
