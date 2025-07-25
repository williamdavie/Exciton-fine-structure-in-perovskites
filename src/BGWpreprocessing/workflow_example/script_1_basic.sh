#!/bin/bash

# Runs all required quantum espresso files for calculation type: Basic -> fastest path to exciton binding energy. 

# This file runs 01-scf and 02-wfn, 03 needs 02 results. 

NP=32 
QEPATH=/lsc/opt/quantum-espresso-7.4.1/bin
BGWPATH=/local/data/public/wd324/BerkeleyGW-4.0/bin
PREFIX='CsPbI_0.0_0.0'
OMP_NUM_THREADS=1 # change as needed 

# ------------------------------------------------------------------------
# SCF
# ------------------------------------------------------------------------

cd ./01-scf
echo "Entered 01 SCF"

mpirun -np $NP $QEPATH/pw.x -inp ${PREFIX}_scf.in > ${PREFIX}_scf.out
wait
mpirun -np $NP $QEPATH/pw2bgw.x -inp pw2bgw.in > pw2bgw.out
wait

echo "Finished 01 SCF"

cd ..

echo "Exited 01 SCF"

# For basic script we only need to run WFN_co, WFNq_co, WFN_fi, WFNq_fi
# ------------------------------------------------------------------------
# WFN_co
# ------------------------------------------------------------------------

cp -r ./01-scf/${PREFIX}.save ./04-wfn_co/
cp -r ./01-scf/${PREFIX}.xml ./04-wfn_co/

cd ./04-wfn_co

echo "Entered 04 WFN_co"

mpirun -np $NP $QEPATH/pw.x -inp ${PREFIX}_wfn_co.in > ${PREFIX}_wfn_co.out
wait
mpirun -np $NP $QEPATH/pw2bgw.x -inp pw2bgw.in > pw2bgw.out
wait

cd .. 
echo "Exited 04 WFN_co"

# ------------------------------------------------------------------------
# WFNq_co
# ------------------------------------------------------------------------


cp -r ./01-scf/${PREFIX}.save ./05-wfnq_co/
cp -r ./01-scf/${PREFIX}.xml ./05-wfnq_co/

cd ./05-wfnq_co

echo "Entered 05 WFNq_co"

mpirun -np $NP $QEPATH/pw.x -inp ${PREFIX}_wfnq_co.in > ${PREFIX}_wfnq_co.out
wait
mpirun -np $NP $QEPATH/pw2bgw.x -inp pw2bgw.in > pw2bgw.out
wait

cd .. 
echo "Exited 05 WFNq_co"

# ------------------------------------------------------------------------
# WFN_fi
# ------------------------------------------------------------------------

cp -r ./01-scf/${PREFIX}.save ./06-wfn_fi/
cp -r ./01-scf/${PREFIX}.xml ./06-wfn_fi/

cd ./06-wfn_fi

echo "Entered 06 WFN_fi"

mpirun -np $NP $QEPATH/pw.x -inp ${PREFIX}_wfn_fi.in > ${PREFIX}_wfn_fi.out
wait
mpirun -np $NP $QEPATH/pw2bgw.x -inp pw2bgw.in > pw2bgw.out
wait

cd .. 
echo "Exited 06 WFN_fi"

# ------------------------------------------------------------------------
# WFNq_fi
# ------------------------------------------------------------------------

cp -r ./01-scf/${PREFIX}.save ./07-wfnq_fi/
cp -r ./01-scf/${PREFIX}.xml ./07-wfnq_fi/

cd ./07-wfnq_fi

echo "Entered 07 WFNq_fi"

mpirun -np $NP $QEPATH/pw.x -inp ${PREFIX}_wfnq_fi.in > ${PREFIX}_wfnq_fi.out
wait
mpirun -np $NP $QEPATH/pw2bgw.x -inp pw2bgw.in > pw2bgw.out
wait

cd .. 
echo "Exited 07 WFNq_fi"

echo "Job Done"
