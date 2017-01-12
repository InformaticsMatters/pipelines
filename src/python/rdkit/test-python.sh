#!/usr/bin/env bash

rm -rf ../../../tmp
mkdir ../../../tmp
msg_fail="\n========================== TEST FAILED ==========================\n"
msg_file_notCreated="\n========================== FILE NOT CREATED =====================\n"


echo "Testing screen.py reading from STDIN and writing to STDOUT input as smiles"
zcat ../../../data/dhfr_3d.sdf.gz | python screen.py\
  --smiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2'\
  --simmin 0.45\
  -if sdf > /dev/null || echo -e $msg_fail

echo "Testing screen.py reading from STDIN and writing to STDOUT input from molfile"
zcat ../../../data/dhfr_3d.sdf.gz | python screen.py\
  --molfile ../../../data/pyrimethamine.mol\
  --simmin 0.7\
  --simmax 0.8\
  -if sdf > /dev/null || echo -e $msg_fail

echo "Testing screen.py reading and writing files using sdf"
python screen.py\
  -i ../../../data/dhfr_3d.sdf.gz\
  -o ../../../tmp/screen2\
  -of sdf\
  --smiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2'\
  --simmin 0.45\
  -if sdf > /dev/null || echo -e $msg_fail
if [ ! -f  ../../../tmp/screen2.sdf.gz ]
then
    echo -e $msg_file_notCreated
fi

echo "Testing screen.py reading and writing files using json"
python screen.py\
  -i ../../../data/nci100.data.gz\
  -o ../../../tmp/screen3\
  -of json\
  --smiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2'\
  --simmin 0.45\
  -if json > /dev/null || echo -e $msg_fail
if [ ! -f  ../../../tmp/screen3.data.gz ]
then
    echo -e $msg_file_notCreated
fi


echo "Testing conformers.py reading from STDIN and writing to STDOUT"
zcat ../../../data/Kinase_inhibs.sdf.gz | python conformers.py -n 2 -if sdf > /dev/null || echo -e $msg_fail

echo "Testing o3dAlign.py reading from STDIN and writing to STDOUT"
zcat ../../../data/Kinase_inhibs.sdf.gz | python o3dAlign.py\
  ../../../data/pyrimethamine.mol -n 2 -t 10 -if sdf > /dev/null || echo -e $msg_fail

echo "Testing filter.py reading from STDIN and writing to STDOUT"
zcat ../../../data/Kinase_inhibs.sdf.gz | python filter.py --hacmin 25 --hacmax 30 -if sdf > /dev/null || echo -e $msg_fail

echo "Testing cluster_butina.py reading from STDIN and writing to STDOUT"
zcat ../../../data/Kinase_inhibs.sdf.gz | python cluster_butina.py -t 0.6 -if sdf > /dev/null || echo -e $msg_fail

echo "Finished"
