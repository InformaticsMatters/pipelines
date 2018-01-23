#!/usr/bin/env bash

rm -rf tmp
mkdir tmp
msg_fail="\n========================== TEST FAILED ==========================\n"
msg_file_notCreated="\n========================== FILE NOT CREATED =====================\n"


echo "Testing screen.py reading from STDIN and writing to STDOUT input as smiles"
gunzip -c data/dhfr_3d.sdf.gz | python src/python/pipelines/rdkit/screen.py\
  --qsmiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2'\
  --simmin 0.45\
  -if sdf > /dev/null || echo -e $msg_fail


echo "Testing screen.py reading from STDIN and writing to STDOUT input from molfile"
gunzip -c data/dhfr_3d.sdf.gz | python src/python/pipelines/rdkit/screen.py\
  --qmolfile data/pyrimethamine.mol\
  --simmin 0.7\
  --simmax 0.8\
  -if sdf > /dev/null || echo -e $msg_fail

#echo "Testing screen.py reading and writing files using sdf"
#python src/python/pipelines/rdkit/screen.py\
#  -i data/dhfr_3d.sdf.gz\
#  -o tmp/screen2\
#  -of sdf\
#  --qsmiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2'\
#  --simmin 0.45\
#  -if sdf > /dev/null || echo -e $msg_fail
#if [ ! -f  tmp/screen2.sdf.gz ]
#then
#    echo -e $msg_file_notCreated
#fi

#echo "Testing screen.py reading and writing files using json"
#python src/python/pipelines/rdkit/screen.py\
#  -i data/nci100.data.gz\
#  -o tmp/screen3\
#  -of json\
#  --qsmiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2'\
#  --simmin 0.45\
#  -if json > /dev/null || echo -e $msg_fail
#if [ ! -f  tmp/screen3.data.gz ]
#then
#    echo -e $msg_file_notCreated
#fi

#echo "Testing screen.py reading and writing files using thin sdf"
#python src/python/pipelines/rdkit/screen.py\
#  -i data/dhfr_3d.sdf.gz\
#  -o tmp/screen2\
#  -of sdf\
#  --qsmiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2'\
#  --simmin 0.45\
#  -if sdf\
#  --thin > /dev/null || echo -e $msg_fail
#if [ ! -f  tmp/screen2.sdf.gz ]
#then
#    echo -e $msg_file_notCreated
#fi

#echo "Testing screen.py reading and writing files using thin json"
#python src/python/pipelines/rdkit/screen.py\
#  -i data/nci100.data.gz\
#  -o tmp/screen3\
#  -of json\
#  --qsmiles 'C1N=C(C2=CC=CC=C2)C2=CC=CC=C2C2=C1C=NC(NC1=CC=CC=C1)=N2'\
#  --simmin 0.45\
#  -if json\
#  --thin > /dev/null || echo -e $msg_fail
#if [ ! -f  tmp/screen3.data.gz ]
#then
#    echo -e $msg_file_notCreated
#fi

echo "Testing screen_multi.py reading taget form sdf file, query as json file and writing to STDOUT"
gunzip -c data/dhfr_3d.sdf.gz | python src/python/pipelines/rdkit/screen_multi.py\
  -if sdf --qjson data/nci100.data.gz --simmin 0.55 > /dev/null || echo -e $msg_fail

echo "Testing conformers.py reading from STDIN and writing to STDOUT"
gunzip -c data/Kinase_inhibs.sdf.gz | python src/python/pipelines/rdkit/conformers.py -n 2 -if sdf > /dev/null || echo -e $msg_fail

echo "Testing constrained_conf_gen.py reading from STDIN and writing to STDOUT"
python src/python/pipelines/rdkit/constrained_conf_gen.py -n 2 -i data/XChemReactionMaker1.sdf.gz -r data/ref_mol.sdf.gz -o tmp/constrained_conf_gen.sdf.gz || echo -e $msg_fail

echo "Testing conformers.py with clustering reading from STDIN and writing to STDOUT"
gunzip -c data/Kinase_inhibs.sdf.gz | python src/python/pipelines/rdkit/conformers.py -n 2 -c RMSD -if sdf > /dev/null || echo -e $msg_fail

echo "Testing o3dAlign.py reading from STDIN and writing to STDOUT"
gunzip -c data/Kinase_inhibs.sdf.gz | python src/python/pipelines/rdkit/o3dAlign.py\
  data/pyrimethamine.mol -n 2 -t 10 -if sdf > /dev/null || echo -e $msg_fail

echo "Testing filter.py reading from STDIN and writing to STDOUT"
gunzip -c data/Kinase_inhibs.sdf.gz | python src/python/pipelines/rdkit/filter.py --hacmin 25 --hacmax 30 -if sdf > /dev/null || echo -e $msg_fail

echo "Testing cluster_butina.py reading from STDIN and writing to STDOUT"
gunzip -c data/Kinase_inhibs.sdf.gz | python src/python/pipelines/rdkit/cluster_butina.py -t 0.6 -if sdf > /dev/null || echo -e $msg_fail

echo "Testing cluster_butina_matrix.py reading from STDIN and writing TSV to STDOUT"
gunzip -c data/Kinase_inhibs.sdf.gz | src/python/pipelines/rdkit/cluster_butina_matrix.py -t 0.6 -if sdf -of tsv > /dev/null || echo -e $msg_fail

echo "Testing cluster_butina_matrix.py reading from STDIN and writing JSON to STDOUT"
gunzip -c data/Kinase_inhibs.sdf.gz | src/python/pipelines/rdkit/cluster_butina_matrix.py -t 0.6 -if sdf -of json > /dev/null || echo -e $msg_fail

echo "Testing rxn_smarts_filter.py reading from STDIN and writing to files using SDF"
gunzip -c data/Kinase_inhibs.sdf.gz | python src/python/pipelines/rdkit/rxn_smarts_filter.py -if sdf -o tmp/rxn_smarts_filter1 || echo -e $msg_fail

echo "Testing rxn_smarts_filter.py reading from STDIN and writing to files using JSON"
gunzip -c data/nci100.data.gz | python src/python/pipelines/rdkit/rxn_smarts_filter.py -if json -o tmp/rxn_smarts_filter2 -of json --meta --thin || echo -e $msg_fail


#echo "Testing rxn_smarts_filter.py reading from sd file and writing to multiple files"
#python src/python/pipelines/rdkit/rxn_smarts_filter.py -i data/Kinase_inhibs.sdf.gz -o tmp/rxn_smarts_filter2 --multi || echo -e $msg_fail

#echo "Testing rxn_maker.py reading from files"
#python src/python/pipelines/rdkit/rxn_maker.py -i data/sulfonyl_chloride.sdf -r Sulfonamide -rl data/sdf-aliphatic-primary-amines-175.sdf.gz  -o tmp/rxnoutput || echo -e $msg_fail

#echo "Testing pbf_ev.py reading from files"
#python src/python/pipelines/rdkit/pbf_ev.py -i data/dhfr_3d.sdf -o tmp/pbf_ev_output || echo -e $msg_fail

echo "cleaning up ..."
rm -f *.metadata

echo "Finished"
