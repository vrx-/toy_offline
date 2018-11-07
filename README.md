# toy_offline
test for ROMS offline biology

1) Compile and run upwelling hidrodynamics 

bin/build-coawst.py --clean upwelling

./coawstS < External/ocean_upwelling.in 

2) Prepare files for offline run 

python make_frc_from_his.py

python make_clm_from_his.py

python make_nud.py

python make_bio_csv.py 

python make_bio_ini.py (this one should be the last, the order of the others doesn't matter)

3) compile and run offline biology

bin/build-coawst.py --clean upwelling_offline

./coawstS < External/ocean_upwelling_offline.in
