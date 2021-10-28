##Automated NMR Chemical Shift Perturbation Assignment


Example of simple peak assignement prediction
```
python estimateShifts.py --free test_in_free_assigned.csv  --with_ligand  test_in_ligand_unassigned.csv  --N_divisor 5 --alg RA          
```



Peak assignement prediction with known target peaks
```
python estimateShifts.py --free test_in_free_assigned.csv  --with_ligand  test_in_ligand_assigned.csv  --N_divisor 5 --alg RA --labeled            
```

Run examples using Docker:
```
 docker  run  -v `pwd`:`pwd` -w `pwd` -i -t cermcirmmp/picasso python /picasso/estimateShifts.py --free /picasso/test_in_free_assigned.csv  --with_ligand  /picasso/test_in_ligand_unassigned.csv  --N_divisor 5 --alg RA
```
```
 docker  run  -v `pwd`:`pwd` -w `pwd` -i -t cermcirmmp/picasso python /picasso/estimateShifts.py --free /picasso/test_in_free_assigned.csv  --with_ligand  /picasso/test_in_ligand_assigned.csv  --N_divisor 5 --alg RA --labeled
 
 ```
