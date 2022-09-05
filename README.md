# pleio-sims
Eidos scripts for SLiM simulations from "Pleiotropy drives repeatability in the genetic basis of adaptation" Battlay, Yeaman &amp; Hodgins 2022
https://www.biorxiv.org/content/10.1101/2021.09.13.459985v1

# QTL5-EFF
# run the simulation 100x for each combination of parameters
for PLEIO in 0 0.25 0.5 0.75 0.99
do
for MIG in 0 0.001 0.005
do
for EFF in 0.1 0.5 1 5
do
for i in {1..100}
do
slim -d PLEIO=$PLEIO -d MIG=$MIG -d EFF=$EFF QTL5-EFF.slim \
| awk -v var="$i" '{print $0, var}' >> QTL5-EFF/QTL#5#pleio#$PLEIO#mig#$MIG#eff#$EFF#.txt
done
done
done
done
