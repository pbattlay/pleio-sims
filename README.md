# pleio-sims
Eidos scripts for SLiM simulations from "Pleiotropy drives repeatability in the genetic basis of adaptation", Battlay, Yeaman &amp; Hodgins 2022
https://www.biorxiv.org/content/10.1101/2021.09.13.459985v1

Example run of veff-bg0.1-qtl5-pheno2.slim
```
# make directory for SLiM output
mkdir veff-bg0.1-qtl5-pheno2/

# run simulations in SLiM
for PLEIO in 0 0.25 0.5 0.75 0.99
do
for MIG in 0 0.001 0.005
do
for EFF in 0.1 0.5 1 5
do
for i in {1..100}
do
slim -d PLEIO=$PLEIO -d MIG=$MIG -d EFF=$EFF veff-bg0.1-qtl5-pheno2.slim \
| awk -v var="$i" '{print $0, var}' >> veff-bg0.1-qtl5-pheno2/QTL#5#pleio#$PLEIO#mig#$MIG#eff#$EFF#.txt
done
done
done
done

# make separate directory for formatted output files
mkdir veff-bg0.1-qtl5-pheno2-clean/

cd veff-bg0.1-qtl5-pheno2/

# format output files for analysis
for i in *
do
cat $i \
| grep -v "initialize" \
| grep -v "//" \
| grep ".* .* .* .*" > ../veff-bg0.1-qtl5-pheno2-clean/$i
done
```

Example repeatability analysis [here](results-analysis.R)
