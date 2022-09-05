#####
# repeatablity analysis for a set of SLiM simulations
#####

library(tidyr)
library(dgconstraint)

# set directory to 'cleaned' output
setwd("veff-bg0.1-qtl5-pheno2-clean/")

# list output files
file.names = dir(pattern = "QTL#")

# convert file names to table
df = separate(as.data.frame(file.names),
	col = file.names,
	sep = "#", convert = T,
	remove = T,
	into = c("QTL0", "QTL", "pleio0", "pleio", "mig0", "mig", "eff0", "eff", "ext"),
	extra = "warn", fill = "warn")[c(2,4,6,8)]

# set generation to analyze
GEN = 40000

# summary stats
df$pheno0.div = NA
df$pheno0.div.sem = NA
df$pheno1.div = NA
df$pheno1.div.sem = NA
df$pheno0.fnf.rat = NA
df$pheno1.fnf.rat = NA
df$pheno0.cchisq = NA
df$pheno1.cchisq = NA

# loop over each parameter combination
for(k in 1:length(file.names)){
	slim = read.table(file.names[k], sep = " ", header = F)

	# two phenotypes
	colnames(slim) = c("Generation", "Location", "pheno0.Effect", "pheno1.Effect", "p1.Freq", "p2.Freq", "SimRun")

	# five phenotypes
	#colnames(slim) = c("Generation", "Location", "pheno0.Effect", "pheno1.Effect", "pheno2.Effect", "pheno3.Effect", "pheno4.Effect", "p1.Freq", "p2.Freq", "SimRun")

	# ten phenotypes
	#colnames(slim) = c("Generation", "Location", "pheno0.Effect", "pheno1.Effect", "pheno2.Effect", "pheno3.Effect", "pheno4.Effect", "pheno5.Effect", "pheno6.Effect", "pheno7.Effect", "pheno8.Effect", "pheno9.Effect", "p1.Freq", "p2.Freq", "SimRun")	

	slim = subset(slim, Generation == GEN)

	# calculate phenotype 0 genetic value for each mutation
	slim$mut.genetic.value.pheno0 = (slim$p1.Freq - slim$p2.Freq) * slim$pheno0.Effect

	# calculate phenotype 1 genetic value for each mutation
	slim$mut.genetic.value.pheno1 = (slim$p1.Freq - slim$p2.Freq) * slim$pheno1.Effect

	# define QTL numbers (five QTL)
	slim$QTL = NA
	slim[which(slim$Location >= 0 & slim$Location < 500), ]$QTL = 1
	slim[which(slim$Location >= 500 & slim$Location < 1000), ]$QTL = 2
	slim[which(slim$Location >= 1000 & slim$Location < 1500), ]$QTL = 3
	slim[which(slim$Location >= 1500 & slim$Location < 2000), ]$QTL = 4
	slim[which(slim$Location >= 2000 & slim$Location < 2500), ]$QTL = 5

	# define QTL numbers (20 QTL)
	#slim$QTL = NA
	#slim[which(slim$Location >= 0 & slim$Location < 500), ]$QTL = 1
	#slim[which(slim$Location >= 500 & slim$Location < 1000), ]$QTL = 2
	#slim[which(slim$Location >= 1000 & slim$Location < 1500), ]$QTL = 3
	#slim[which(slim$Location >= 1500 & slim$Location < 2000), ]$QTL = 4
	#slim[which(slim$Location >= 2000 & slim$Location < 2500), ]$QTL = 5
	#slim[which(slim$Location >= 2500 & slim$Location < 3000), ]$QTL = 6
	#slim[which(slim$Location >= 3000 & slim$Location < 3500), ]$QTL = 7
	#slim[which(slim$Location >= 3500 & slim$Location < 4000), ]$QTL = 8
	#slim[which(slim$Location >= 4000 & slim$Location < 4500), ]$QTL = 9
	#slim[which(slim$Location >= 4500 & slim$Location < 5000), ]$QTL = 10
	#slim[which(slim$Location >= 5000 & slim$Location < 5500), ]$QTL = 11
	#slim[which(slim$Location >= 5500 & slim$Location < 6000), ]$QTL = 12
	#slim[which(slim$Location >= 6000 & slim$Location < 6500), ]$QTL = 13
	#slim[which(slim$Location >= 6500 & slim$Location < 7000), ]$QTL = 14
	#slim[which(slim$Location >= 7000 & slim$Location < 7500), ]$QTL = 15
	#slim[which(slim$Location >= 7500 & slim$Location < 8000), ]$QTL = 16
	#slim[which(slim$Location >= 8000 & slim$Location < 8500), ]$QTL = 17
	#slim[which(slim$Location >= 8500 & slim$Location < 9000), ]$QTL = 18
	#slim[which(slim$Location >= 9000 & slim$Location < 9500), ]$QTL = 19
	#slim[which(slim$Location >= 9500 & slim$Location < 10000), ]$QTL = 20

	# summarize by simulation rep
	slim.by.sim = slim %>%
		group_by(SimRun) %>%
		summarise(sim.genetic.value.pheno0 = 2 * sum(mut.genetic.value.pheno0),
			sim.genetic.value.pheno1 = 2 * sum(mut.genetic.value.pheno1))

	# output to main data frame
	df$pheno0.div[k] = mean(slim.by.sim$sim.genetic.value.pheno0)
	df$pheno0.div.sem[k] = sd(slim.by.sim$sim.genetic.value.pheno0) / sqrt(dim(slim.by.sim)[1])
	df$pheno1.div[k] = mean(slim.by.sim$sim.genetic.value.pheno1)
	df$pheno1.div.sem[k] = sd(slim.by.sim$sim.genetic.value.pheno1) / sqrt(dim(slim.by.sim)[1])

	# individual QTL stats
	SimRuns = unique(slim$SimRun)
	QTLs = unique(slim$QTL)
	
	# PHENOTYPE 0
	slim.by.simQTL = matrix(ncol = length(SimRuns), nrow = length(QTLs) + 1)

	# loop through each QTL
	for (j in QTLs) {

		# loop through each SLiM run
		for (i in SimRuns) {
			# get genetic value for each QTL in each rep
			slim.by.simQTL[j,i] = 2 * sum(slim[which(slim$SimRun == i & slim$QTL == j), ]$mut.genetic.value.pheno0)
		}
	}

	# get genetic value for all QTL
	for (i in SimRuns) {
		slim.by.simQTL[length(QTLs) + 1, i] = 2 * sum(slim[which(slim$SimRun == i), ]$mut.genetic.value.pheno0)
	}

	df$pheno0.fnf.rat[k] = mean(slim.by.simQTL[1, ] / slim.by.simQTL[6, ])
	df$pheno0.cchisq[k] = pairwise_c_chisq(slim.by.simQTL[1:length(QTLs), ], num_permute = 1000, na.rm = F)

	# PHENOTYPE 1
	slim.by.simQTL = matrix(ncol = length(SimRuns), nrow = length(QTLs) + 1)

	# loop through each QTL
	for (j in QTLs) {

		# loop through each SLiM run
		for (i in SimRuns) {
			# get genetic value for each QTL in each rep
			slim.by.simQTL[j,i] = 2 * sum(slim[which(slim$SimRun == i & slim$QTL == j), ]$mut.genetic.value.pheno1)
		}
	}

	# get genetic value for all QTL
	for (i in SimRuns) {
		slim.by.simQTL[length(QTLs) + 1, i] = 2 * sum(slim[which(slim$SimRun == i), ]$mut.genetic.value.pheno1)
	}

	df$pheno1.fnf.rat[k] = mean(slim.by.simQTL[1, ] / slim.by.simQTL[6, ])
	df$pheno1.cchisq[k] = pairwise_c_chisq(slim.by.simQTL[1:length(QTLs), ], num_permute = 1000, na.rm = F)

}

# save output
write.table(df, file = paste("veff-bg0.1-qtl5-pheno2-", GEN, ".df", sep = ""), quote = F, row.names = F)

