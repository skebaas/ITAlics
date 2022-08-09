#!/usr/bin/python3

# Script checks outputs from PhysioCorrection.sh for excessive 0s
# If a regressor column has non-zero data in less than 75% of its cells,
# it will be deleted

import sys
import pandas as pd

file = str(sys.argv[1])

df = pd.read_csv(file, header=None, index_col=None, sep="\t")
df = df.drop(columns=24)
bad = 0

for column in df:
	count = 0
	for value in df[column]:
		if float(value) == 0.0:
			count = count+1
		else:
			continue
	if (count / len(df[column])) > .25:
		print("Too many zeros found in regressor column " + str(column) + " for " + file + ", deleting column")
		df = df.drop(columns=column)
		with open('physio_qc_output.txt', 'a') as out:
			out.write("Too many zeros found in regressor column " + str(column) + " for " + file + ", deleting column\n\n")
		bad = bad + 1
	else:
		continue

if bad > 0:
	df.to_csv(file, sep="\t")
else:
	print("Regressors in " + file + " are ok")