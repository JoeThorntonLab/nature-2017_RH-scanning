# nature-2017_RH-scanning

Input data and scripts to recreate all figures and generate underlying data for Starr, Picton, Thornton "Alternative evolutionary histories in the sequence space of an ancient protein", Nature 2017.

For each script, create an output subdirectory "x_script-name_out/"

Input data gives the read counts mapped to each RH combination coming from each library/barcode combination.

Run scripts in order by file name. All scripts are preceded by a pointer to the source file location, which needs to be changed accordingly. Be sure all packages referenced in the beginning of each script are downloaded:

1_calc-meanF_11P.R : calculates the ML mean fluorescence for each RH variant on ERE and SRE in the AncSR1+11P background

2_calc-meanF_SR1.R : calculates the ML mean fluorescence for each RH variant on ERE and SRE in the AncSR1 background

3_classify-variant-distributions.R : generates null sampling distributions for classifying RH variants. (takes a while to generate positive sampling distributions)

4_classify-variants.R : classifies each variant in each background as null, weak, or strong on each RE

5_assess-replicates.R : generates some quality summary statistics; applies minimum # cells sampling cutoff to data based on these statistics

6_predict-missing.R : fits continuation ratio ordinal logistic regression models to the data to predict missing genotypes. Cross-validation is used to select the lambda penalization parameter. (Takes a long time!)

7_assess-connectivity.R : prepares force-directed graph files for analysis in gephi, and computes characteristics of trajectories available in each DBD background using igraph

8_bchem-determ.R : prepare alignments for logo generation; infer logistic regression models to reveal biophysical determinants of RE-specificity; analyze foldX data

9_quant-scale-landscape : look into error structure of quantitative mean fluorescence estimates; consider alternative evolutionary regime where SRE mean fluorescence must increase with each step 

10_alt-class.R : test robustness of connectivity results to alternate schemes for classifying "specific" and functional variants
