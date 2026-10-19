#!/usr/bin/env bash
# 02_generate_beast_xml.sh
# This script does not generate the XML automatically because our model is highly custom. 
# It prints the exact BEAUti recipe for the hybrid model.

cat << 'EOF'
=== BEAST v10.5.0 XML GENERATION RECIPE (BEAUti) ===

1. TAXA & DATES:
   - Import the TempEst rerooted tree or alignment.
   - Import sampling dates (decimal years).

2. SITE MODEL:
   - Substitution Model: HKY (Empirical base frequencies).

3. CLOCK MODEL (The Hybrid Approach):
   - Select "Fixed Local Clock".
   - We will use a dummy rate (1.0) in BEAUti and script the exact priors later.
   - Alternatively, set the global clock rate to 1x10^-3.

4. TREE PRIOR (The Two-Phase Coalescent):
   - Select "Coalescent: Bayesian SkyGrid".
   - Set the number of parameters appropriately for your dataset.

5. MCMC:
   - Chain Length: 100,000,000
   - Log Every: 10,000
   - Echo state to screen: every 10,000

6. SAVE XML:
   - Save as: marburg_hybrid_model.xml

7. POST-BEAUti XML MANUAL EDITS (Crucial for your hybrid method):
   - Open marburg_hybrid_model.xml in a text editor (nano/VS Code).
   - Locate the <branchRateModel> tag. 
     -> Ensure the outbreak clade rate is fixed to 1e-3.
     -> Ensure the relaxed clock prior is LogNormal (M=5e-4, S=0.5).
   - Locate the <treeModel> tag.
     -> Ensure the prior is a mixture of ExponentialGrowth (for the Ethiopian lineage) and ConstantPopulation (for the deeper nodes).
   - Save the file.
EOF
