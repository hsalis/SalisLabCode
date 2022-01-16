# Synthesis_Success_Calculator
Python modules for predicting the synthesis success rate (1st pass) for DNA fragments. Additional python modules for training Random Forest classifiers.

When using this code, remember to cite its publication: Halper, Sean M., Ayaan Hossain, and Howard M. Salis. "Synthesis success calculator: predicting the rapid synthesis of DNA fragments with machine learning." ACS Synthetic Biology 9, no. 7 (2020): 1563-1571.
https://pubs.acs.org/doi/abs/10.1021/acssynbio.9b00460
Correspondence should be addressed to H.M.S. (salis@psu.edu)

Abstract:

The synthesis and assembly of long DNA fragments has greatly accelerated synthetic biology and biotechnology research. However, long turnaround times or synthesis failures create unpredictable bottlenecks in the design-build-test-learn cycle. We developed a machine learning model, called the Synthesis Success Calculator, to predict whether a long DNA fragment can be readily synthesized with a short turnaround time. The model also identifies the sequence determinants associated with the synthesis outcome. We trained a random forest classifier using biophysical features and a compiled dataset of 1076 DNA fragment sequences to achieve high predictive performance (F1 score of 0.928 on 251 unseen sequences). Feature importance analysis revealed that repetitive DNA sequences were the most important contributor to synthesis failures. We then applied the Synthesis Success Calculator across large sequence datasets and found that 84.9% of the Escherichia coli MG1655 genome, but only 34.4% of sampled plasmids in NCBI, could be readily synthesized. Overall, the Synthesis Success Calculator can be applied on its own to prevent synthesis failures or embedded within optimization algorithms to design large genetic systems that can be rapidly synthesized and assembled.
