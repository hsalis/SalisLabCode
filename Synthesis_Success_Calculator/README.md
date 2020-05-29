# Synthesis_Success_Calculator
Python modules for predicting the synthesis success rate (1st pass) for DNA fragments. Additional python modules for training Random Forest classifiers.

Relevant publication:

The Synthesis Success Calculator: Predicting the Rapid Synthesis of DNA Fragments with Machine Learning

Sean M. Halper, Ayaan Hossain, and Howard M. Salis

Correspondence should be addressed to H.M.S. (salis@psu.edu)

Abstract:

The synthesis and assembly of long DNA fragments has greatly accelerated synthetic biology and biotechnology research. However, long turnaround times or synthesis failures create unpredictable bottlenecks in the design-build-test-learn cycle. We developed a machine learning model, called the Synthesis Success Calculator, to predict whether a long DNA fragment can be readily synthesized with a short turnaround time. The model also identifies the sequence determinants associated with the synthesis outcome. We trained a random forest classifier using biophysical features and a compiled dataset of 1076 DNA fragment sequences to achieve high predictive performance (F1 score of 0.928 on 251 unseen sequences). Feature importance analysis revealed that repetitive DNA sequences were the most important contributor to synthesis failures. We then applied the Synthesis Success Calculator across large sequence datasets and found that 84.9% of the Escherichia coli MG1655 genome, but only 34.4% of sampled plasmids in NCBI, could be readily synthesized. Overall, the Synthesis Success Calculator can be applied on its own to prevent synthesis failures or embedded within optimization algorithms to design large genetic systems that can be rapidly synthesized and assembled.
