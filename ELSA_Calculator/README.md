# ELSA Calculator
Python modules for designing highly non-repetitive Extra Long sgRNA Arrays (ELSAs) for multiplexed CRISPR applications.

When using this code, remember to cite its publication: Reis, Alexander C.*, Sean M. Halper*, Grace E. Vezeau, Daniel P. Cetnar, Ayaan Hossain, Phillip R. Clauer, and Howard M. Salis. "Simultaneous repression of multiple bacterial genes using nonrepetitive extra-long sgRNA arrays." Nature biotechnology 37, no. 11 (2019): 1294-1301.

*These authors contributed equally to this work.

https://www.nature.com/articles/s41587-019-0286-9

Correspondence should be addressed to H.M.S. (salis@psu.edu).

Abstract:
Microbial engineering often requires the regulation of many gene expression levels, for example, to redirect metabolic flows or create complex phenotypes. While CRISPR-based systems enable cross-species gene regulation, it remains difficult to express many single-guide RNAs within the same organism without triggering genetic instability, due to the presence of repetitive DNA sequences. Here we stably co-expressed up to 22 single-guide RNAs within highly compact extra-long sgRNA arrays (ELSAs) to simultaneously knock down the expression of many targeted genes by up to 3500-fold. To do this, we developed toolboxes of highly non-repetitive genetic parts, such as sgRNA handles, by combining biophysical modeling, biochemical characterization, and machine learning to identify sequence-function relationships. We then developed an automated algorithm to design ELSAs that regulate desired sets of genes, utilizing the developed toolboxes of non-repetitive genetic parts and 23 design rules quantifying DNA synthesis complexity, sgRNA expression levels, sgRNA target selection, and genetic stability. As demonstrations, we introduced designed ELSAs into the E. coli genome to create and stably maintain highly selective phenotypes, characterized by long-term growth, RT-qPCR, RNA-Seq, and LC-MS measurements. A 15-sgRNA ELSA targeting hisD, proC, lysA, tyrA, aroF, pheA, leuA, ilvD, and argH resulted in multiple amino acid auxotrophy. A 20-sgRNA ELSA targeting poxB, sdhC, sdhD, ackA, pta, and iclR increased succinic acid production by 150-fold. A 22-sgRNA ELSA targeting yncG, plsB, dkgA, yncE, ansP, narQ, yncH, adiA, iclR, ycfS, marR, and wzb disrupted stress response, inhibited membrane biosynthesis, and reduced persister cell formation by 21-fold after antibiotic treatment. Altogether, we show that ELSAs enable simultaneous regulation of many genes for diverse metabolic engineering and synthetic biology applications.
