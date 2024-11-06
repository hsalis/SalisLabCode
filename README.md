**Salis Lab Source Code and Relevant Publications**

_If you use our code, please cite our work!_

This source code is available under MIT or GPLv3 licenses. See LICENSE files for more information.

**Promoter Calculator v1.0**

Transcription rates are regulated by the interactions between RNA polymerase, sigma factor, and promoter DNA sequences in bacteria. However, it remains unclear how non-canonical sequence motifs collectively control transcription rates. Here, we combine massively parallel assays, biophysics, and machine learning to develop a 346-parameter model that predicts site-specific transcription initiation rates for any σ<sup>70</sup> promoter sequence, validated across 22132 bacterial promoters with diverse sequences. We apply the model to predict genetic context effects, design σ<sup>70</sup> promoters with desired transcription rates, and identify undesired promoters inside engineered genetic systems. The model provides a biophysical basis for understanding gene regulation in natural genetic systems and precise transcriptional control for engineering synthetic genetic systems.

Relevant Article: LaFleur, Travis L., Ayaan Hossain, and Howard M. Salis. "Automated model-predictive design of synthetic promoters to control transcriptional profiles in bacteria." Nature communications 13, no. 1 (2022): 5159.

Link to Article: https://www.nature.com/articles/s41467-022-32829-5

**mRNA Stability Calculator**

mRNA degradation is a central process that affects all gene expression levels, though it remains challenging to predict the stability of a mRNA from its sequence, due to the many coupled interactions that control degradation rate. Here, we carried out massively parallel kinetic decay measurements on 62120 bacterial mRNAs, using a learn-by-design approach to develop and validate a predictive sequence-to-function model of mRNA stability. mRNAs were designed to systematically vary translation rates, secondary structures, sequence compositions, G-quadruplexes, i-motifs, and RppH activity, resulting in mRNA half-lives from about 20 seconds to 20 minutes. We combined biophysical models and machine learning to develop steady-state and kinetic decay models of mRNA stability with high accuracy and generalizability, utilizing transcription rate models to identify mRNA isoforms and translation rate models to calculate ribosome protection. Overall, the developed model quantifies the key interactions that collectively control mRNA stability in bacterial operons and predicts how changing mRNA sequence alters mRNA stability, which is important when studying and engineering bacterial genetic systems.

Relevant Article: Cetnar, Daniel P., Ayaan Hossain, Grace E. Vezeau, and Howard M. Salis. "Predicting synthetic mRNA stability using massively parallel kinetic measurements, biophysical modeling, and machine learning." Nature communications 15, no. 1 (2024): 9601.

Link to Article: https://www.nature.com/articles/s41467-024-54059-7

**Gene Expression Model Test System**

Gene expression models greatly accelerate the engineering of synthetic metabolic pathways and genetic circuits by predicting sequence-function relationships and reducing trial-and-error experimentation. However, developing models with more accurate predictions remains a significant challenge. Here we present a model test system that combines advanced statistics, machine learning, and a database of 9862 characterized genetic systems to automatically quantify model accuracies, accept or reject mechanistic hypotheses, and identify areas for model improvement. We also introduce model capacity, a new information theoretic metric for correct cross-data-set comparisons. We demonstrate the model test system by comparing six models of translation initiation rate, evaluating 100 mechanistic hypotheses, and uncovering new sequence determinants that control protein expression levels. We then applied these results to develop a biophysical model of translation initiation rate with significant improvements in accuracy. Automated model test systems will dramatically accelerate the development of gene expression models, and thereby transition synthetic biology into a mature engineering discipline.

Relevant Article: Reis, Alexander C., and Howard M. Salis. "An automated model test system for systematic development and improvement of gene expression models." ACS synthetic biology 9, no. 11 (2020): 3145-3156.

Link to Article: https://pubs.acs.org/doi/abs/10.1021/acssynbio.0c00394

**Non-Repetitive Parts Calculator**

Engineered genetic systems are prone to failure when their genetic parts contain repetitive sequences. Designing many non-repetitive genetic parts with desired functionalities remains a difficult challenge with high computational complexity. To overcome this challenge, we developed the Nonrepetitive Parts Calculator to rapidly generate thousands of highly non-repetitive genetic parts from specified design constraints, including promoters, ribosome-binding sites and terminators. As a demonstration, we designed and experimentally characterized 4,350 non-repetitive bacterial promoters with transcription rates that varied across a 820,000-fold range, and 1,722 highly nonrepetitive yeast promoters with transcription rates that varied across a 25,000-fold range. We applied machine learning to explain how specific interactions controlled the promoters’ transcription rates. We also show that using non-repetitive genetic parts substantially reduces homologous recombination, resulting in greater genetic stability.

Relevant Article: Hossain, Ayaan, Eriberto Lopez, Sean M. Halper, Daniel P. Cetnar, Alexander C. Reis, Devin Strickland, Eric Klavins, and Howard M. Salis. "Automated design of thousands of nonrepetitive parts for engineering stable genetic systems." Nature biotechnology 38, no. 12 (2020): 1466-1475.

Link to Article: https://www.nature.com/articles/s41587-020-0584-2

**Synthesis Success Calculator v1.0**

The synthesis and assembly of long DNA fragments has greatly accelerated synthetic biology and biotechnology research. However, long turnaround times or synthesis failures create unpredictable bottlenecks in the design–build–test–learn cycle. We developed a machine learning model, called the Synthesis Success Calculator, to predict whether a long DNA fragment can be readily synthesized with a short turnaround time. The model also identifies the sequence determinants associated with the synthesis outcome. We trained a random forest classifier using biophysical features and a compiled data set of 1076 DNA fragment sequences to achieve high predictive performance (F1 score of 0.928 on 251 unseen sequences). Feature importance analysis revealed that repetitive DNA sequences were the most important contributor to synthesis failures. We then applied the Synthesis Success Calculator across large sequence data sets and found that 84.9% of the Escherichia coli MG1655 genome, but only 34.4% of sampled plasmids in NCBI, could be readily synthesized. Overall, the Synthesis Success Calculator can be applied on its own to prevent synthesis failures or embedded within optimization algorithms to design large genetic systems that can be rapidly synthesized and assembled.

Relevant Article: Halper, Sean M., Ayaan Hossain, and Howard M. Salis. "Synthesis success calculator: predicting the rapid synthesis of DNA fragments with machine learning." ACS Synthetic Biology 9, no. 7 (2020): 1563-1571.

Link to Article: https://pubs.acs.org/doi/abs/10.1021/acssynbio.9b00460

**ELSA Calculator (Extra Long sgRNA Arrays)**

Engineering cellular phenotypes often requires the regulation of many genes. When using CRISPR interference, coexpressing many single-guide RNAs (sgRNAs) triggers genetic instability and phenotype loss, due to the presence of repetitive DNA sequences. We stably coexpressed 22 sgRNAs within nonrepetitive extra-long sgRNA arrays (ELSAs) to simultaneously repress up to 13 genes by up to 3,500-fold. We applied biophysical modeling, biochemical characterization and machine learning to develop toolboxes of nonrepetitive genetic parts, including 28 sgRNA handles that bind Cas9. We designed ELSAs by combining nonrepetitive genetic parts according to algorithmic rules quantifying DNA synthesis complexity, sgRNA expression, sgRNA targeting and genetic stability. Using ELSAs, we created three highly selective phenotypes in <i>Escherichia coli</i>, including redirecting metabolism to increase succinic acid production by 150-fold, knocking down amino acid biosynthesis to create a multi-auxotrophic strain and repressing stress responses to reduce persister cell formation by 21-fold. ELSAs enable simultaneous and stable regulation of many genes for metabolic engineering and synthetic biology applications.

Relevant Article: Reis, Alexander C., Sean M. Halper, Grace E. Vezeau, Daniel P. Cetnar, Ayaan Hossain, Phillip R. Clauer, and Howard M. Salis. "Simultaneous repression of multiple bacterial genes using nonrepetitive extra-long sgRNA arrays." Nature biotechnology 37, no. 11 (2019): 1294-1301.

Link to Article: https://www.nature.com/articles/s41587-019-0286-9
