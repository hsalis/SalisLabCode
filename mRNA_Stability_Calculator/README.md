<b>mRNA Stability Calculator</b>

mRNA degradation is a central process that affects all gene expression levels, though it remains challenging to predict the stability of a mRNA from its sequence, due to the many coupled interactions that control degradation rate. Here, we carried out massively parallel kinetic decay measurements on 62120 bacterial mRNAs, using a learn-by-design approach to develop and validate a predictive sequence-to-function model of mRNA stability. mRNAs were designed to systematically vary translation rates, secondary structures, sequence compositions, G-quadruplexes, i-motifs, and RppH activity, resulting in mRNA half-lives from about 20 seconds to 20 minutes. We combined biophysical models and machine learning to develop steady-state and kinetic decay models of mRNA stability with high accuracy and generalizability, utilizing transcription rate models to identify mRNA isoforms and translation rate models to calculate ribosome protection. Overall, the developed model quantifies the key interactions that collectively control mRNA stability in bacterial operons and predicts how changing mRNA sequence alters mRNA stability, which is important when studying and engineering bacterial genetic systems.

Relevant Article: Cetnar, Daniel P., Ayaan Hossain, Grace E. Vezeau, and Howard M. Salis. Predicting synthetic mRNA stability using massively parallel kinetic measurements, biophysical modeling, and machine learning. in review

Link to Article:

<b>Instructions to Repeat the Train-Test Procedure</b>

Install Python Dependencies:
```python
pip3 install numpy pandas matplotlib seaborn scipy sklearn lightgbm biopython dm-sonnet sonnet graphs mpi4py openpyxl zipfile
```

Run Train-Test Procedure:
```python
python3 ModelTrainingTesting.py
```
