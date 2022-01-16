# SynBioMTS
An automated model test system for synthetic biology models of gene expression and regulation.
When using this code, remember to cite its publication: Reis, Alexander C., and Howard M. Salis. "An automated model test system for systematic development and improvement of gene expression models." ACS Synthetic Biology 9, no. 11 (2020): 3145-3156.
https://pubs.acs.org/doi/abs/10.1021/acssynbio.0c00394
Correspondence should be addressed to H.M.S. (salis@psu.edu)

Abstract:
Gene expression models greatly accelerate the engineering of synthetic metabolic pathways and genetic circuits by predicting sequence-function relationships and reducing trial-and-error experimentation. However, developing models with more accurate predictions remains a significant challenge. Here we present a model test system that combines advanced statistics, machine learning, and a database of 9862 characterized genetic systems to automatically quantify model accuracies, accept or reject mechanistic hypotheses, and identify areas for model improvement. We also introduce model capacity, a new information theoretic metric for correct cross-data-set comparisons. We demonstrate the model test system by comparing six models of translation initiation rate, evaluating 100 mechanistic hypotheses, and uncovering new sequence determinants that control protein expression levels. We then applied these results to develop a biophysical model of translation initiation rate with significant improvements in accuracy. Automated model test systems will dramatically accelerate the development of gene expression models, and thereby transition synthetic biology into a mature engineering discipline.

`synbiomts` uses a database of over 16000 unique characterized genetic systems to run Python-wrapped sequence-function models, quantify model accuracy, accept or reject proposed mechanistic hypotheses, and identify sources of model error. This package is easily modifiable to expand the genetic system database, calculate additional statistical test metrics, and test new and improved gene expression models implemented in nearly any programming language.

## Getting Started

### Dependencies
Python packages used are listed below. You can install the first three packages together as part of the [SciPy Stack](https://www.scipy.org/install.html).
* pandas - Database management
* scipy - Statistics calculations
* numpy - General purpose numerical computing
* [scikit-learn](http://scikit-learn.org/stable/install.html) - Machine learning

#### Optional

* [WebLogo](https://github.com/WebLogo/weblogo) - We've wrapped the WebLogo python API so that you can easily generate weblogos for your sequence datasets.

[ViennaRNA](https://www.tbi.univie.ac.at/RNA/) - A C code library for the prediction and comparison of RNA secondary structures. ViennaRNA is wrapped with /models/PyVRNA.py for use in modeling and machine learning analysis.

### Installing
Install with the following:
```
git clone https://github.com/reisalex/SynBioMTS
cd SynBioMTS
sudo python setup.py install
```
The model test system can then be imported in Python:
```python
import synbiomts
```

### Usage
If you would like to use the provided genetic system database, the best way is to navigate to /synbiomts, and run the database initialization module (initdb.py):
```
cd synbiomts/examples/RBS
python initdb.py
```

To use the model test system:
1. Wrap the model with a Python function.
2. Create a models `Container` object and pass the wrapped functions with the `add` method.
3. Specify the functional form between the model predictor and the system function with the `setform` method.
4. Create a `ModelTest` object and pass the models Container.
5. Run model calculations and statistics with `run`.

```python
import synbiomts

# Wrap the model with a function
import RBS_Calculator_v2
def RBSCalcv2(sequence,temperature):
    rRNA = 'ACCTCCTTA'
    model = RBS_Calculator_v2.RBS_Calculator(mRNA=sequence,rRNA=rRNA)
    model.temp = temperature
    model.run()
    RBS = model.output() # simplified for the example

    # Results should be returned as a dictionary
    # The keys will become labels in the resulting pandas dataframe
    results = {
    'TIR': RBS.tir,
    'dG_total': RBS.dG_total,
    'dG_mRNA_rRNA': RBS.dG_mRNA_rRNA',
    'dG_mRNA': RBS.dG_mRNA
    }
    return results

if __name__ == "__main__":

    # create models Container object
    models = synbiomts.interface.Container()

    # add the model(s)
    models.add(RBSCalcv2)

    # specify the form of each model
    # RBS_Calculator is a thermodynamic model where: Protein ~ K*exp(-0.45*dG_total)
    models.setform(['RBSCalcv2'],x='dG_total',y='PROT.MEAN',yScale='ln',a1=-0.45)

    # create test system object
    testsystem = synbiomts.analyze.ModelTest(models,'geneticsystems.db',add_data=True,verbose=True)

    # run model predictions and statistics calculations
    testsystem.run()

    # if you want to shelve the model calculations
    # testsystem.run(calcsFilename='savedcalcs.db')
```

When you add models to the `Containers` object, you can specify arguments of the wrapped function. This comes in handy when you want to vary a parameter and test which is most accurate:
```python
for s in range(0,16):
    name = "RBSCalcv2-s={}.format(s)
    models.add(RBSCalcv2,optimal_spacing=s)
```

You can specify filters to run predictions on a subset of genetic systems with shared properties:
```python
filters = { 'ORGANISM': ['Escherichia coli'],
            'DATASET' : ['Beck_PLoS_2016',
                         'Salis_NBT_2009',
                         'Tian_NAR_2015']
}
testsystem = synbiomts.analyze.ModelTest(models,'geneticsystems.db',filters)
```

The model test system uses multiprocessing, with the number of available CPUs by default, to run model predictions. You can specify the number of processes to force single process or specify a desired number:
```python
testsystem = synbiomts.analyze.ModelTest(models,'geneticsystems.db',nprocesses=1)
```

The `run` method calculates both model predictions and calculates statistics. If you only want to run model predictions, you can use `predict`:
```python
testsystem.predict()
```

See /examples for more detailed examples.

### Statistics
By default, the model test system will run statistics assuming the model predictor and the system function share a linear relationship. Specifically, `analyze.statistics()` calls a custom function `linear_complete` from the stats module to compute the following:
* Fitted slope and y-intercept with outliers removed (via MAD method)
* Relative model error (Apparent Value/Predicted Value)
* Pearson & Spearman correlation coefficients
* One-sided model error cummulative distribution function
* Kullback-Leibler divergence

If you want to run futher statistics, you can import the stats module:
```python
from synbiomts import stats
```
You can always add additional stats functions as needed.

### Exporting
Export to Excel is as simple as:
```python
testsystem.run()
testsystem.to_excel('filename')
```
By default pandas exports with the labels (columns) alphabetized. The model test system overrides the default export if you specify the labels. See /examples/labels for the ones I use:
```python
test.run()

with open("labels/labels1.txt","r") as f:
    predictLabels = [x.strip('\n') for x in f.readlines()]

with open("labels/labels_stats.txt","r") as f:
    statsLabels = [x.strip('\n') for x in f.readlines()]

test.to_excel('filename',predictLabels,statsLabels)
```

## Acknowledgements

Thanks to Howard M Salis (Penn State), Iman Farasat (Merck), Amin Espah Borujeni (MIT), Tian Tian (JBEI), Daniel Goodman (Harvard), Sri Kosuri (UCLA), Robert Egbert (Berkeley), Mark Mimee (MIT), and Heather Beck (Vienna) for providing high quality characterization data. A special thanks to Daniel Goodman for discussion on Flow-seq and for providing additional information on the 2013 Flow-seq datasets.

If you use `synbiomts`, please cite:

Alexander C. Reis, and Howard M. Salis. An automated model test system for systematic development and improvement of gene expression models, In Preparation (2017).
