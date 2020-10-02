"""
Custom statistics module for SynBioMTS

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

"""

from __future__ import division
import re
import numpy as np
import scipy
import scipy.stats

def correlation(x,y,name="Pearson"):
    # Linear or rank correlation
    # Returns:
    #   r,rho, or tau (float) = the test statistic
    #   pvalue (float) = the p-value for the test

    assert len(x)==len(y), "Arrays x & y must be equal length."

    if name == "Pearson":
        (r,pvalue) = scipy.stats.pearsonr(x,y)
        return r,pvalue

    elif name == "Spearman":
        (rho,pvalue) = scipy.stats.spearmanr(x,y)
        return rho,pvalue

    elif name == "Kendall":
        (tau,pvalue) = scipy.stats.kendalltau(x,y)
        return tau,pvalue

    else:
        error = ("The {} correlation is not available."
                 "Please use 'Pearson', 'Spearman', or 'Kendall.'")
        error.format(name)
        raise ValueError(error)

def vartest2(x,y,logNormal=False,alpha=0.05,test="F"):
    '''Two-sample F-test/Barlett-test/Levene-test for equal variances
    Returns:
        h (bool)  = True if the test rejects the null hypothesis
        S (float) = the test statistic
        pvalue (float) = the p-value for significance of test decision'''

    if logNormal:
        varx = np.exp(np.var(np.log(x))) # an estimate
        vary = np.exp(np.var(np.log(y)))
        # varx = np.var(np.log(x))
        # vary = np.var(np.log(y))
    else:
        varx = np.var(x)
        vary = np.var(y)

    if test == "F":
        statistic = varx/vary
        if statistic < 1:
            df1 = len(y)-1
            df2 = len(x)-1
            F = 1/statistic
        else:
            df1 = len(x)-1
            df2 = len(y)-1            
            F = statistic
        pvalue = scipy.stats.f.sf(F, df1, df2)
        if F == 1.0:
            pvalue = 0.0

    elif test == "Barlett":
        (statistic,pvalue) = scipy.stats.bartlett(x,y)

    elif test == "Levene":
        (statistic,pvalue) = scipy.stats.levene(x,y)
    else:
        error = "The {}-test for equal variances is not available. \
        Please use 'F', 'Barlett', or 'Levene'.".format(test)
        raise ValueError(error)

    h = pvalue < alpha
    # h = True, if we reject the null hypothesis that x & y come from
    # normal distributions with the same variance
    return (h,statistic,pvalue)

def ttest2(x,y,alpha=0.05):
    '''Two sample t-test
    Returns:
        h (bool)  = True if the test rejects the null hypothesis
        t (float) = the test statistic
        pvalue (float) = the p-value for significance of test decision'''

    (t,pvalue) = scipy.stats.ttest_ind(x,y,equal_var=False,nan_policy='omit')
    h = pvalue < alpha
    return (h,t,pvalue)

def fit_linear_model(x,y,slope=None):
    '''Linear least squares (LSQ) linear regression (polynomial fit of degree=1)
    Returns:
        m (float) = slope of linear regression line of form (y = m*x + b)
        b (float) = intercept of linear regression line'''

    assert len(x)==len(y), ("Arrays x & Y must be equal length to fit "
                            "linear regression model.")
    if slope == None:
        (m,b) = scipy.polyfit(x,y,deg=1)
    else:
        LSQ = lambda b: np.sum( (y-(slope*x+b))**2.0 )
        res = scipy.optimize.minimize(LSQ,x0=1,bounds=None)
        (m,b) = (slope,res.x[0])
    return (m,b)

def mad(x):
    '''Median absolute deviation
    Returns:
        MAD (float) = the median absolute deviation of the values in X
        diff (list) = a list of deviations of the values in X from the median of X'''
    median = np.median(x)
    diff = np.absolute(x - median)
    MAD = np.median(diff)
    return (MAD,diff)

''' Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
    Handle Outliers", The ASQC Basic References in Quality Control:
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.'''  
def find_outliers(x,threshold=2):
    '''Outlier detection method using Iglewicz & Hoaglin's modified Z-score
    Returns a list of bools, where True is an outlier'''
    (MAD,diff) = mad(x)
    M = 0.6745 * diff / MAD
    outliers = M > threshold
    return outliers

def empirical_cdf(data,bins,rangemin=None,rangemax=None):
    '''Calculate the empirical cumulative distribution function for data
    Inputs:
        data (array) = data to calculate ecdf of
        bins (int/array) = same behavior as numpy.histogram (int defines
                           number of bins, array specifies bin edges)
        rangemin/rangemax (float) = must be specified if bins is of type=int,
                                    and rangemin < rangemax
    Returns:
        ecdf (array) = empirical cumulative distribution function
        edges (array) = edges used during binning of histogram'''
    
    if isinstance(bins,np.ndarray):
        edges = bins
    else:        
        message = "bins must either be an int or a list corresponding to the edges"
        assert any(isinstance(bins,t) for t in [int,list]), message

        if isinstance(bins,int):
            assert bins > 1, "number of bins must be greater than 1."
            assert rangemin < rangemax, "rangemin must be less than rangemax"
            assert rangemin is not None and rangemax is not None, ("rangemin and "
                "rangemax must be specified if bins is an integer.")

            edges = np.linspace(rangemin,rangemax,bins+1)

        elif isinstance(bins,list):
            edges = np.array(bins)
    
    h = np.histogram(data,edges)
    ecdf = np.cumsum(h[0])/len(data)

    return ecdf,edges

def area_under_ROC_curve(predictions,observations,cutoff=None,n=50):
    '''Calculate ROC curve and the area under the curve (auc)

    Returns:
        auroc (float) = area under the ROC curve
        fpr (array) = false positive rates
        tpr (array) =  true positive rates
        thresholds (array) = thresholds from ROC curve'''

    # define outcomes as TRUE/FALSE defined by cutoff
    if cutoff == None: cutoff = np.median(outcomes)
    outcomes = (observations > cutoff)

    # preallocate arrays and define n-many threshodls
    fpr = np.zeros(n)
    tpr = np.zeros(n)
    thresholds = np.linspace(0.9*min(observations),1.1*max(observations),n)

    for i in xrange(0,n):

        threshold = thresholds[i]
        scores = predictions > threshold

        # calculate confusion matrix
        a = sum(scores & outcomes)
        b = sum(scores & ~outcomes)
        c = sum(~scores & outcomes)
        d = sum(~scores & ~outcomes)            

        # calculate fpr and tpr for each threshold
        fpr[i] = b/(b+d)
        tpr[i] = a/(a+c)

    # Remove nans
    nans = np.isnan(fpr) + np.isnan(tpr)
    x = np.array(fpr)
    y = np.array(tpr)
    x = x[~nans]
    y = y[~nans]

    # Sort by increasing x
    indx = np.argsort(x)
    x = x[indx]
    y = y[indx]

    # Remove repeated elements
    x,indx = np.unique(x,return_index=True)
    y = y[indx]

    # Add (0,0) or (1,1) if not present!
    if x[0]!=0.0:
        x = np.append([0.0],x)
        y = np.append([0.0],y)

    if x[-1]!=1.0:
        x = np.append(x,1.0)
        y = np.append(y,1.0)

    # Calculate area under the ROC curve
    fpr = x
    tpr = y
    auroc = np.trapz(tpr,fpr)

    return auroc,fpr,tpr,thresholds

def normKLdiv(data,b):
    '''Calculate normalized Kullback-Leibler divergence

    normKLdiv uses data and calculates the normKLdiv relative to a
    random model (Q) over bin range (-b,b).
    Returns:
        NKLdiv (float)   = normalized Kullback-Leibler divergence
        KLdiv (float)    = KL-divergence
        KLdivmax (float) = KL-divergence for a perfect model (dirac delta fxn)'''

    (P,x) = pdist(data,b)
    Pmax = np.zeros(len(x))
    Pmax[0] = 1.0
    Q = np.ones(len(x))/len(x)
    KLdiv = entropy(P,Q)
    KLdivmax = entropy(Pmax,Q)
    NKLdiv = KLdiv/KLdivmax

    return (NKLdiv,KLdiv,KLdivmax)

def pdist(data,b,nbins=100,make_outliers_rand=True):
    '''Calculate the discrete probability distribution function
    of a set of data defined over the range (-b,b) with nbins
    Inputs:
        ignore_outliers (bool) = if True ignore outliers,
                                 if False add random information to pdf
    Returns:
        pdf (array) = discrete pdf with probability of event i as pdf[i]
        x (array) = midpoints of the bins of the histogram that defines pdf'''

    # remove outliers outside of range (-b,b)
    data2 = data[~(np.abs(data) > b)]

    # calculate Pset (the pdf of the set of data points within the range)
    edges = np.linspace(-b,b,nbins+1)
    Pset,edges = np.histogram(data2,edges,density=True)
    fset = len(data2)/len(data)
    x = (edges[0:-1]+edges[1:])/2

    if make_outliers_rand:
        # calculate Poutliers (pdf corresponding to random guesses)
        Pout = np.ones(len(x))/len(x)
        fout = 1 - fset
        pdf = fset*Pset + fout*Pout
    else:
        pdf = Pset

    return (pdf,x)

def entropy(pk,qk=None,base=None):
    '''Calculate the entropy of a distribution for given probability values

    Inputs:
        pk (array)   = probability distribution of a discrete distribution
        qk (array)   = pdf of a second sequence with which pk is compared against
        base (float) = the logarithmic base to use (default=e)
    Returns:
        S (float) = Shannon entropy if qk is None else Kullback-Leibler divergence'''
    return scipy.stats.entropy(pk,qk,base)

def model_capacity(Hseq,mu,sigma,error):
    '''Calculate model capacity (MC)

    Inputs:
        Hseq (float)       = total Shannon sequence entropy for the genetic system
        mu, sigma (arrays) = mean and standard deviations of system outcomes with mu[i]
                             and sigma[i] as the mean and std of outcome i
        error (array)      = model error for system outcomes (of equal length of mu & sigma)
    Returns:
        MC (float)      = model capacity in bits
        Hmodel (float)  = Shannon entropy of the model error distribution
        Hrandom (float) = Shannon entropy of a random model error distribution
        N (float)       = number of distinguishable outcomes
        CV (float)      = coefficient of variation of system outcomes'''

    # calculate N-distinguishable outcomes
    CV = mean( sigma/mu )
    if CV == 0.0:
        N = 1000
    else:
        N = np.round((max(np.log10(mu))-min(np.log10(mu)))/np.log10(1+CV))

    # calculate Hmodel
    pdf,_ = pdist(error,b=4,nbins=N)
    Hmodel = entropy(pdf)
    
    # calculate Hrandom
    pdf,_ = pdist(np.array([np.inf]),b=4,nbins=N)
    Hrandom = entropy(pdf)

    # calculate MC
    MC = Hseq*(N-1)*(1-Hmodel/Hrandom)

    return (MC,Hmodel,Hrandom,N,CV)

def sequence_entropy(sequences,align="left",positions=None):
    '''Calculate total Shannon sequence entropy

    Inputs:
        sequences (list) = sequences in a list
        align (string)   = "left"/"right"
                           "positions" indicates to use positions for alignment
        positions (list) = list of integers to do alignment
                            (e.g. start codon positions)
    Returns:
        hseq (float) = total Shannon sequence entropy summed over all positions
        S (array)    = Shannon sequence entropy at each position'''

    sequences = [seq.upper().replace("U","T") for seq in sequences[:]]
    exp = re.compile('[ATGC]')
    if any([exp.match(seq) is None for seq in sequences]):
        raise ValueError("Invalid letters found in sequences. Only ATGCU accepted.")

    nseqs = len(sequences)
    maxseqlen = max([len(seq) for seq in sequences])
    samelen   = all([len(sequences[0])==len(seq) for seq in sequences[1:]])

    if   align == "left" and not samelen:
        # align sequences to left and buffer right with Xs
        sequences = [seq+"X"*(maxseqlen-len(seq)) for seq in sequences[:]]

    elif align == "right" and not samelen:
        # align sequences to the right by buffering left with Xs
        sequences = ["X"*(maxseqlen-len(seq))+seq for seq in sequences[:]]

    elif not positions is None:
        # align sequences at alignment_positions and buffer both ends with Xs
        lefts     = [len(seq[0:pos]) for seq,pos in zip(sequences,positions)]
        rights    = [len(seq[pos:]) for seq,pos in zip(sequences,positions)]
        maxleft   = max(lefts)
        maxright  = max(rights)
        sequences = ["X"*(maxleft-L)+seq+"X"*(maxright-R) for L,seq,R in zip(lefts,sequences[:],rights)]
        
    else:
        raise ValueError("align cannot be {}. Please use 'left','right',or positions".format(align))

    # Code test
    assert all(len(sequences[0])==len(seq) for seq in sequences[1:])

    # tabulate frequency of each nt at each position
    if samelen: alphabet = 'ATCG'
    else:       alphabet = 'ATCGX'
    counts = {key: np.zeros(maxseqlen) for key in alphabet}

    for seq in sequences:
        for i,c in enumerate(seq):
            counts[c][i] += 1
    
    pkList = [[counts[c][i]/nseqs for c in alphabet] for i in xrange(maxseqlen)]

    # calculate Shannon entropy at each position and total Shannon entropy (hseq)
    S = [entropy(pk,base=2) for pk in pkList]
    Hseq = sum(S)

    return Hseq,S

def linear_complete(xVals,yVals,ystd,xScale='linear',yScale='linear',slope=None):

    # Useful lambda functions
    calc_x = lambda a0,a1,y: (y-a0)/a1
    calc_y = lambda a0,a1,x: a1*x+a0

    if   yScale == 'log10':  yVals = np.log10(yVals)
    elif yScale == 'ln':     yVals = np.log(yVals)
    elif yScale == 'linear': pass
    else: raise ValueError("Invalid input in ModelTest._linear_model_stats for yScale: {}".format(yScale))

    if   xScale == 'log10':  xVals = np.log10(xVals)
    elif xScale == 'ln':     xVals = np.log(xVals)
    elif xScale == 'linear': pass
    else: raise ValueError("Invalid input in ModelTest._linear_model_stats for xScale: {}".format(xScale))    

    # determine outliers with initial fit
    (a1,a0) = fit_linear_model(xVals,yVals,slope=slope)
    app_xVals = calc_x(a0,a1,yVals)
    abs_delta_x = np.absolute(xVals - app_xVals)
    outliers = find_outliers(abs_delta_x)
    keepers = np.invert(outliers)

    # reapply fit with trimmed dataset & calculate error
    xVals1 = xVals[keepers]
    yVals1 = yVals[keepers]
    (a1,a0) = fit_linear_model(xVals1,yVals1,slope=slope)
    app_xVals = calc_x(a0,a1,yVals)
    y_predicted = calc_y(a0,a1,xVals)

    residuals = yVals - y_predicted

    if yScale == 'log10': yError = 10**(residuals)
    elif yScale == 'ln':  yError = np.exp(residuals)
    else:                 yError = yVals/y_predicted

    # Pearson/Spearman correlation coefficients
    (R,Pearson_p) = correlation(y_predicted,yVals)
    (rho,Spearman_p) = correlation(y_predicted,yVals)

    # Root Mean Square Error (RMSE)
    # n = len(yVals)
    # SSE = np.sum(res**2.0 for res in residuals)
    # RMSE = np.sqrt(SSE/n)
    RMSE = np.sqrt(1-R**2.0)*np.std(yVals)

    # One-sided model error cdfs
    yError1 = np.array([1/val if val < 1 else val for val in yError])
    bins = np.concatenate((np.linspace(1,10,10),np.linspace(20,100,9),np.linspace(200,1000,9)))
    onesided_cdf,_ = empirical_cdf(yError1,bins)

    # Kullback-Leibler divergence
    (NKLdiv,KLdiv,KLdivmax) = normKLdiv(yError,b=4)

    # AUC ROC
    threshold = (max(yVals) + min(yVals))/2.0
    auroc,fpr,tpr,thresholds = area_under_ROC_curve(y_predicted,yVals,cutoff=threshold)

    # Relative information gain (RIG) over uniform model
    # remove nans from ystd
    # if len(ystd)==0, then skip information theory analysis
    ystd = np.array(ystd, dtype=np.float64)
    nans = np.isnan(ystd)
    ystd = ystd[~nans]
    yVals = yVals[~nans]

    if len(ystd)==0:
        CV = 0.0
        N = 0.0
        RIG = 0.0
    else:
        if   yScale == 'log10': ystd = np.log10(ystd)
        elif yScale == 'ln':    ystd = np.log(ystd)
        else: pass

        nonzero = np.nonzero(ystd*yVals)
        ystd = ystd[nonzero]
        yVals = yVals[nonzero]

        negative = (yVals < 0.0) + (ystd < 0.0)
        CV = np.mean(ystd[~negative]/yVals[~negative])
        N = int(np.floor((max(yVals) - min(yVals))/CV))
        edges = np.linspace(-4,4,N+1)
        # filter out residuals that don't fall in the bins
        rmv = (residuals < edges[0]) + (residuals > edges[-1])
        residuals = residuals[~rmv]

        [dist,_] = np.histogram(residuals,edges,density=True)
        dist_model = dist[dist!=0.0]
        Hmodel = entropy(dist_model)
        dist = np.ones(N)/N;
        Hrandom = entropy(dist)
        RIG = 1 - Hmodel/Hrandom
    # WARNING yVals and y_predicted have been modified after MC calculations

    results = {
    "Pearson R-squared": R**2.0,
    "Pearson p": Pearson_p,
    "Spearman R-squared": rho**2.0,
    "Spearman p": Spearman_p,
    "N-outliers": np.sum(outliers),
    "slope": a1,
    "intercept": a0,
    "2-fold Error": onesided_cdf[1],
    "5-fold Error": onesided_cdf[4],
    "10-fold Error": onesided_cdf[9],
    "Normalized KL-divergence": NKLdiv,
    "KL-divergence": KLdiv,
    "Max KL-divergence": KLdivmax,
    "RMSE": RMSE,
    "RIG": RIG,
    "N-states": N,
    "CV": CV,
    "AUROC": auroc
    }

    return results,yError

def linear_simple(xVals,yVals,ystd,xScale='linear',yScale='linear'):

    yError = yVals/xVals

    if yScale == 'log10':    yVals = np.log10(yVals)
    elif yScale == 'ln':     yVals = np.log(yVals)
    elif yScale == 'linear': pass
    else: raise ValueError("Invalid input in ModelTest._linear_model_stats for yScale: {}".format(yScale))

    if xScale == 'log10':    xVals = np.log10(xVals)
    elif xScale == 'ln':     xVals = np.log(xVals)
    elif xScale == 'linear': pass
    else: raise ValueError("Invalid input in ModelTest._linear_model_stats for xScale: {}".format(xScale))    

    # Pearson/Spearman correlation coefficients
    (R,Pearson_p) = correlation(xVals,yVals)
    (rho,Spearman_p) = correlation(xVals,yVals)

    # Root Mean Square Error (RMSE)
    residuals = yVals - xVals
    # n = len(yVals)
    # SSE = np.sum(res**2.0 for res in residuals)
    # RMSE = np.sqrt(SSE/n)
    RMSE = np.sqrt(1-R**2.0)*np.std(yVals)

    # One-sided model error cdfs
    yError1 = np.array([1/val if val < 1 else val for val in yError])
    bins = np.concatenate((np.linspace(1,10,10),np.linspace(20,100,9),np.linspace(200,1000,9)))
    onesided_cdf,_ = empirical_cdf(yError1,bins)

    # Kullback-Leibler divergence
    (NKLdiv,KLdiv,KLdivmax) = normKLdiv(yError,b=4)

    # AUC ROC
    threshold = (max(yVals) + min(yVals))/2.0
    auroc,fpr,tpr,thresholds = area_under_ROC_curve(xVals,yVals,cutoff=threshold)

    # Relative entropy gain over uniform model
    # Relative information gain (RIG) over uniform model
    # remove nans from ystd
    # if len(ystd)==0, then skip information theory analysis
    ystd = np.array(ystd, dtype=np.float64)
    nans = np.isnan(ystd)
    ystd = ystd[~nans]
    yVals = yVals[~nans]
    
    if len(ystd)==0:
        CV = 0.0
        N = 0.0
        RIG = 0.0
    else:
        if   yScale == 'log10': ystd = np.log10(ystd)
        elif yScale == 'ln':    ystd = np.log(ystd)
        else: pass

        nonzero = np.nonzero(ystd*yVals)
        ystd = ystd[nonzero]
        yVals = yVals[nonzero]

        negative = (yVals < 0.0) + (ystd < 0.0)
        CV = np.mean(ystd[~negative]/yVals[~negative])
        N = int(np.floor((max(yVals) - min(yVals))/CV))
        edges = np.linspace(-4,4,N+1)
        # filter out residuals that don't fall in the bins
        rmv = (residuals < edges[0]) + (residuals > edges[-1])
        residuals = residuals[~rmv]

        [dist,_] = np.histogram(residuals,edges,density=True)
        dist_model = dist[dist!=0.0]
        Hmodel = entropy(dist_model)
        dist = np.ones(N)/N;
        Hrandom = entropy(dist)
        RIG = 1 - Hmodel/Hrandom
    # WARNING yVals and y_predicted have been modified after MC calculations

    # Number of outliers
    outliers = find_outliers(yError)

    results = {
    "Pearson R-squared": R**2.0,
    "Pearson p": Pearson_p,
    "Spearman R-squared": rho**2.0,
    "Spearman p": Spearman_p,
    "N-outliers": np.sum(outliers),
    "2-fold Error": onesided_cdf[1],
    "5-fold Error": onesided_cdf[4],
    "10-fold Error": onesided_cdf[9],
    "Normalized KL-divergence": NKLdiv,
    "KL-divergence": KLdiv,
    "Max KL-divergence": KLdivmax,
    "RMSE": RMSE,
    "RIG": RIG,
    "N-states": N,
    "CV": CV,
    "AUROC": auroc
    }

    return results,yError

if __name__ == "__main__":
    pass