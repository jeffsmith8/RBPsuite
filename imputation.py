# -*- coding: utf-8 -*-
"""
Redistributed Code

Copyright (c) 2015, Martin Fitzpatrick
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import numpy as np
try:
    import sklearn
except ImportError:
    sklearn = False
else:
    from sklearn.cross_decomposition import PLSRegression


def gaussian(df, width=0.3, downshift=-1.8, prefix=None):
    """
    Impute missing values by drawing from a normal distribution

    :param df:
    :param width: Scale factor for the imputed distribution relative to the standard deviation of measured values. Can be a single number or list of one per column.
    :param downshift: Shift the imputed values down, in units of std. dev. Can be a single number or list of one per column
    :param prefix: The column prefix for imputed columns
    :return:
    """

    df = df.copy()

    imputed = df.isnull()  # Keep track of what's real

    if prefix:
        mask = np.array([l.startswith(prefix) for l in df.columns.values])
        mycols = np.arange(0, df.shape[1])[mask]
    else:
        mycols = np.arange(0, df.shape[1])


    if type(width) is not list:
        width = [width] * len(mycols)

    elif len(mycols) != len(width):
        raise ValueError("Length of iterable 'width' does not match # of columns")

    if type(downshift) is not list:
        downshift = [downshift] * len(mycols)

    elif len(mycols) != len(downshift):
        raise ValueError("Length of iterable 'downshift' does not match # of columns")

    for i in mycols:
        data = df.iloc[:, i]
        mask = data.isnull().values
        mean = data.mean(axis=0)
        stddev = data.std(axis=0)

        m = mean + downshift[i]*stddev
        s = stddev*width[i]

        # Generate a list of random numbers for filling in
        values = np.random.normal(loc=m, scale=s, size=df.shape[0])

        # Now fill them in
        df.iloc[mask, i] = values[mask]

    return df, imputed


def pls(df):
    """
    A simple implementation of a least-squares approach to imputation using partial least squares
    regression (PLS).

    :param df:
    :return:
    """

    if not sklearn:
        assert('This library depends on scikit-learn (sklearn) to perform PLS-based imputation')

    df = df.copy()
    df[np.isinf(df)] = np.nan

    dfo = df.dropna(how='any', axis=0)
    dfo = dfo.astype(np.float64)
    
    dfi = df.copy()
    imputed = df.isnull() #Keep track of what's real

    # List of proteins with missing values in their rows
    missing_values = df[ np.sum(np.isnan(df), axis=1) > 0 ].index
    ix_mask = np.arange(0, df.shape[1])
    total_n = len(missing_values)

    #dfi = df.fillna(0)

    plsr = PLSRegression(n_components=2)

    for n, p in enumerate(missing_values.values):
        # Generate model for this protein from missing data
        target = df.loc[p].values.copy().T

        ixes = ix_mask[ np.isnan(target) ]

        # Fill missing values with row median for calculation
        target[np.isnan(target)] = np.nanmedian(target)
        plsr.fit(dfo.values.T, target)

        # For each missing value, calculate imputed value from the column data input
        for ix in ixes:
            imputv = plsr.predict(dfo.iloc[:, ix])[0]
            dfi.ix[p, ix] = imputv

        print("%d%%" % ((n/total_n)*100), end="\r")


    return dfi, imputed