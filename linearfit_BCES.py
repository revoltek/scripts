#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""
Bivariate, Correlated Errors and intrinsic Scatter (BCES)
    translated from the FORTRAN code by Christina Bird and Matthew Bershady
    (Akritas & Bershady, 1996)

Do linear regression including error on both variables and intrinsic scatter.

"""

import scipy

def bces(x1, x2, x1err=None, x2err=None, cerr=None, nsim=5000, model='yx', \
         bootstrap=5000, verbose='normal', full_output=True):
  """
  Bivariate, Correlated Errors and intrinsic Scatter (BCES)
    translated from the FORTRAN code by Christina Bird and Matthew Bershady
    (Akritas & Bershady, 1996)

  Linear regression in the presence of heteroscedastic errors on both
  variables and intrinsic scatter

  Parameters
  ----------
    x1        : array of floats
                Independent variable, or observable
    x2        : array of floats
                Dependent variable
    x1err     : array of floats (optional)
                Uncertainties on the independent variable
    x2err     : array of floats (optional)
                Uncertainties on the dependent variable
    cerr      : array of floats (optional)
                Covariances of the uncertainties in the dependent and
                independent variables
    nsim      : int (default 1000)
                Number of bootstrap samples for uncertainties on best-fit
                parameters
    model     : {'yx', 'xy', 'bi', 'orth'}
                BCES model with which to calculate regression. See Notes
                below for details.
    bootstrap : False or int (default False)
                get the errors from bootstrap resampling instead of the
                analytical prescription? if bootstrap is an int, it is the
                number of bootstrap resamplings
    verbose   : str (default 'normal')
                Verbose level. Options are {'quiet', 'normal', 'debug'}
    full_output : bool (default True)
                If True, return also the covariance between the normalization
                and slope of the regression.

  Returns
  -------
    a         : tuple of length 2
                Best-fit normalization and its uncertainty (a, da)
    b         : tuple of length 2
                Best-fit slope and its uncertainty (b, db)

  Optional outputs
  ----------------
    cov       : 2x2 array of floats
                covariance between a and b. Returned if full_output is set to
                True.

  Notes
  -----
    If verbose is normal or debug, the results from all the BCES models will
    be printed (still, only the one selected in *model* will be returned).

    the *model* parameter:
      -'yx' stands for BCES(Y|X)
      -'xy' stands for BCES(X|Y)
      -'bi' stands for BCES Bisector
      -'orth' stands for BCES Orthogonal

  """
  
  def _bess(npts, x1, x2, x1err, x2err, cerr):
    """
    Do the entire regression calculation for 4 slopes:
      OLS(Y|X), OLS(X|Y), bisector, orthogonal
    """

    # calculate sigma's for datapoints using length of confidence intervals
    sig11var = sum(x1err ** 2) / npts
    sig22var = sum(x2err ** 2) / npts
    sig12var = sum(cerr) / npts

    # calculate means and variances
    x1av = scipy.average(x1)
    x1var = scipy.std(x1) ** 2
    x2av = scipy.average(x2)
    x2var = scipy.std(x2) ** 2
    covar_x1x2 = sum((x1 - x1av) * (x2 - x2av)) / npts

    # compute the regression slopes for OLS(X2|X1), OLS(X1|X2), 
    # bisector and orthogonal
    b = scipy.zeros(4)
    b[0] = (covar_x1x2 - sig12var) / (x1var - sig11var)
    b[1] = (x2var - sig22var) / (covar_x1x2 - sig12var)
    b[2] = (b[0] * b[1] - 1 + scipy.sqrt((1 + b[0] ** 2) * \
           (1 + b[1] ** 2))) / (b[0] + b[1])
    b[3] = 0.5 * ((b[1] - 1 / b[0]) + scipy.sign(covar_x1x2) * \
           scipy.sqrt(4 + (b[1] - 1 / b[0]) ** 2))

    # compute intercepts for above 4 cases:
    a = x2av - b * x1av

    # set up variables to calculate standard deviations of slope and intercept
    xi = []
    xi.append(((x1 - x1av) * (x2 - b[0] * x1 - a[0]) + b[0] * x1err ** 2) / \
              (x1var - sig11var))
    xi.append(((x2 - x2av) * (x2 - b[1] * x1 - a[1]) + x2err ** 2) / \
              covar_x1x2)
    xi.append((xi[0] * (1 + b[1] ** 2) + xi[1] * (1 + b[0] ** 2)) / \
              ((b[0] + b[1]) * scipy.sqrt((1 + b[0] ** 2) * (1 + b[1] ** 2))))
    xi.append((xi[0] / b[0] ** 2 + xi[1]) * b[3] / \
              scipy.sqrt(4 + (b[1] - 1 / b[0]) ** 2))
    zeta = []
    for i in range(4):
      zeta.append(x2 - b[i] * x1 - x1av * xi[i])

    # calculate  variance for all a and b
    bvar = scipy.zeros(4)
    avar = scipy.zeros(4)
    for i in range(4):
      bvar[i] = scipy.std(xi[i]) ** 2 / npts
      avar[i] = scipy.std(zeta[i]) ** 2 / npts

    return a, b, avar, bvar, xi, zeta

  def _bootspbec(npts, x, y, xerr, yerr, cerr):
    """
    Bootstrap samples
    """
    j = scipy.random.randint(npts, size = npts)
    xboot = x[j]
    xerrboot = xerr[j]
    yboot = y[j]
    yerrboot = yerr[j]
    cerrboot = cerr[j]
    return xboot, yboot, xerrboot, yerrboot, cerrboot

  # ----  Main routine starts here  ---- #

  models = [['yx', 'xy', 'bi', 'orth'],
            ['BCES(Y|X)', 'BCES(X|Y)', 'BCES Bisector', 'BCES Orthogonal']]
  # which to return?
  j = models[0].index(model)

  npts = len(x1)
  # are the errors defined?
  if x1err is None:
    x1err = scipy.zeros(npts)
  if x2err is None:
    x2err = scipy.zeros(npts)
  if cerr is None:
    from scipy import random
    cerr = scipy.zeros(npts)
    #cerr = scipy.cov(x1err, x2err)[1][0] * scipy.ones(npts)

  if verbose == 'debug':
    print('x1 =', x1)
    print('x1err =', x1err)
    print('x2 =', x2)
    print('x2err =', x2err)
    print('cerr =', cerr)
    print('\n ** Returning values for', models[1][j], '**')
    if bootstrap is not False:
      print('    with errors from %d bootstrap resamplings' %bootstrap)
    print('')

  # calculate nominal fits
  bessresults = _bess(npts, x1, x2, x1err, x2err, cerr)
  (a, b, avar, bvar, xi, zeta) = bessresults
  # covariance between normalization and slope
  if full_output:
    covar_ab = scipy.cov(xi[j], zeta[j])

  if bootstrap is not False:
    # make bootstrap simulated datasets, and compute averages and
    # standard deviations of regression coefficients
    asum = scipy.zeros(4)
    assum = scipy.zeros(4)
    bsum = scipy.zeros(4)
    bssum = scipy.zeros(4)
    sda = scipy.zeros(4)
    sdb = scipy.zeros(4)
    for i in range(nsim):
      samples = _bootspbec(npts, x1, x2, x1err, x2err, cerr)
      (x1sim, x2sim, x1errsim, x2errsim, cerrsim) = samples
      besssim = _bess(npts, x1sim, x2sim, x1errsim, x2errsim, cerrsim)
      (asim, bsim, avarsim, bvarsim, xi, zeta) = besssim
      asum += asim
      assum += asim ** 2
      bsum += bsim
      bssum += bsim ** 2

    aavg = asum / nsim
    bavg = bsum / nsim
    for i in range(4):
      sdtest = assum[i] - nsim * aavg[i] ** 2
      if sdtest > 0:
        sda[i] = scipy.sqrt(sdtest / (nsim - 1))
      sdtest = bssum[i] - nsim * bavg[i] ** 2
      if sdtest > 0:
        sdb[i] = scipy.sqrt(sdtest / (nsim - 1))

  if verbose in ('normal', 'debug'):
    print('%s   B          err(B)' %('Fit'.ljust(19)), end=' ')
    print('         A          err(A)')
    for i in range(4):
      print('%s  %9.2e +/- %8.2e    %10.3e +/- %9.3e' \
            %(models[1][i].ljust(16), b[i], 
              scipy.sqrt(bvar[i]), a[i], scipy.sqrt(avar[i])))
      if bootstrap is not False:
        print('%s  %9.2e +/- %8.2e    %10.3e +/- %9.3e' \
              %('bootstrap'.ljust(16), bavg[i], sdb[i], aavg[i], sda[i]))
      print('')
    if verbose == 'debug':
      print('cov[%s] =' %models[model])
      print(covar_ab)

  if bootstrap is not False:
    if full_output:
      return (a[j], sda[j]), (b[j], sdb[j]), covar_ab
    else:
      return (a[j], sda[j]), (b[j], sdb[j])

  if full_output:
    return (a[j], scipy.sqrt(avar[j])), (b[j], scipy.sqrt(bvar[j])), covar_ab
  else:
    return (a[j], scipy.sqrt(avar[j])), (b[j], scipy.sqrt(bvar[j]))

