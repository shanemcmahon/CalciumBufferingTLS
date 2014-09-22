CalciumBufferingTLS
===================

Source code to fit calcium buffering data to a 1-3 buffer saturable buffer model using total least squares.

Original source code is included as old.source.R, but the functions have been rewritten for clarity and to lessen the burden of making updates in the future. The newer source is also allows for fitting arbitrary, user specefied models, while the original version works only for the built in calcium buffering model.

The current version seemed to performed adequately on a number of test data sets. However, after applying it to a larger number of data sets it has become clear that the current version fails to find the minimum with significantly greater frequency than the original code. We are currently debugging this issue and hope to have it resolved soon.

Source code to accompany the paper "In Situ Ca2+ Titration in the Fluorometric Study of Intracellular Ca2+ Binding"

Includes a library for total least squares regression on arbitrary models using a stochastic gradient descent algorithm, an interface for performing the regression on a calcium buffering model, and an example data set.

To obtain the code, click on the "download zip" button in the right panel of the website https://github.com/shanemcmahon/CalciumBufferingTLS

The zip file contains:

A detailed user manual, "CalciumBufferingTLSUserManual.docx"

3 R files:

Setup.Ca.Buffering.TLS.R: installs required packages

Ca.Buffering.TLS.Sub.R: contains subroutines for performing TLS fitting

Run.Ca.Buffering.TLS.R: provides a user interface for fitting calcium buffering models

4 data files:

f.csv: initial and final fluorescence measurements

f.sem.csv: estimated standard errors for fluorescence measurements

ca.csv: change in total calcium concentration estimates obtained from integrated calcium currents

ca.sem.csv: estimated standard errors for calcium increments
