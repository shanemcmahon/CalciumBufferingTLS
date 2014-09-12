CalciumBufferingTLS
===================

Source code to fit calcium buffering data to a 1-3 buffer saturable buffer model using total least squares

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
