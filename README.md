# Score-Bootstrap-for-Cox-Regression-Models

Here we provide R code to perform Cox regressino where we have explicitly account for time-dependent weights.  Coxtdw_tests.pdf provides the logic and notation used in the code for this and provides a few tests to verify that the code is correct. File coxtdw.R is the main code to perform the Cox regression and coxtdw.utils.R contains utility code for computing terms needed in calculations. Have the explicit score and information matrix facilitates implementing the score bootstrap for computing standard errors. File scoreBSE.r contains code for computing bootstrapped standard errors using teh score bootstrap.  This method may have marginally lower coverage than using a non-parametric bootstrap, but it is so much faster within the context of a large simultaion study.
