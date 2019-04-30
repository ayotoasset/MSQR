# MSQR 
## Markov-Switching-(non-)linear-Quantile-Regression

With permission of Prof. Thomas Kneib (Uni Göttingen) and Timo Adams (Uni Bielefeld), I, Tim Ruhkopf (Uni Göttingen) publish my work at the Chair of Statistics and Econometrics of the Georg-August-University Göttingen.

Please take a look at the Use-cases first, to get an impression on the supported interface. 
Note that MS.drawy is a very generic function, supporting diverse data generating scenarios, that go far beyond linear data. To generate non linear data from the GAMLSS family, consider the nested state_formula structure, that holds the state regimes and describes for each state the functions/constants on the GAMLSS distribution's parameters of the supplied yfamily.

One major drawback of this work is the RQSS packages' incapability to handle the weights properly, prohibiting Timo Adam's (refactored) Algorithm MS.fitquantreg to be enhanced to handle non-linear quantile regression. However, this package is non-linear-qr ready, in both the data generating process and the fiteval. Should RQSS or any other package handle the weights - required by the Expectation Maximization - properly, only the rq call in MS.fitquantreg needs to be adjusted. The next release will use natural cubic basis splines with penalization & cross validation on the number of knots to create a designmatrix that is fed to rq in MS.fitquantreg, enabling non-linear quantile fitting in a workaround fashion. This will also require a Cross Validation wrapper function, that could run in parallel hypothetically.
Further note, that rq also sometimes runs into singularity issues - as does RQSS inevitably - so an extension to recover data even in these cases will be favourable. Particularly if the testing scenarios become exhaustive. In these cases, it becomes fairly frequent, that at least one rq call suffers from singularity, causing all data to be lost.

Further detail, particulary concerning the possible testing scenarios, can be taken from the ROXY documentation in the functions' headers.

For usage of the provided code and errors exceptions please contact <tim.ruhkopf@outlook.de>.
