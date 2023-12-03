# Uncertain DFA

<a href="https://www.mozilla.com](https://www.overleaf.com/project/6560a7b71a2908eaada8eb67"> Overleaf paper </a>

Detrended Fluctuation Analysis (DFA) in the perspective of uncertainties regarding among others the interpretation of the alpha coefficient.

According to wikipedia (and the scientific litterature), the interpretation of the alpha coefficient (output of the DFA applied to any signal):
- $\alpha < 1/2$: anti-correlated
- $\alpha \approx 1/2$: uncorrelated, white noise
- $\alpha > 1/2$: correlated
- $\alpha \approx 1$: 1/f-noise, pink noise
- $\alpha > 1$: non-stationary, unbounded
- $\alpha \approx 3/2$: Brownian noise

The first objective of this work is to extend this crisp interpretation partition into a fuzzy partitioning (e.g. a single alpha value could correspond to different types of noise). To do so, tests are carried out based on the simulation of different type of noise (eventually start from a simple sinusoid signal) and on the computation of alpha coefficient associated with a standard DFA. The known classification of the signal should help to define an empirical (i.e. based on simulations) uncertainty model of the fuzzy partition.

The second objective of this work is to make a sensibility analysis on the DFA in order to asses the impact of different parameters (number of points, duration, segments definition, type/order of regression model, etc.). To do so, the standard DFA is being recoded and experiments will be planned to measure the impact of each parameters (multiple effects will be considered as well). 

The third/main objective of this work is to identify all the uncertainties of the DFA process in order to model them and to extend the DFA to an "uncertain" or "imprecise" version. To do so the computed impact of the different parameters could be used, e.g. m(\alpha)=1-1/n where n is the number of considered points. 

The targeted journal is "Applied Soft Computing".

 ### Questions:
 - the first experiment on a sinusoid noised signal show some model mixture in the regression plot (simple sinusoid -> 2 linear models could be mixed) -> could the number of models be extracted from the spectral analysis of the plot (simple sinusoid -> 2 frequential peaks -> 2 linear models?)?
 - what is the real meaning of the alpha-DFA? is it a real complexity level?
 - what is the best complexity indicators among entropy and DFA based known methods?

### TO DO LIST:
#### Implementation:
- Francis: simulation of many noise signal -> computation of the corresponding alpha coefficients (for fuzzyfying the interpretation partitioning the alpha values)
- Nicolas: manually code the whole DFA process -> parametrisation for sensitivity analysis
#### Thinking:
- Francis: try to establish a fuzzy interpretation partitioning the alpha values (from simulations)
- Nicolas: if models mixture is used (2 models -> 2 slopes -> 2 fdc/probas to be fused? -> use of the E2M algorithm? what type of combination?
Reading:
- both: bibliography on DFA, uncertain/imprecise regression

## UPDATE (24/11/2023)

### Research trails
- extend the geometrical growth of fluctuations depending on a complexity level (F=n^{alpha}) to other growing models
- imprecision modelling (fuzzy-possibility-belief) in case of low quality final DFA regression -> imprecise alpha
- Detrended Distributionnal Analysis (DDA): instead of focusing on fluctuations (ie error between the signal and its linear trend), focus on the distribution based on histogram. The underlying assumption = some (natural) pattern might be at different times of the signal and at different scales.

### Open questions
- why a gaussian noise could have memory?
- what is the type(s) of signal the DFA can handle (not only in terms of size)?

### TODOLIST
- sensitivity analysis of initial window choice (or on initial pattern length?
