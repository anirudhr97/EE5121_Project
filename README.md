# Interleaved ISTA-CG (iiCG)

Project for EE5121 Optimisation (July-Nov 2021 Semester) based on :
Solntsev, Stefan, Jorge Nocedal, and Richard H. Byrd. "An algorithm for quadratic ℓ1-regularized optimization with a flexible active-set strategy." Optimization Methods and Software 30, no. 6 (2015): 1213-1237.
[doi link](https://doi.org/10.1080/10556788.2015.1028062)

The original implementation of the algorithm is available [here](https://github.com/stefanks/Ql1-Algorithm), which we built over for our simulations. [alg_ql1.m](./alg_ql1.m) is the main algorithm, for which [Auxiliary](./Auxiliary) contains helper functions.

## Implementations of ISTA and FISTA
[FISTA.m](./FISTA.m) and [ISTA.m](./ISTA.m) contain our own implementations of FISTA and ISTA.
Check out a few lectures of this brilliant series of [lectures](https://youtu.be/JRerBpNggN0) by Prof. Constantine Caramanis (UT Austin) to help you understand the algorithm and theory behind it.

## 2D visualisations
[visualise.m](./visualise.m) is a script that creates 2D solution path visualisations for iiCG, ISTA and FISTA for the same starting point, to compare their behaviour and get a feel for how they work.

## Performance analysis
### Spectra Problem
The spectra dataset can be loaded using [spectra.mat](./spectra.mat).

-- If you do not have the spectra dataset on your MATLAB version, use the 'spectra.mat' file given in this repository. Use the following command 
to import the dataset into your MATLAB workspace:
`spectra = importdata('spectra.mat')`

-- Use the following command to get an idea on what the dataset is about:
`spectra.Description`

Once this dataset has been imported into your MATLAB workspace, the plots shown in the slides can be obtained by running [test_spectra_data.m](./test_spectra_data.m).
The plots generated will be saved as .png files in [plots](./plots).

### Random Matrix Problem
We test the working of the algorithms on randomly generated matrices(The random matrices are generated from the normal distribution).

The plots shown in the slides can be obtained by running [test_random_data.m](./test_random_data.m). The plots generated will be saved as .png files in [plots](./plots).
