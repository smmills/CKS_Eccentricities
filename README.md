# CKS_Eccentricities
Analysis code for determining eccentricities from Kepler transit durations.

This package contains the code associated with the paper Mills et al. 2019 (link here). Please cite this work if you make use of this code or its results. Included are scripts for running the analysis and producing some useful figures, as well as the data used for the analysis.

Each python script includes a short description at the top, and require various numerical python packages to run including:
- numpy
- matplotlib
- scipy

The data used is contained in the `resources` directory, including the stellar data and relevant Kepler data. Please cite the relevant papers if you use this data, e.g., Morton et al. 2016 (https://arxiv.org/abs/1605.02825) and Fulton & Petigura 2018 (https://arxiv.org/abs/1805.01453).

Also included is a text file `koipreflist.txt`, which contains a list of every singly transiting KOI which passed the quality cuts described in the paper, it's log likelihood preference for high eccentricity, measured duration, expected circular edge-on duration, and planet radius.



