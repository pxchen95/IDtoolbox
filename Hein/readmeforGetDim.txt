Intrinsic Dimensionality Estimation
===================================

This code estimates the intrinsic dimension of a dataset. It calculates three estimates. The estimator proposed by the authors in

M. Hein, J-Y. Audibert. Intrinsic dimensionality estimation of submanifolds in Euclidean space. Proceedings of the 22nd ICML, 289-296, Eds. L. de Raedt and S. Wrobel, 2005.

as well as two classical estimators: the correlation dimension and the Takens estimator (see the above paper for references).

Two possible implementations are available. The first one is C-Code and the other one is a mex-file for the usage in MATLAB (see https://www.ml.uni-saarland.de/code/IntDim/IntDim.htm).

The usage of this code is allowed for scientific purposes. If you use this code, please cite the above paper.

Matlab Mex-Files
================

GetDim:
This is a mex-file for Matlab. Go to the directory where you have copied the downloaded files. Then use the command
     mex GetDim.cpp
to compile the file.

Let X be the matrix of the data points, where one column corresponds to one data point. Then use
     GetDim(X)
to estimate the intrinsic dimension. The matrix can be sparse or full.

Output are three numbers:
The estimate for the intrinsic dimension, the correlation dimension, the Takens estimator.

IMPORTANT NOTE: To save time the maximal number of dimension is limited to 500 for the intrinsic dimensionality estimator.