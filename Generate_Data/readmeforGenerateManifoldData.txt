Intrinsic Dimensionality Estimation
===================================

This code generates datasets (and more) used in

M. Hein, J-Y. Audibert. Intrinsic dimensionality estimation of submanifolds in Euclidean space. Proceedings of the 22nd ICML, 289-296, Eds. L. de Raedt and S. Wrobel, 2005.

Two possible implementations are available. The first one is C-Code and the other one is a mex-file for the usage in MATLAB (see https://www.ml.uni-saarland.de/code/IntDim/IntDim.htm).

The usage of this code is allowed for scientific purposes. If you use this code, please cite the above paper.


Matlab Mex-Files
================

GenerateManifoldData:
This is a mex-file for Matlab. Go to the directory where you have copied the downloaded files. Then use the command

mex generateManifoldData

generateManifoldData has three inputs
1) The number of the dataset 
2) The dimension of the dataset/extrinsic dimension (D)
3) The number of points (N)

e.g. generateManifoldData(0,200,5000)
generates 5000 data points of dataset 0 of dimension 3 (the dimension 200 is ignored since the dimension is fixed for this dataset)

Here is a description of the possible datasets:

0:  1d sinusoid over the circle with high curvature in 3d space (dimension fixed)
1:  a D-1 dimensional sphere in D space (uniformly sampled) 
2:  a 3d affine space in 5d space                               (dimension fixed)
3:  a strange 4d figure in 6d (but very concentrated so it is ess 3d) (dim fixed)
4:  a 4d manifold in 8d space                                   (dimension fixed)
5:  a 2d helix in 3d space                                      (dimension fixed)
6:  a 6-dim manifold in 36d space                               (dimension fixed)
7:  2d "swiss roll" in 3d space                                 (dimension fixed)
8:  12d manifold in 72d space                                   (dimension fixed)
9:  a 20d affine space in 20d space                             (dimension fixed)
10: D-1 dimension hypercube, uniformly sampled in D space (padded with 0s)
11: 2d Moebius band, 10 times twisted in 3d space               (dimension fixed)
12: D dimensional Multivariate Gaussian in D space
13: 1d curve in D space