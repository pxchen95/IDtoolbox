Intrinsic Dimensionality Toolbox (MATLAB)
Author: Paula Chen, Division of Applied Mathematics, Brown University
Last Updated: September 10th, 2019

KEY:

N = sample size
D = extrinsic/embedding dimension
d = intrinsic dimension
X = DxN input data matrix (each column = 1 data pt)

NOTE: .cpp/.c files must be compiled, i.e. use command: mex subfolder\___.cpp 
      (currently contains a total of 4 mex files)
      Algorithms (i.e. neighbor-based methods) using pdist2 may crash for N > ~10^6.


=====================================================================================
DATA GENERATION (Folder: Generate_Data)
=====================================================================================

GenerateManifoldData.cpp: mex Generate_Data\GenerateManifoldData.cpp
                          X = GenerateManifoldData(#,D,N)

Adapted from: Hein, Matthias, and Jean-Yves Audibert. "Intrinsic dimensionality estimation of submanifolds in R^d." Proceedings of the 22nd international conference on Machine learning. ACM, 2005.

Source link: https://www.ml.uni-saarland.de/code/IntDim/IntDim.htm

Generates X of manifold type #

#  | d   | D         | Manifold                                           
-------------------------------------------------------------------------------------
0    1     3 (fixed)   Sinusoid over circle w/ high curvature:                 
                          t\in [0,2pi) --> (sin t, cos t, 1/10 sin(150 t))
1    D-1   D           Uniformly sampled sphere                     
2    3     5 (fixed)   Affine space                                 
3    4     6 (fixed)   Strange figure (highly concentrated so ~3d)  
4    4     8 (fixed)   Manifold                                     
5    2     3 (fixed)   Helix                                        
6    6     36 (fixed)  Manifold                                     
7    2     3  (fixed)  Swiss roll                                   
8    12    72 (fixed)  Manifold, effectively lies in a 24D subspace
                          x(t): [0,1]^12 --> R^72
                          x_{2i-1}(t) = t_{i+1} cos(2pi t_i), i = 1, ..., 11
                            x_{2i}(t) = t_{i+1} sin(2pi t_i), i = 1, ..., 11
                            x_{23}(t) = t_1 cos(2pi t_{12}
                            x_{24}(t) = t_1 sin(2pi t_{12}
                          x_{j+24}(t) = x_{j+48} = x_j,       j = 1, ..., 24
9    20    20 (fixed)  Affine space                                 
10   D-1   D           Hypercube (padded with 0s)                   
11   2     3 (fixed)   Moebius band, twisted 10 times           
                          (u,v)\in [-1,1]x[0,2pi) --> ((1+u/2 cos(k/2 v)) cos v,
                                                       (1+u/2 cos(k/2 v)) sin v,
                                                       u/2 sin(k/2 v))
12   D     D           Multivariate Gaussian, i.e. normally sample from D-plane
13   1     D           Curve

-------------------------------------------------------------------------------------

gen_plane.m: X = gen_plane(d,D,N,pad_zero)

Generates X uniformly sampled from a d-dimensional plane, embedded in D-dimensions either by padding with zeros (pad_zero = TRUE) or via a random distance-preserving linear transform (pad_zero = FALSE)

-------------------------------------------------------------------------------------

transform.m: X = transform(X,D)

Embeds X into D dimensions via a random distance-preserving linear transform


=====================================================================================
ID ESTIMATION METHODS
=====================================================================================

Folder: PCA

dim_PCA.m: d = dim_PCA(X,T)

Performs global PCA on X, s.t. d dimensions account for (1-T)*100% of the variance of the data.

Notes:
      + efficient & fast
      - linear method (fit d-dimensional hyperplane to entire dataset)

-------------------------------------------------------------------------------------

Folder: localPCA

Note that this method is not meant to be used to estimate ID, although it theoretically could be used to do so (in the no noise case, the estimate would be the d for which the T vs. d plot hits the x-axis; in the noise case, we would be looking for some type of elbow point in this plot). In developing this algorithm we were hoping to answer the following questions:
      (1) How linear is a dataset?
      (2) If a different algorithm estimates the ID as d, can we use the results
          from local PCA to extrapolate how much (variance) of the data is explained 
          by d dimensions? --- I'm not sure this is possible since the ID does not 
          emcompass the shape of the manifold approximated by the algorithm, so 
          a piecewise linear estimate (i.e. the result of local PCA) is not
          necessarily comparable to however the other algorithms may fit the data to a 
          manifold.

local_pca.m: d = local_pca(X,k,n)
     k = # of neighbors (including the center point)
     n = # of center points to consider

Do PCA on n neighborhoods (n center pts chosen uniformly randomly from input data, without replacement) of size k (including the center pt). Projection error is calculated on each of the neighbors in the neighborhood (incl. the center pt).

local_pca_2k.m: d = local_pca_2k(X,k,n)
     k = # of points to perform PCA on (incl. center pt)/# of points to fit
     n = # of center points to consider

Pick n neighborhoods (center pts picked uniformly randomly w/o replacement) of size 2k. Do PCA on every other neighbor (incl. the center pt). Projection error is calculated using the other k neighbors.

Currently, the projection error is the (sum of distance b/t neighbor and its projection onto d dimensions)/(sum of variance of neighborhood to global data mean). Output (T) is the average over all center points considered. See localPCA\readmeforlocalpca for more details on the projection error.

-------------------------------------------------------------------------------------

Folder: kNN

nearneighbor.m: d = nearneighbor(X,k,tol,maxiter)

Implementation of basic near neighbor algorithm as described in: Pettis, Karl W., et al. "An intrinsic dimensionality estimator from near-neighbor information." IEEE Transactions on pattern analysis and machine intelligence 1 (1979): 25-37.

     k           = neighborhood size (# of neighbors)
     tol/maxiter = tolerance/max # of iterations for the iterative method for         
                   calculating the ID estimate

Derived from the expected value of the average distance to the kth nearest neighbor.

Notes: 
      + iterative step usu. requires few iterations (usu. maxiter = 4 suffices)
      + consistent for k = 2
      - tends to underestimate (d > ~10) - why? (this occurs ind. of k)
      - slow with large data sets (N > 10^5)
      - difficult to analyze even with simple underlying distribs

-------------------------------------------------------------------------------------

Folder: 2NN

twoNN.m: d = twoNN(X)

Implementation of two nearest neighbor algorithm as described in: Facco, Elena, et al. "Estimating the intrinsic dimension of datasets by a minimal neighborhood information." Scientific reports 7.1 (2017): 12140.

Derived from the cdf of the ratio between the distance to the 2nd nearest neighbor and the distance to the nearest neighbor.

Notes:
      + (paper claims) using few neighbors should avoid effects of curvature/density
        variation
      - involves arbitrary outlier cutoff (already coded in)
      - requires assumption of local uniformity

-------------------------------------------------------------------------------------

Folder: MLE

mledim.m: d = mledim(X,k1,k2)

Maximum Likelihood Estimator algorithm as described in: Levina, Elizaveta, and Peter J. Bickel. "Maximum likelihood estimation of intrinsic dimension." Advances in neural information processing systems. 2005.

Source link: http://www.stat.lsa.umich.edu/~elevina/mledim.m

Derived by considering the log-likelihood of the observed process N(t) (= # of observations within distance t from any of the observed pts). N(t) is approximated as a Poisson process. Estimate d is the MLE associated w/ the described log-likelihood function. Author replaces epsilon-neighborhoods with k-neighborhoods. Average over different k values to try to avoid effects of neighborhood size.

Notes: 
      + existence and uniqueness of MLE asymptotically guaranteed (pf uses 
        exponential families)
      | as with all MLEs, tradeoff b/t bias and variance
      - unclear how to choose parameters k1, k2 (MLE is k-dependent, author uses k1
        = 10, k2 = 20)
      - performed worse than all other algorithms tested for ~all test manifolds
      - requires very large N for high d (d >~ 20)

-------------------------------------------------------------------------------------

Folder: DANCo

DANCo.m: d = DANCo(X)

Dimensionality from Angle and Norm Concentration algorithm as described in: Ceruti, Claudio, et al. "DANCo: Dimensionality from angle and norm concentration." arXiv preprint arXiv:1206.3881 (2012).

Source link: https://www.mathworks.com/matlabcentral/fileexchange/40112-intrinsic-dimensionality-estimation-techniques

Estimator minimizes the sum of the Kullback-Leibler divergences (relative entropies) of the distribution of normalized nearest-neighbor distances and the distribution of pairwise angles. Using KL divergence is meant to avoid edge effects. Using both of these particular KL divergences is meant to reduce tendency to underestimate when using just nearest-neighbor distances alone.

Notes:
      + a faster variant exists (FastDANCo, see source link)
      - empirically seems less consistent than other methods on toy manifolds (high 
        variance b/t trials)
      - should theoretically be able to provide a non-integer value (i.e. fractal dim
        estimate), but this is implemented using fminsearch, which can get stuck on   
        local minima, which leads to wildly off ID estimates

-------------------------------------------------------------------------------------

Folder: Hein

GetDim.cpp: mex Hein\GetDim.cpp
            d = GetDim(X); 
                d(1) = intrinsic dimension estimator (as described by Hein, et al), 
                       smoothed version of CD estimator, replaces step function w/ 
                       kernel (1-x)_+ (chosen for comp. efficiency, authors claim 
                       choice of kernel does not change result much)
                d(2) = correlation dimension (CD) estimator, CD is a type of fractal 
                       dim that describes amt of space a random set of pts from the 
                       manifold will take up; I recommend using this instead of 
		       intrinsic_dim(X,'CorrDim')
                d(3) = Takens estimator, alt formulation of CD

ID estimator as described in: Hein, Matthias, and Jean-Yves Audibert. "Intrinsic dimensionality estimation of submanifolds in R^d." Proceedings of the 22nd international conference on Machine learning. ACM, 2005.

Source link: https://www.ml.uni-saarland.de/code/IntDim/IntDim.htm

Note that the correlation dimension and Takens estimates can also be considered as ID estimates.

CD source: Grassberger, Peter, and Itamar Procaccia. "Measuring the strangeness of strange attractors." Physica D: Nonlinear Phenomena 9.1-2 (1983): 189-208.

Takens estimate source: Takens, Floris. "On the numerical determination of the dimension of an attractor." Dynamical systems and bifurcations. Springer, Berlin, Heidelberg, 1985. 99-106.

Notes: 
  Hein:
      - highly concentrated probability measure leads to underestimation
      - high curvature leads to overestimation
  CD:
      - reported to require large amount of data for accurate estimates
      - severely underestimates if data is nonuniformly distributed

-------------------------------------------------------------------------------------

Folder: Other

intrinsic_dim.m: d = intrinsic_dim(X,method)

Dimension Reduction toolbox. It is unclear exactly which sources he was following in writing the algorithms. Below, I cite my best guesses at which algorithm I think he was following. The full toolbox also includes dimension reduction and some denoising techniques.

     method = 'EigValue'
              Global PCA, should be identical to PCA\dim_PCA.m

     method = 'MLE'
              Presumably the same as MLE\mledim.m?

     method = 'CorrDim' 
              Correlation Dimension (see Folder: Hein, for more details), presumably
              the same as output d(2) of Hein\GetDim.cpp?
              Implementation doesn't use the least squares fitting that all sources
              describe...

     method = 'NearNbDim'
              k-Nearest Neighbor algorithm. Some sort of fitting over a range of
              neighborhood sizes.

     method = 'PackingNumbers'
              Approximation of the Capacity Dimension, which does not depend on the
              data distribution on the manifold

     Source? Kégl, Balázs. "Intrinsic dimension estimation using packing numbers." 
             Advances in neural information processing systems. 2003.

     method = 'GMST'
              Geodesic Minimum Spanning Tree - connect all vertices while minimizing
              total edge weight, which are weighted via an approximation of geodesic 
              distances - global version of kNNG?

     Source? Costa, Jose A., and Alfred O. Hero. "Geodesic entropic graphs for
             dimension and entropy estimation in manifold learning." IEEE Transactions 
             on Signal Processing 52.8 (2004): 2210-2221.
      
Source link: https://lvdmaaten.github.io/drtoolbox/

Notes:
      - I have NOT verified any of his code - the above sources may not be accurate,
        as it was unclear exactly which algorithm/paper he was following
      - method = 'GMST' is very slow for N >= 10^5

-------------------------------------------------------------------------------------

Folder: ISOMAP

IsoMap.m: IsoMap(X,n_fcn,n_size)
     n_fcn  = use 'epsilon' or 'k' neighborhoods
     n_size = neighborhood size (radius or # of neighbors)

Method as described in: Tenenbaum, Joshua B. "Mapping a manifold of perceptual observations." Advances in neural information processing systems. 1998.

Source link: https://www.mathworks.com/matlabcentral/fileexchange/62449-isomap-d-n_fcn-n_size-options

Use geodesic instead of Euclidean distances b/t pts to define neighborhoods. ID estimate is defined to be dimension at which the "elbow" of the residual variance vs. ISOMAP dimensionality graph occurs (1st output figure).

Notes:
      - not yet automated to generate ID estimate, must visually find "elbow pt"

-------------------------------------------------------------------------------------

Folder: kNNG
     
        mex kNNG\kNNlengthmex.c
        mex kNNG\kNNgraphmex.c

knn_graph_estim_1.m: [d,~,~,~] = knn_graph_estim_1(X,k,gamma,M,N,samp)
knn_graph_estim_2.m: [d,~,~,~] = knn_graph_estim_2(X,k,gamma,M,N,samp) 
                     (supposedly more efficient for high D, o/w should be identical)
     k     = # of neighbors
     gamma = edge weighting parameter
     M     = # of independent least-square runs (?)
     N     = # of resampling trials per LS dwell (?)
     samp  = range of # of sample sizes to use for boostrap estimate of mean length
             fxn

swiss_roll_example.m: demo of above two algorithms using a swiss roll and a sphere.

Algorithm as described in: J. A. Costa and A. O Hero, "Entropic Graphs for Manifold Learning", Proc. of IEEE Asilomar Conf. on Signals, Systems, and Computers, Pacific Groove, CA, November, 2003.

Source link: http://web.eecs.umich.edu/~hero/IntrinsicDim/

Creates kNN graph, where edges b/t points are weighted by some distance function (approx geodesic distances?).

Notes:
      | for below reasons, I don't recommend using this algorithm...
      - VERY slow even for reasonable sample sizes (i.e. N = 2000)
      - not sure what parameters M,N are for...
      - unclear how input parameters should be picked
      - results for toy manifolds were VERY off for some manifolds - somehow got 
        negative values for some?


=====================================================================================
SCRIPTS FOR RUNNING/COMPARING ALGORITHMS (Folder: Scripts)
=====================================================================================

compare.m

Runs algorithms on toy manifold problems. Plots results either as a scatter plot of sample size vs. ID estimate or as a bar plot of difference between true and estimated ID.

-------------------------------------------------------------------------------------

run_local_pca.m

Runs and generates plots for local_pca.m and local_pca_2k.m. Adds old results from compare.m onto a plot of T vs. d for local PCA using toy manifold data. Also plots T vs. d for local PCA and 2k local PCA, as run on toy manifolds.


=====================================================================================
PAST RESULTS (Folder: Sample_Results)
=====================================================================================

compare_results.mat

Results from running compare.m. Note that in compare_results.mat, the CD estimate is calculated using intrinsic_dim(X,'CorrDim') (Other/intrinsic_dim.m), although I have now modified compare.m to use [~, d, ~] = GetDim(X) (Hein/GetDim.m) as the CD estimate.

-------------------------------------------------------------------------------------

Folder: compare_algorithms

Plots generated using compare.m (ID estimate vs. sample size -- each figure represents a different manifold; results resummarized in a bar plot). The data used for these plots is stored in compare_results.mat.

Notes:
      * the nonlinear algorithms all underestimate the ID of linear manifolds
      * since in practice, ID estimates are rounded to nearest integer, any estimate 
        in [d-0.5,d+0.5] would be considered "correct"

-------------------------------------------------------------------------------------

Folder: run_local_pca

Plots for different toy manifolds of T vs. d for local PCA and for 2k local PCA using different neighborhood sizes. Input data used N = 5000, n = N, and no noise.

Notes: 
      * for local PCA, the resulting ID estimate d <= k, so for k too small, local PCA 
        will never be able to get the correct estimate (will maximally estimate d = k)
      * for k = N, local PCA = global PCA
      * note that for the 12d manifold in 72d, local PCA appears to predict 24d, which 
        is consistent with how the data on this manifold is generated; I postulate 
        that if we sample densely enough, i.e. increase N enough, then theoretically 
        local PCA should be able to estimate d correctly


=====================================================================================
OTHER NOTES
=====================================================================================

PUBLICLY AVAILABLE KINEMATIC DATA SETS:

(1) Note that a large publicly available data set (Dryad Data Package) is available here: https://doi.org/10.5061/dryad.1k84r

Source (Original Publication): Atzori M, Gijsberts A, Castellini C, Caputo B, Mittaz Hager A, Elsig S, Giatsidis G, Bassetto F, Müller H (2014) Electromyography data for non-invasive naturally controlled robotic hand prostheses. Scientific Data 1:140053. https://doi.org/10.1038/sdata.2014.53

Source (Dryad Data Package): Atzori M, Gijsberts A, Castellini C, Caputo B, Mittaz Hager A, Elsig S, Giatsidis G, Bassetto F, Müller H (2014) Data from: Electromyography data for non-invasive naturally controlled robotic hand prostheses. Dryad Digital Repository. https://doi.org/10.5061/dryad.1k84r


(2) Another data set (Ninapro Dataset 4) is also available here: https://doi.org/10.5281/zenodo.1000138

(3) and (Ninapro Dataset 5) here: https://doi.org/10.5281/zenodo.1000116

Source (Original Publication): Pizzolato et al., Comparison of Six Electromyography Acquisition Setups on Hand Movement Classification Tasks, Plos One 2017 (accepted).