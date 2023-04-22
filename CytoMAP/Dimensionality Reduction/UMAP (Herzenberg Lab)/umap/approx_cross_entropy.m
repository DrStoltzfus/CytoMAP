function result = approx_cross_entropy(dists, weights, a, b)
%APPROX_CROSS_ENTROPY Given a distance for each 1-simplex in low-dimensional space
% and the original weights of the 1-simplices in high-dimensional space,
% compute the approximation to the cross-entropy between the two simplicial
% complexes. This calculation uses the modified smooth formula Phi for
% low-dimensional weight that is used in the stochastic gradient descent.
%
% result = APPROX_CROSS_ENTROPY(dists, weights, a, b)
%
% Parameters
% ----------
% dists: array of size (n_1_simplices, 1)
%     The current distance between the two endpoints of the 1-simplex in
%     low-dimensional Euclidean space.
%
% weights: array of size (n_1_simplices, 1)
%     The original weights assigned to the 1-simplices in high-dimensional
%     space.
%
% a: double
%     Parameter of differentiable approximation of right adjoint functor.
% 
% b: double
%     Parameter of differentiable approximation of right adjoint functor.
% 
% Returns
% -------
% result: double
%     The total approximated cross entropy between the two simplicial complexes.
%
% See also: CROSS_ENTROPY
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <cgmeehan@alumni.caltech.edu>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

Phi = ones(size(dists,1),1)./(1 + a*(dists.^(2*b)));
Phi_summands = -(weights.*log(Phi) + (1-weights).*log(1-Phi));

Phi_summands(isnan(Phi_summands)) = 0;
result = sum(Phi_summands);
end