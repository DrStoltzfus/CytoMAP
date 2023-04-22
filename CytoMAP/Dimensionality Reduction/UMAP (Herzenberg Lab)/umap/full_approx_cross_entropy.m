function result = full_approx_cross_entropy(head_embedding, tail_embedding, head, tail, weights, a, b, same_embedding)
%FULL_APPROX_CROSS_ENTROPY Given a distance for each 1-simplex in low-dimensional space
% and the original weights of the 1-simplices in high-dimensional space,
% compute the approximation to the cross-entropy between the two simplicial
% complexes. This calculation uses the modified smooth formula Phi for
% low-dimensional weight that is used in the stochastic gradient descent.
%
% result = FULL_APPROX_CROSS_ENTROPY(dists, weights, a, b)
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
    n1 = size(head_embedding, 1);
    n2 = size(tail_embedding, 1);

    full_dists = pdist2(head_embedding, tail_embedding);
    full_weights = full(sparse(head, tail, weights, n1, n2));
    if same_embedding
        full_weights = full_weights + eye(n1);
    end
    Phi = ones(size(full_weights))./(1 + a*(full_dists.^(2*b)));
    fw0 = full_weights == 0;
    fw1 = full_weights == 1;
    summands_0 = fw0.*log(1-Phi);
    summands_1 = fw1.*log(Phi);
    summands_2 = (~fw0 & ~fw1).*(full_weights.*log(Phi) + (1-full_weights).*log(1-Phi));
    summands_0(isnan(summands_0)) = 0;
    summands_1(isnan(summands_1)) = 0;
    summands_2(isnan(summands_2)) = 0;
    Phi_summands = -(summands_0 + summands_1 + summands_2);

    result = sum(sum(Phi_summands));
end