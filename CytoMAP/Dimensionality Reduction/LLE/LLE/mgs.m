function [Q, R] = mgs(A);
%
%   modified Gram-Schmidt. this is a more stable way to compute a
%   qr factorization
%
   [m, n] = size(A);
%
%   we assume that m>= n.
%   
   V = A;
   R = zeros(n,n);
   for i=1:n
      R(i,i) = norm(V(:,i));
      V(:,i) = V(:,i)/R(i,i);
      if (i < n)
         for j = i+1:n
            R(i,j) = V(:,i)' * V(:,j);
            V(:,j) = V(:,j) - R(i,j) * V(:,i);
         end
      end
   end
   Q = V;