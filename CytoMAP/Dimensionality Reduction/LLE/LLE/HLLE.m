function [Y, mse] = HLLE(X,k,d)

%HLLE Runs the standard Hessian LLE implementation of Hessian Eigenmaps
%   https://github.com/gpeyre/matlab-toolboxes/blob/master/toolbox_dimreduc/hlle.m
%   X is the high-dimensional data to be processed
%   k is the number of nearest neighbor points to be used 
%   if k is a scalar, same size used at all points. if k is a vector of
%   length N, neighborhood N(i) will be assigned k(i) nearest neighbors
%   
%   d is the number of dimensions to embed X in
%   
%   Y is the output embedded data
%   mse is the sum (at each neighborhood used) of the eigenvalues(d+2:end)
%   of the local coordinate representation. used for adaptive neighborhood
%   restructuring
%   Example:
%   N=1000; k=12; d=2;
%   tt = (3*pi/2)*(1+2*rand(1,N));  height = 21*rand(1,N);
%   X = [tt.*cos(tt); height; tt.*sin(tt)];
%   [Y, mse] = HLLE(X,k,d);

%   C. Grimes and D. Donoho, March 2003
%   Last Revision: 
%   Version 1.0

%get data size
N = size(X,2);
%check for constant neighborhood size
if max(size(k)) ==1
    kvec = repmat(k,N,1);
elseif max(size(k)) == N
    kvec=k;
else
    error('Neighborhood Vector Size does not match data');
end

disp(['-->Running HLLE for ', num2str(N), ' points']);
disp('-->Computing HLLE neighbors');

%Compute Nearest neighbors
D1 = L2_distance(X,X,1);

dim = size(X,1);
nind = repmat(0, size(D1,1), size(D1,2));
%extra term count for quadratic form
dp = d*(d+1)/2;
W = repmat(0,dp*N,N);

if(mean(k)>d) 
  disp('[note: k>d; regularization will be used]'); 
  tol=1e-3; % regularlizer in case constrained fits are ill conditioned
else
  tol=0;
end;

for i=1:N
    tmp = D1(:,i);
    [ts, or] = sort(tmp);
%take k nearest neighbors
    nind(or(2:kvec(i)+1),i) = 1;
    thisx = X(:,or(2:kvec(i)+1));
    %center using the mean 
    thisx = thisx - repmat(mean(thisx')',1,kvec(i));

    %compute local coordinates
    [U,D,Vpr] = svd(thisx);
    V = Vpr(:,1:d);
    
    %Neighborhood diagnostics
    vals = diag(D);
    mse(i) = sum(vals(d+1:end));
    
    
%build Hessian estimator
    clear Yi; clear Pii;
    ct = 0;
    for mm=1:d
        startp = V(:,mm);
        for nn=1:length(mm:d)
            indles = mm:d;
            Yi(:,ct+nn) = startp.*(V(:,indles(nn)));
        end
        ct = ct+length(mm:d);
    end
    Yi = [repmat(1,kvec(i),1), V, Yi];
%orthogonalize linear and quadratic forms
    [Yt, Orig] = mgs(Yi);
    Pii = Yt(:,d+2:end)';
 %double check weights sum to 1
    for j=1:dp
        if sum(Pii(j,:)) >0.0001
            tpp = Pii(j,:)./sum(Pii(j,:)); 
        else
            tpp = Pii(j,:);
        end
        %fill weight matrix
       W((i-1)*dp+j, or(2:kvec(i)+1)) = tpp;
    end
end

%%%%%%%%%%%%%%%%%%%%Compute eigenanalysis of W
disp('-->Computing HLLE embedding');

G=W'*W;
G = sparse(G);

options.disp = 0; 
options.isreal = 1; 
options.issym = 1;

%tol=1e-3; %sometimes useful for pathological fits
tol=0;
% [Yo,eigenvals] = eigs(G,d+1,tol,options); % old code
[Yo,eigenvals] = eigs(G,d+1,'SM',options);
Y = Yo(:,1:d)'*sqrt(N); % bottom evect is [1,1,1,1...] with eval 0


disp('-->Orienting Coordinates');
%compute final coordinate alignment
R = Y'*Y;
R2 = R^(-1/2);
Y = Y*R2;