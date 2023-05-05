function [t,rndt,zsc] = riemann_dist_rnd(trueNetwork,estNetwork,flagsin)

% Compute the Riemannian distance between two non-negative definite matrices.
% Matrices are assumd to be symmetric.
%
% FC matrix can have either all ones diagonal or all zeros.
% For SC matrices, the symmetric normalized Laplacian is first computed for both networks.
%
% Input:
% trueNetwork is the "True" network,
% estNetwork is the estimated network -- the input order matters.
% Need to specify in flagsin the type of network, flagsin.FC=1 for FC, or flagsin.SC=1 for SC.
%
% Functions called:
% squareform_sp
% LaplacianMtx
%
% Outputs:
% t: Riemann distance between trueNetwork and estNetwork, scalar
% rndt: Vector of Riemann distances between estNetowrk and 1000 randomized
% versions of trueNetwork
% zsc: struct,
% zsc.zvalue: z-score of obtained Riemann distances
% zsc.zmean: z-mean
% zsc.zsigma: z-sigma

%% Main program

A = trueNetwork;
B = estNetwork;

if(flagsin.SC == flagsin.FC)
    error('Error, flags FC and SC cannot be both 1 or 0');
end

szrnd = 1000;

t = riemann_dist(A,B,flagsin);

rndt = zeros(szrnd,1);

% Here randomize "true" network A
for ii=1:szrnd
    U = rnd_graph(A);
    rndt(ii) = riemann_dist(U,B,flagsin);
end

% Compute zscore:
[zsc.zvalue,zsc.zmean,zsc.zsigma] = zscore([rndt' t]);


%% Compute Riemann distance 

function t = riemann_dist(A,B,flagsin)

epsilon = 0.000001;

% Make sure diagonal elements are one for FC so as to avoid negative eigenvalues 

if(flagsin.FC == 1)
    if(trace(A) == 0)
        A = A + eye(size(A,1));
    end

    if(trace(B) == 0)
        B = B + eye(size(B,1));
    end
    [Va,Da] = eig(A);
end

% Compute Riemann distance for SC networks
if(flagsin.SC == 1)
    
    A = LaplacianMtx(A);
    B = LaplacianMtx(B);

    [Va,Da] = eig(A);
end

Da = diag(Da); % Vectorize eigenvalues 
idx = (Da > epsilon); % Exclude eigenvalues below threshold epsilon
diagDa = 1./sqrt( Da(idx) );
diagDa = diag(diagDa);

A = Va(:,idx) * diagDa * Va(:,idx)';

M = real(A * B * A); % Since matrix A is real and symmetric then its eigenvalues are real
M = (M + M')/2; % M is symmetric

clear A B Va Da diagDa idx;

[~,Dm] = eig(M);
Dm = real(Dm);
Dm = diag(Dm); % Extract diagonal elements
idx = (Dm > epsilon); % Exclude eigenvalues below threshold epsilon

t = sqrt( sum(log(Dm(idx)).^2) );

clear idx

%% Generate random networks

function U = rnd_graph(A)

% All randomized networks preserve symmetry.
% Zero edges in the original network are preserved. This is useful mainly
% for the SC networks.
%
% Required functions:
% squareform_sp

U = triu(A,1);
u = vect_upper_tri(U); % Vectorize uppper triangular part of the matrix

idx = (u~=0);
idx2 = randperm(length(idx));
idx3 = idx2(idx);
u(idx) = u( idx3 );
%U = squareform_sp(u); % Rebuild u into a matrix
U = squareform(u); % Diagonal elements are zero
U = U + eye(size(U,1));

function vA2 = vect_upper_tri(A)

% Vectorize an upper triangular matrix

At = A.';
m  = (1:size(At,1)).' > (1:size(At,2));
vA2  = At(m);
clear m At;
