function [g_dist] = geodesic_dist(Q1, Q2)


% Q1 = rand(4, 4);
% Q2 = rand(4,4);

Q = Q1^(-0.5) * Q2 * Q1^(-0.5);
% Q = sqrt( inv(Q1) ) * Q2 * sqrt( inv(Q1) );
% th = 0.001;
% 
% [U, S, V] = svd(Q1);
% for i = 1: length(S)
%     if S(i, i) < th
%        S(i,i) = th;
%     end
% end
% 
% Q1_tilde = U * S^(-0.5) * V';
% 
% Q = Q1_tilde * Q2 * Q1_tilde;


[~, d] = eig(Q);
log_d = log2(diag(d));
g_dist = sqrt(log_d' * log_d );

end