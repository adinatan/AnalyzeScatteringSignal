
function baseline = arPLS_baseline(signal, smoothness_param, min_diff)
% implements Baseline correction using adaptive iteratively reweighted penalized least squares (10.1039/B922045C)

% defaults
if nargin<2
    smoothness_param = 1e3;
    min_diff = 1e-6;
end
signal = signal(:);

order = 2; % Difference filter order
N = 100; % Maximum # of iterations

L = numel(signal);
difference_matrix = diff(speye(L), order);
minimization_matrix = (smoothness_param*difference_matrix')*difference_matrix;
penalty_vector = ones(L,1);

for count = 1:N
    penalty_matrix = spdiags(penalty_vector, 0, L, L);
    % Cholesky decomposition
    C = chol(penalty_matrix + minimization_matrix);
    baseline = C \ (C'\(penalty_vector.*signal));
    d = signal - baseline;
    % make d-, and get penalty_vector^t with mean and std
    dn = d(d<0);
    penalty_vector_temp = 1./(1+exp(2*(d-(2*std(dn)- mean(dn)))/std(dn)));
    % check exit condition and backup
    if norm(penalty_vector-penalty_vector_temp)/norm(penalty_vector) < min_diff
%         count
        break;
    end
    penalty_vector = penalty_vector_temp;
end
