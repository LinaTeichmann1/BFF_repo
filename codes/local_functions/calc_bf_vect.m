% Usage:
%
% X             PxR data for P tests and R observations
% nullvalue is the what ever you want the null mu to be (e.g., 0.5)
%
% nullinterval is an array with the lower and upper bounds (e.g., [0.5, Inf]
%
% rscale is the cauchy scale factor
% e.g., 
% medium = sqrt(2) / 2
% wide = 1
% ultrawide = sqrt(2)
%

% some dummy inputs
function [bfs_a, bfs_b] = calc_bf_vect(X, nullvalue, nullinterval, rscale)
% nullvalue = 0.5
% rscale = sqrt(2)
% nullinterval = [0.5, Inf]
% first we'll work out the summary stats
n = size(X, 2);
mu = mean(X - nullvalue, 2);
sd = std(X, 0, 2);
se = sd / sqrt(n);
t = mu ./ se;

LIM = @(T) 10 * sqrt(((2 * n) / ((n) * (n)) + (((T * sqrt(n)) * (T * sqrt(n))) / (2 * (2 * n)))));

likelihood_fun = @(theta, t) nctpdf(t, n - 1, theta);
prior_full = @(x) tpdf(x / (rscale * sqrt(n)), 1) / (rscale * sqrt(n));

M0 = arrayfun(@(T) likelihood_fun(0, T), t);

if  all(isinf(abs(nullinterval)))
    M1 = arrayfun(@(T) integral(@(x) likelihood_fun(x, T) .* prior_full(x), -1 * LIM(T), LIM(T)), t);
    bfs_a = M1_a ./ M0;
    bfs_b = zeros(size(t));
else 

    boundpoint = nullinterval(~isinf(nullinterval)) * sqrt(n);
    lowerbound = boundpoint;
    upperbound = Inf;

    % normalising constant for the BF
    K_a = integral(prior_full, lowerbound, upperbound);

    % normalising contant for the complement
    K_b = 1 - K_a;

    prior_a = @(x) prior_full(x) / K_a;
    prior_b = @(x) prior_full(x) / K_b;

    % you might want to change this into a parallel for loop but just be weary of
    % overheads, because it might not speed it up
    M1_a = arrayfun(@(T) integral(@(x) likelihood_fun(x, T) .* prior_a(x), boundpoint, LIM(T)), t);
    M1_b = arrayfun(@(T) integral(@(x) likelihood_fun(x, T) .* prior_b(x), -1 * LIM(T), boundpoint), t);
    bfs_a = M1_a ./ M0;
    bfs_b = M1_b ./ M0;
end

