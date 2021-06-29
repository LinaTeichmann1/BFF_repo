
%%% FIRST TRY RUNNING IT THROUGH R

%% Compare Bayes Factors with different width priors for colour data
%   Figure 3 contains Bayes Factors for three different prior widths
%   Uses Bayes Factor wrapper to calculate Bayes Factors (needs to be installed in R)   

%   April 2021
%% setup & load
clearvars;
% add toolboxes
addpath(genpath('./local_functions'));

% load data
load('../data_colour/ds_stacked_realcolour.mat')
load('../data_colour/significant_timepoints_realcolour_new.mat')

% colormaps
tv = ds_stacked_realcolour.a.fdim.values{1}*1000;
load('color_bwr.mat');
load('color_gray.mat');

% calculate Bayes Factors
X = ds_stacked_realcolour.samples';
prior_widths = {'ultrawide';'wide';'medium'};
bfs = cell(1,size(prior_widths,1));
tic; 
for i = 1:size(bfs,2)
    bfs{i} =bayesfactor_R_wrapper(X,'args',['mu=0.5,rscale="' prior_widths{i} '",nullInterval=c(0.5,Inf)'],'returnindex',1);
end
toc;

% a single bayes factor
tic;
rscale = [sqrt(2), 1, 1 / sqrt(2)];
bfs_a = cell(1,size(prior_widths,1));
bfs_a_cmp = cell(1,size(prior_widths,1));
tic; 
for i = 1:size(bfs,2)
    [bfs_a{i}, bfs_a_cmp{i}]  = calc_bf_vect(X, 0.5, [0.5, Inf], rscale(i));
end
toc;

nullvalue = 0.5
n = size(X, 2);
mu = mean(X - nullvalue, 2);
sd = std(X, 0, 2);
se = sd / sqrt(n);
t = mu ./ se;

% compare relative and absolute difference
% return 1 if one of them is within tol range
rel_diff = @(a, b, tol) any([(abs(a - b) / min(abs(a), abs(b))) <= tol, ...
            abs(a - b) <= tol]);


%{

Here's a list of the one's that don't match within a 2 dp abs/rel
tolerance.

The reason why they don't match is because the numbers are being 
truncated before they're being sent to R. So the Matlab version is
closer to what you'd get in R if you didn't truncate the numbers.

This is happens for mid-szied *t* values
When the *t* values are very big there's also a deviance
but this is then very small in relative terms. 

%}

t(arrayfun(@(a,b) rel_diff(a, b, 1e-2), bfs{1}, bfs_a{1}) == 0)
t(arrayfun(@(a,b) rel_diff(a, b, 1e-2), bfs{2}, bfs_a{2}) == 0)
t(arrayfun(@(a,b) rel_diff(a, b, 1e-2), bfs{3}, bfs_a{3}) == 0)

bfs{1}(arrayfun(@(a,b) rel_diff(a, b, 1e-2), bfs{1}, bfs_a{1}) == 0)
bfs_a{1}(arrayfun(@(a,b) rel_diff(a, b, 1e-2), bfs{1}, bfs_a{1}) == 0)

t(arrayfun(@(a,b) rel_diff(a, b, 1e-2), bfs{2}, bfs_a{2}) == 0)
t(arrayfun(@(a,b) rel_diff(a, b, 1e-2), bfs{3}, bfs_a{3}) == 0)

