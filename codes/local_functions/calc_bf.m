% Usage:
%
% X is an array of values
%
% nullvalue is the what ever you want the null mu to be (e.g., 0.5)
%
% nullinterval is an array with the lower and upper bounds (e.g., [0.5, Inf]
%
% rscale is the cauchy scale factor
% e.g., 
% medium = sqrt(2) / 2
% wide = 1
% ultrawide = sqrt(2)

function bf = calc_bf(X, nullvalue, nullinterval, rscale)




    if nargin == 1
        nullvalue = 0.5;
        nullinterval = [0.5, Inf];
        rscale = sqrt(2); % 'ultrawide'
    end
    
    % Take the input and convert it to a T stat
    X_mean = mean(X) - nullvalue;
    X_sd = std(X);
    X_n = length(X);
    X_df = X_n - 1;
    
    X_se = X_sd / sqrt(X_n);
    X_t = X_mean / X_se;
    
    % define the likelihood function
    % which is a noncentral t distribution
    likelihood_fun = @(theta) nctpdf(X_t, X_df, theta);
    
    % define the alternative prior
    % which is a cauchy prior with the following parameters
    % location (x_0) = 0
    % scale (gamma) = rscale * sqrt(n)
    % and it's truncated with a lower bound at
    % nullinterval_lowerbound * sqrt(n)
    % note: A cauchy prior is just a student t distribution with df = 1
    
    
    % scale the uppper and lower bound
        
    lowerbound = nullinterval(1) * sqrt(X_n);
    upperbound = nullinterval(2);
 
    alternative_prior = cauchypdf(rscale * sqrt(X_n),...
                                  lowerbound,...
                                  upperbound);
    
    
    % define the null prior
    % which is just a point at 0
    M0 = likelihood_fun(0);
 

    
    alt_model = @(theta) likelihood_fun(theta) .* alternative_prior(theta);
    % this is a workaround because the Matlab algorithm for gaussian quadrate
    % is seemingly a bit shit
    A_BIG_NUMBER = 100; 
    M1 = integral(alt_model,lowerbound,A_BIG_NUMBER);
    
    bf = M1 / M0;
end

function prior_normalised = cauchypdf(scale, lowerbound, upperbound)
        

% NOTE: There is a step here where I have to normalise the prior
% This happens for each BF that you calculate...
% But because it's normalised the same way each time it would be more
% efficient to compute the normalisation factor once and then pass that 
% value around, because the normalisation step is costly. 
% there's a couple of ways to do this
% one way might be write out a global the first time it's called
% or you could work it out at another point and pass it into the BF
% function

prior_full = @(x) tpdf(x / scale, 1) / scale;

prior_trunc = @(x) ifelse(x < lowerbound, 0, prior_full(x));

K = integral(prior_trunc, lowerbound, upperbound);

prior_normalised = @(x) prior_trunc(x) / K;


end

function output = ifelse(condition, yes, no)

if condition
    output = yes;
else
    output = no;
end

        
end