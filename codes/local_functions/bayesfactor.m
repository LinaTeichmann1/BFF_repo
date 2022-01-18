function [bf, bf_complement] = bayesfactor(X, varargin)
    % [bf, bf_complement] = bayesfactor(X, varargin)
    % X                 PxR data for P tests and R observations (subjects)
    % Optional inputs:
    %   'nullvalue'     Mean under null hypothesis (default: 0)
    %   'interval'      Interval for alternative hypothesis
    %                   For a full cauchy (two sided), use [-Inf Inf]
    %                   For a half-cauchy (one-sided), use [0 Inf]
    %                   For a half-cauchy + interval, use [0.5 Inf]
    %                   default: [-Inf Inf]
    %   'rscale'        Cauchy scale factor, default: 0.7071
    %                     medium = sqrt(2) / 2
    %                     wide = 1
    %                     ultrawide = sqrt(2)
    % Returns:
    %   bf              Px1 vector with Bayes factors
    %   bf_complement   Px1 vector with Bayes factors for complementary
    %                   interval (in case an interval was specified)
    % Example:
    %   X = randn(5,10);
    %   bayesfactor(X) % standard
    %   bayesfactor(X,'interval',[0.5 Inf]) % with interval on HA

    %% deal with input arguments
    opt=struct();
    opt.verbose=false;
    opt.nullvalue = 0;
    opt.interval = [-Inf Inf];
    opt.rscale = 0.7071;
    % read input key-value pairs and overwrite the default values
    fnames = varargin(1:2:end);
    fvalues = varargin(2:2:end);
    assert(numel(fnames)==numel(fvalues),'invalid input: number of keys must match number of values')
    for f=1:numel(fnames)
        assert(isfield(opt,fnames{f}),sprintf('invalid input: %s',fnames{f}))
        opt.(fnames{f}) = fvalues{f};
    end
    
    % for now, throw an error with a bounded interval
    assert(isinf(opt.interval(2)),'Upper bound of interval must be infinite');
    
    % summary stats
    n = size(X, 2);
    X = X - opt.nullvalue;
    mu = mean(X, 2);
    sd = std(X, 0, 2);
    se = sd / sqrt(n);
    t = mu ./ se;

    LIM = @(T) 10 * sqrt(((2 * n) / ((n) * (n)) + (((T * sqrt(n)) * (T * sqrt(n))) / (2 * (2 * n)))));

    likelihood_fun = @(theta, t) nctpdf(t, n - 1, theta);
    prior_full = @(x) tpdf(x / (opt.rscale * sqrt(n)), 1) / (opt.rscale * sqrt(n));

    % H0
    M0 = likelihood_fun(0, t);

    % HA
    if all(isinf(abs(opt.interval))) % without interval
        M1_a = arrayfun(@(T) integral(@(x) likelihood_fun(x, T) .* prior_full(x), -1 * LIM(T), LIM(T)), t);
        bf = M1_a ./ M0;
        bf_complement = bf;
    else % with interval
        boundpoint = opt.interval(1) * sqrt(n);
        lowerbound = boundpoint;
        upperbound = Inf;
        
        % normalising constant for the BF
        K_a = integral(prior_full, lowerbound, upperbound);

        % normalising contant for the complement
        K_b = 1 - K_a;

        prior_a = @(x) prior_full(x) / K_a;
        prior_b = @(x) prior_full(x) / K_b;

        M1_a = arrayfun(@(T) integral(@(x) likelihood_fun(x, T) .* prior_a(x), boundpoint, LIM(T)), t);
        M1_b = arrayfun(@(T) integral(@(x) likelihood_fun(x, T) .* prior_b(x), -1 * LIM(T), boundpoint), t);

        bf = M1_a ./ M0;
        bf_complement = M1_b ./ M0;
    end
end