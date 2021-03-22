function bf = bayesfactor_R_wrapper(data,varargin)
    %% Function that calls the standard ttestBF function in R from R Bayesfactor Package
    % Requires a working installation of R with the Bayesfactor package installed
    % This function works by writing data to a temporary textfile, creating
    % and calling a temporary R script, and reading in its outputs.
    %
    % @ Tijl Grootswagers, 2020
    %
    % Usage:
    %
    % bf = bayesfactor_R_wrapper(data, ...)
    %
    % Inputs:
    %   data            PxR data for P tests and R observations (subjects)
    % Optional inputs:
    %   'Rpath',p       a string with the path to Rscript. Use 'which Rscript'
    %                   in a unix terminal to find this path. This should
    %                   not be the path to R (default '/usr/local/bin/Rscript')
    %   'returnindex',r which BF to return, can be either 1 or 2 (default 1).
    %                   a value of 2 returns the complementary bf (e.g. for
    %                   use with a null-interval)
    %   'verbose',v     if true, print the output from R, recommended to check
    %                   that everything works as intended (default false).
    %   'args',a        a string with arguments to pass to ttestBF in R
    %                   must be R code (default ''). Examples:
    %                   - standard JZS medium-width cauchy prior (same as default):
    %                     'mu=0,rscale="medium"'
    %                   - medium-width half-cauchy with a null-interval on 
    %                     the effect size (use with returnindex 2):
    %                     'mu=0,rscale="medium",nullInterval=c(-Inf,0.5)'
    %                   - two-sided medium-width half-cauchy with a null-interval 
    %                     on the effect size (use with returnindex 2):
    %                     'mu=0,rscale="medium",nullInterval=c(-0.5,0.5)'
    % Returns:
    %   bf              Px1 vector with Bayes factors
    % 
    % For more example uses and interpretation, see ttestBF documentation at
    % https://www.rdocumentation.org/packages/BayesFactor/versions/0.9.12-4.2/topics/ttestBF
    % https://richarddmorey.github.io/BayesFactor/#onesample
    %
    % To cite the R package ‘BayesFactor’ in publications use:
    %   Richard D. Morey and Jeffrey N. Rouder (2018). BayesFactor: Computation of
    %   Bayes Factors for Common Designs. https://CRAN.R-project.org/package=BayesFactor
    %   
    
    %% deal with input arguments
    opt=struct();
    opt.Rpath='/usr/local/bin/Rscript';
    opt.returnindex=1;
    opt.args='';
    opt.verbose=false;
    % read input key-value pairs and overwrite the default values
    fnames = varargin(1:2:end);
    fvalues = varargin(2:2:end);
    assert(numel(fnames)==numel(fvalues),'invalid input: number of keys must match number of values')
    for f=1:numel(fnames)
        assert(isfield(opt,fnames{f}),sprintf('invalid input: %s',fnames{f}))
        opt.(fnames{f}) = fvalues{f};
    end
    
    %% check inputs
    assert(ischar(opt.Rpath),'value for ''Rpath'' must be a character array')
    assert(ischar(opt.args),'value for ''args'' must be a character array')
    assert(ismember(opt.returnindex,[1 2]),'value for ''returnindex'' must be 1 or 2')
    
    %% create some temporary filenames
    bfstatstempdir = tempdir;
    bfstatsinfn = [bfstatstempdir 'in.csv'];
    bfstatsoutfn = [bfstatstempdir 'out.csv'];
    bfscriptfn = [bfstatstempdir 'bf_fun.r'];
        
    % write data to csv file
    %writematrix(data,bfstatsinfn);
    writetable(array2table(data),bfstatsinfn,'WriteVariableNames',0);
    
    %%
    fid = fopen(bfscriptfn,'w');
    fprintf(fid,['library("BayesFactor")\n'...
        'X = as.matrix(read.table("%s",header=FALSE,sep=","))\n'...
        'bf10 = numeric(length = nrow(X))\n'...
        'for (row in 1:nrow(X)) {data = as.numeric(X[row,])\n'...
        'bf = ttestBF(x=data,%s)\n'...
        'print(bf)\n',...
        'bf10[row] = as.vector(bf[%i])}\n'...
        'write.table(bf10,"%s",sep=",",row.names=FALSE,col.names="BF10")\n'],...
        bfstatsinfn,opt.args,opt.returnindex,bfstatsoutfn);
    fclose(fid);

    %% test if calling Rscript
    cmd = sprintf('which %s',opt.Rpath);
    [status,output]=system(cmd);
    if opt.verbose
        disp(output)
    end
    if status
        error('%s \n R not found',output)
    end
    
    %%
    cmd = sprintf('%s %s',opt.Rpath,bfscriptfn);
    [status,output]=system(cmd);
    if opt.verbose
        disp(output)
    end
    if status
        error('%s \n non-zero exit status (%i) in calling R',output,status)
    else
        x = readtable(bfstatsoutfn);
        bf = x.BF10;
    end
    