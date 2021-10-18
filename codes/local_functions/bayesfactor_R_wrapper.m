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
    %   'Rpath',p       String with the path to Rscript. Use 'which Rscript'
    %                   in a unix terminal to find this path. This should
    %                   not be the path to R (default '/usr/local/bin/Rscript')
    %   'verbose',v     If true, print the output from R, recommended to check
    %                   that everything works as intended (default false).
    %   'args',a        String with arguments to pass to ttestBF in R.
    %                   Must be valid R code (default '').
    %                   Examples:
    %                   - standard JZS medium-width cauchy prior (default):
    %                       'mu=0,rscale="medium"'
    %                   - medium-width one-sided half-cauchy:
    %                       'mu=0,rscale="medium",nullInterval=c(0,Inf)'
    %                   - medium-width one-sided half-cauchy specifying an 
    %                     interval defined on the alternative:
    %                       'mu=0,rscale="medium",nullInterval=c(0.5,Inf)'
    %                   - two-sided medium-width half-cauchy with an interval 
    %                     defined on the alternative (use with returnindex 2):
    %                       'mu=0,rscale="medium",nullInterval=c(-0.5,0.5)'
    %   'tempdir'       A string defining the path to store temporary files.
    %                   This function uses temporary files to pass to R.
    %                   Use this when calling this function in parallel.
    %   'returnindex',r BF to return, can be either 1 or 2 (default 1).
    %                   a value of 2 returns the bf for the complementary 
    %                   interval when specifying an interval hypothesis
    % Returns:
    %   bf              Px1 vector with Bayes factors
    % 
    % Example use of this function:
    % X = randn(1,20) % random data (20 observations)
    % bf = bayesfactor_R_wrapper(X,'args','nullInterval=c(0.5,Inf)','verbose',1)
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
    opt.verbose=true;
    opt.tempdir=tempdir;
    % read input key-value pairs and overwrite the default values
    fnames = varargin(1:2:end);
    fvalues = varargin(2:2:end);
    assert(numel(fnames)==numel(fvalues),'invalid input: number of keys must match number of values')
    for f=1:numel(fnames)
        assert(isfield(opt,fnames{f}),sprintf('invalid input: %s',fnames{f}))
        opt.(fnames{f}) = fvalues{f};
    end
    
    %% check inputs
    assert(ischar(opt.Rpath),'value for ''Rpath'' must be a string')
    assert(ischar(opt.tempdir),'value for ''tempdir'' must be a string')
    assert(ischar(opt.args),'value for ''args'' must be a string')
    assert(ismember(opt.returnindex,[1 2]),'value for ''returnindex'' must be 1 or 2')
    assert(size(data,2)>1,'the input data should have multiple observations')
    assert(all(any(data')),'the input data contains rows that are all zero or nan')
    assert(all(std(data,[],2)),'the input data contains rows with zero variance')
    
    %% test if path to Rscript exists
    if ~isfile(opt.Rpath)
        error('Rscript not found at: %s \n Try setting opt.Rpath',opt.Rpath)
    end
    
    %% create some temporary filenames
    bfstatstempdir = fullfile(opt.tempdir,tempname);
    mkdir(bfstatstempdir)
    bfstatsinfn = [bfstatstempdir 'in.csv'];
    bfstatsoutfn = [bfstatstempdir 'out.csv'];
    bfscriptfn = [bfstatstempdir 'bf_fun.r'];
        
    % write data to temporary csv file
    %writematrix(data,bfstatsinfn);
    writetable(array2table(data),bfstatsinfn,'WriteVariableNames',0);
    
    % write temporary R script that runs the BF with a cauchy prior
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
    
    %% call R and read results from the temporary file that R created
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
    assert(size(bf,1)==size(data,1),'R output did not match the data input size. Something went wrong.')
    