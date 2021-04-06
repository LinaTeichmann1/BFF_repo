%% Subsampling from the colour data
%   Subsampling from the colour data & re-running 
%   Output: all_res_subsampled, sig & bfs (saved as subsampling_results.mat):
%       first index is for different number of participants (6, 9, 12, 15, 18) 
%       second index is for different number of trials (100%, 75%, 50%, 25%, 10%)
%   This code requires the Bayes Factor R package & CoSMoMVPA
%   Raw data can be downloaded here: https://osf.io/tcqjh/ -- Chapter3_ImpliedColour -- data_MEG_preprocessed

% April 2021
%% setup 
clearvars;rng(1);
addpath(genpath('/System/Volumes/Data/misc/data16/teichmanna2/Toolboxes/CoSMoMVPA'))

participants = [1:6,8:10,12:20]; % exclude outliers as defined in paper
res75 = {};res50 = {};res25 = {};res10 = {};

%% loop over participant and re-run classification
for p = participants
    % load the cdata containing the preprocessed data for each participant
    load(['../data_colour/cosmofiles/P' num2str(p) '_ds'])
    dsa = ds_nopca;

    % select the relevant trials ("real" colour decoding only)
    dsa = cosmo_slice(dsa,dsa.sa.trialtype==1);
    
    % set targets as colour (classifier is distinguishing between red and green)
    dsa.sa.targets = dsa.sa.colourcategory;

    % calculate # of trials for each subsampling iteration, for each stimulus
    n_stims = length(find(dsa.sa.stimulus==1));
    n_trials = [n_stims*0.75,n_stims*0.5,n_stims*0.25,n_stims*0.1];

    % loop over subsampling iterations (number of trials) 
    for ii = 1:size(n_trials,2)
        stim_i = [];
        % loop over all stimuli and randomly select a certain number of trials
        for i = 1:max(dsa.sa.stimulus)
            stim_i= [stim_i;datasample(find(dsa.sa.stimulus==i),n_trials(ii),'Replace',false)];
        end
        
        % slice dataset to retain the sub-selection only
        ds = cosmo_slice(dsa,stim_i);
        
        % re-run classification
        ds.sa.chunks=1+mod(ceil(ds.sa.stimulus/2)-1,5);
        partitions = cosmo_nchoosek_partitioner(ds,1);
        partitions = cosmo_balance_partitions(partitions,ds);

        measure=@cosmo_crossvalidation_measure;
        measure_args=struct(); 
        measure_args.classifier=@cosmo_classify_lda;
        measure_args.partitions=partitions;
        measure_args.nproc=4;
        nbrhood=cosmo_interval_neighborhood(ds,'time','radius',0);
        
        % save classification results for each subsampling iteration
        switch ii
            case 1
                res75{end+1}=cosmo_searchlight(ds,nbrhood,measure,measure_args);
            case 2
                res50{end+1}=cosmo_searchlight(ds,nbrhood,measure,measure_args);
            case 3
                res25{end+1}=cosmo_searchlight(ds,nbrhood,measure,measure_args);
            case 4
                res10{end+1}=cosmo_searchlight(ds,nbrhood,measure,measure_args);
        end
    end

end

% stack & save
res10_all = cosmo_stack(res10);
res25_all = cosmo_stack(res25);
res50_all = cosmo_stack(res50);
res75_all = cosmo_stack(res75);

save('../data_colour/less_trials_results.mat','res10_all','res25_all','res50_all','res75_all')


%% sample from the results to reflect different sample sizes 
% loading the subsampling results
load('../data_colour/less_trials_results.mat','res10_all','res25_all','res50_all','res75_all')
load('../data_colour/ds_stacked_realcolour.mat')

% defining different sample sizes
n_participants = [6,9,12,15,18];

% combine all results
all_res = {ds_stacked_realcolour,res75_all,res50_all,res25_all,res10_all};

% loop over the data and select data from subgroup of participants only
all_res_subsampled = [];
for p = 1:length(n_participants)
    for i  = 1:length(all_res)
        all_res_subsampled{p}{i} = all_res{i};
        all_res_subsampled{p}{i} = cosmo_slice(all_res_subsampled{p}{i},1:n_participants(p));
        all_res_subsampled{p}{i}.sa.targets = ones(size(all_res_subsampled{p}{i}.samples,1),1);
        all_res_subsampled{p}{i}.sa.chunks=(1:size(all_res_subsampled{p}{i}.samples,1))';
    end
end

%% Statistics
%   frequentist: sign-flip permutations, TFCE, monte-carlo cluster
%   correction (CoSMoMVPA)
%   bayes factors: half-cauchy, point-null, range = [0.5,Inf], rscale = 0.707
%   (Bayes Factor R package)
sig = [];bfs=[];
for p = 1:length(n_participants)
    for i  = 1:length(all_res)
        % permutations & p values
        neighborhood = cosmo_cluster_neighborhood(all_res_subsampled{p}{i}); 
        opt = struct();
        opt.niter = 10000;
        opt.h0_mean=1/2;
        opt.nproc = 8;

        tfce_ds = cosmo_montecarlo_cluster_stat(all_res_subsampled{p}{i},neighborhood,opt);
        sig{p}{i}=tfce_ds.samples>1.6449; 

        % Bayes Factors
        X = all_res_subsampled{p}{i}.samples';
        bfs{p}{i} = bayesfactor_R_wrapper(X,'args','mu=0.5,rscale="medium",nullInterval=c(0.5,Inf)','returnindex',1);
    end
end

% save
save('../data_colour/subsampling_results.mat','all_res_subsampled','sig','bfs')




