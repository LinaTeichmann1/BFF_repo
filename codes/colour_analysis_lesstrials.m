%% subsampling from the data
% THIS IS THE RESULT: all_res_subsampled, sig & bfs (saved as subsampling_results.mat):
%   first index is for different number of participants (6, 9, 12, 15, 18) 
%   second index is for different number of trials (100%, 75%, 50%, 25%, 10%)

clearvars;rng(1);
addpath(genpath('/System/Volumes/Data/misc/data16/teichmanna2/Toolboxes/CoSMoMVPA'))

participants = [1:6,8:10,12:20];
res75 = {};res50 = {};res25 = {};res10 = {};
for p = participants
    load(['../data_colour/cosmofiles/P' num2str(p) '_ds'])
    dsa = ds_nopca;

    % colour-decoding for shapes == real colour decoding
    dsa = cosmo_slice(dsa,dsa.sa.trialtype==1);
    dsa.sa.targets = dsa.sa.colourcategory;

    n_stims = length(find(dsa.sa.stimulus==1));
    n_trials = [n_stims*0.75,n_stims*0.5,n_stims*0.25,n_stims*0.1];

    for ii = 1:size(n_trials,2)
        stim_i = [];
        for i = 1:max(dsa.sa.stimulus)
            stim_i= [stim_i;datasample(find(dsa.sa.stimulus==i),n_trials(ii),'Replace',false)];
        end

        ds = cosmo_slice(dsa,stim_i);

        ds.sa.chunks=1+mod(ceil(ds.sa.stimulus/2)-1,5);
        partitions = cosmo_nchoosek_partitioner(ds,1);
        partitions = cosmo_balance_partitions(partitions,ds);

        measure=@cosmo_crossvalidation_measure;
        measure_args=struct(); 
        measure_args.classifier=@cosmo_classify_lda;
        measure_args.partitions=partitions;
        measure_args.nproc=4;
        nbrhood=cosmo_interval_neighborhood(ds,'time','radius',0);
        
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

save('less_trials_results.mat','res10_all','res25_all','res50_all','res75_all')


%% stats
load('less_trials_results.mat','res10_all','res25_all','res50_all','res75_all')
load('../data_colour/ds_stacked_realcolour.mat')

n_participants = [6,9,12,15,18];

all_res = {ds_stacked_realcolour,res75_all,res50_all,res25_all,res10_all};

all_res_subsampled = [];
for p = 1:length(n_participants)
    for i  = 1:length(all_res)
        all_res_subsampled{p}{i} = all_res{i};
        all_res_subsampled{p}{i} = cosmo_slice(all_res_subsampled{p}{i},1:n_participants(p));
        all_res_subsampled{p}{i}.sa.targets = ones(size(all_res_subsampled{p}{i}.samples,1),1);
        all_res_subsampled{p}{i}.sa.chunks=(1:size(all_res_subsampled{p}{i}.samples,1))';
    end
end

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

        % BFs
        X = all_res_subsampled{p}{i}.samples';
        bfs{p}{i} = bayesfactor_R_wrapper(X,'args','mu=0.5,rscale="medium",nullInterval=c(0.5,Inf)','returnindex',1);
    end
end

save('subsampling_results.mat','all_res_subsampled','sig','bfs')




