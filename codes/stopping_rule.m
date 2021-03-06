%% sample from the results to reflect different sample sizes 
%% setup 
clearvars;rng(1);
addpath(genpath('/System/Volumes/Data/misc/data16/teichmanna2/Toolboxes/CoSMoMVPA'))
addpath(genpath('./local_functions'))
% loading the subsampling results
load('../data_colour/less_trials_results.mat','res10_all','res25_all','res50_all','res75_all')
load('../data_colour/ds_stacked_realcolour.mat')

% defining different sample sizes
n_participants = 5:18;

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
%   bayes factors: half-cauchy, point-null, range = [0.5,Inf], rscale = 0.707
%   (Bayes Factor R package)
bfs=[];
for p = 1:length(n_participants)
    for i  = 1:length(all_res)
        % Bayes Factors
        X = all_res_subsampled{p}{i}.samples';
        bfs{p}{i} = bayesfactor_R_wrapper(X,'args','mu=0.5,rscale="medium",nullInterval=c(0.5,Inf)','returnindex',1);
    end
end

bf_all = [];
for s = 1:size(bfs,2)
    bf_all(:,:,s)=[bfs{s}{:}]';
end

save('../data_colour/bfs_stoppingrule.mat','bf_all')

%%

cols=(magma(5));

toplot=[];
criterion = 6;
for t = 1:size(bf_all,1)
    y = bf_all(t,:,:);
    toplot(:,t)=squeeze(mean(y<1/criterion | y>criterion,2));
end

figure(1);clf;hold on
plot(n_participants,n_participants*0+80,'k--','LineWidth',1)



col_idx = [1,3,5];
col_idx=1:5
for i = 1:length(col_idx)
    a(i)=plot(n_participants,toplot(:,i)*100,'Color',[cols(col_idx(i),:),0.8],'LineWidth',2);
end

xlim([n_participants(1),n_participants(end)]);
ylim([50,100]);

% l=legend(a,{'1600 trials','800 trials','160 trials'});
l=legend(a,{'1600 trials','1200 trials','800 trials','400 trials','160 trials'},'Location','northwest');

ax= gca();

ax.YLabel.String = 'Timepoints fulfilling BF criterion (%)';
ax.XLabel.String = '# of participants';

ax.FontSize = 14;
ax.Box='off';

text(5.1,79,'pre-defined cut-off','FontSize',12)