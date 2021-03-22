clearvars;
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

%% plot
tv = ds_stacked_realcolour.a.fdim.values{1}*1000;
num_cols = length(tv); 
co= getPyPlot_cMap('bwr_r',num_cols);
title_names = flipud([{'160 trials (10%)'};{'400 trials (25%)'};{'800 trials (50%)'};{'1200 trials (75%)'};{'1600 trials (100%)'}]);

for p = 1:length(n_participants)
    figure(p);clf;
    for i  = 1:5
        a=subplot(5,1,i);hold on
        exponential_minmax=10;
        val_col_map = logspace(-exponential_minmax,exponential_minmax,num_cols);
        plot(tv(find(sig{p}{i})),repmat(10^9,1,length(find(sig{p}{i}))),'k*'); hold on

        for t = 1:length(tv)
            [~,idx] = min(abs(val_col_map-bfs{p}{i}(t)));  
            st = stem(a,tv(t),bfs{p}{i}(t),'Clipping','off','basevalue',1,'Color','k','MarkerFaceColor',co(idx,1:3),'MarkerEdgeColor','k','LineWidth',1,'MarkerSize',6);
            hold on;
        end
        a.YScale = 'log';
        a.XLim = [tv(1),tv(end)];
        a.YLim = [10^-exponential_minmax, 10^exponential_minmax];
        a.YTick = [10^(-exponential_minmax), 10^(-exponential_minmax/2),10^0, 10^(exponential_minmax/2),10^exponential_minmax];
        xlabel('time (ms)')
        ylabel('BF (log scale)')
        a.FontSize = 14;
        title([title_names{i} ', ' num2str(n_participants(p)) ' participants'])
        
        ax = axes();
        ax.Position = a.Position;
        ax.Position(4) = a.Position(4)/4;
        mu = mean(all_res_subsampled{p}{i}.samples*100);
        se = std(all_res_subsampled{p}{i}.samples*100)./sqrt(size(all_res_subsampled{p}{i}.samples,1));
        fill([tv fliplr(tv)],[mu-se fliplr(mu+se)],[102,0,204]/255,...
            'FaceAlpha',.2,'LineStyle','none');hold on
        plot(tv,mu,'k-','LineWidth',2,'Color',[102,0,204]/255,'Marker','None');
        plot(tv,tv*0+50,'k--')


        ax.YColor=[102,0,204]/255;
        ax.YAxisLocation='right';
        ax.YLabel.String = 'Accuracy (%)';
        ax.Color = 'none';
        ax.Box = 'off';
        ax.XColor='none';
        ax.YLim = [48,60];
    end
    
    f = gcf();
    f.Position = [-1500 230 1070 1310];
    drawnow()
    saveas(f,['../figures/subsampling_n' num2str(n_participants(p)) '.png'])


end


%% make cooler plots
% all_res_subsampled:
%   first 

all_res_subsampled{1}



