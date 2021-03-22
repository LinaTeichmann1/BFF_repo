%%
% addpath('~/CoSMoMVPA/mvpa/')
% addpath('~/Repository/CommonFunctions/matplotlib/')

%% load data 
clearvars
load('../data_colour/ds_stacked_realcolour.mat')
load('../data_colour/permutations_stacked.mat')

addpath(genpath('./local_functions'))

%%
timevec = ds_stacked_realcolour.a.fdim.values{1};
X = ds_stacked_realcolour.samples' - 0.5;
bf_args = 'mu=0,rscale="medium",nullInterval=c(0.5,Inf)';
bf_interval = bayesfactor_R_wrapper(X,'args',bf_args,'returnindex',1 );

%% run bf for permutations
for p = 1:100
    disp([' running perm ' num2str(p)])
    X_p = perm_null_stacked{p}.samples' - 0.5;
    bf_interval_perm(:,p) = bayesfactor_R_wrapper(X_p,'args',bf_args,'returnindex',1 );
end

save('bf_permutations.mat','bf_interval_perm')


% bootstrapping DO I NEED TO DO THIS?
nboot=10000;rng(1);
rand_idx=datasample(1:100,nboot);
bf_boot = [];
for i = 1:nboot
    bf_boot(:,i) = bf_interval_perm(:,rand_idx(i));  
end



cft = 10;
clusters_observed = bwboundaries(bf_interval>cft);

for b=1:size(bf_boot,2)
    c = bwboundaries(bf_boot(:,b)>cft);
    % maxsize
    if isempty(c)
        clusterstat(b) = 1;
    else
        clusterstat(b) = max(cellfun(@(x) size(x,1), c));
    end
end

clusters_corrected = clusters_observed(cellfun(@(x) size(x,1), clusters_observed) > prctile(clusterstat,95));

f=figure(1);clf
f.Position = [f.Position(1:2) 700 500];
a=subplot(2,2,1);
plot(timevec,mean(X'))
legend('decoding data (colours)')
a.XLim=timevec([1 end]);

a=subplot(2,2,2);
plot(sort(clusterstat))
a.YScale = 'log';
line([.95 .95].*numel(clusterstat),a.YLim,'Color','k')
text(.95.*numel(clusterstat),prctile(clusterstat,95),sprintf('95%% = %i  ',prctile(clusterstat,95)),'HorizontalAlignment','right')
xlabel('permutation (sorted)')
ylabel('maximum cluster size')

a=subplot(2,1,2);
stem(timevec,bf_interval,'k','basevalue',1,'DisplayName','BF');hold on
idx = bf_interval>cft;
plot(timevec(idx),1/20+0*bf_interval(idx),'ko','MarkerFaceColor','k','DisplayName',sprintf('BF>%i (uncorrected)',cft))
a.YScale = 'log';a.XLim=timevec([1 end]);
co = parula
for c=1:numel(clusters_corrected)
    idx = unique(clusters_corrected{c}(:,1)');
    plot(timevec(idx),1./(100*c) + 0*bf_interval(idx),'o','MarkerFaceColor',co(c,:),'Color',co(c,:),'DisplayName',sprintf('BF>%i cluster %i (corrected)',cft,c));
end
legend()


% save
fn = '../figures/figure_cluster_corr';
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');


