%% Compare results for different data sizes
%   Figure contains Bayes Factors for different sample sizes and trial
%   numbers

%   April 2021
%% Setup
load('../data_colour/subsampling_results.mat')
colors = tab10;
co = colors([3,5,10],:);

prct = fliplr([0.1 .25 .50 .75 1.00]);
n_trials = 1600;
participants=[6,12,18];
tv = 1000*all_res_subsampled{1}{1}.a.fdim.values{1};

% restructure data to plot
toplot = [];p_toplot =[];
for s = 1:5
    toplot(:,:,s)=log([bfs{s}{:}])';
    p_toplot(:,:,s)=vertcat(sig{s}{:});
end

% only plotting 6,12 & 18 participants --> get rid of 9 & 15 participants
toplot(:,:,[2,4])=[];
p_toplot(:,:,[2,4])=[];

% make plot titles
plottitles = prct*n_trials;
plottitles=arrayfun( @(x) [num2str(x) ' trials'],plottitles,'UniformOutput',false);

%% make the plot comparing different data sizes
f=figure(1);clf
f.Position = [f.Position(1:2) 800 800];

% make 5 subplots, each showing Bayes Factors for different # of trials
% used
plotnr=0;
for s = 1:size(toplot,1)
    plotnr=plotnr+1;
    a=subplot(5,1,plotnr); hold on
    % loop over the three different sample sizes
    for i = 1:size(toplot,3)
        y = toplot(s,:,i);
        p_idx = logical(p_toplot(s,:,i));
        bf_idx = logical(toplot(s,:,i)>0);
        
        stem3(tv,i+0*tv,y,'LineWidth',1,'Color',[0.7,0.7,0.7],'LineWidth',1,'Marker','None');
        stem3(tv(bf_idx),i+0*tv(bf_idx),y(bf_idx),'LineWidth',1,'Color',co(i,:),'LineWidth',1,'Marker','None');
        plot3(tv,i+0*tv,y,'k.','MarkerFaceColor',[0.7,0.7,0.7],'MarkerSize',6)
        plot3(tv(bf_idx),i+0*tv(bf_idx),y(bf_idx),'k.','MarkerFaceColor',[0.7,0.7,0.7],'MarkerEdgeColor',co(i,:),'MarkerSize',6)
        plot3(tv(p_idx),i+0*tv(p_idx),tv(p_idx)*0-3,'k^','MarkerSize',2,'MarkerFaceColor',co(i,:),'MarkerEdgeColor',co(i,:))
        plot3(tv,i+0*tv,tv*0,'k-') 
        
        a.View=[-10,30];
    end
    
    % add tick labels 
    a.XLabel.String = 'time (ms)';
    a.ZLabel.String = 'BF';a.ZLabel.Position(3) = 15;a.ZLabel.Position(1) = -180;
    a.ZTick = [0 10 20]; 
    a.ZTickLabel=arrayfun(@(x) ['10^{' num2str(x) '}'], a.ZTick, 'UniformOutput', false);
    a.ZLim = [-10,20];
    a.YTickLabel = participants;
    for i = 1:length(a.YTickLabel)
        ticklabels_new{i} = [['\color[rgb]{' num2str(co(i,:)) '}'],'n = ',num2str(participants(i))];
    end
    a.YTickLabel=ticklabels_new;
    a.FontSize=14;  
    % add title
    a.Title.String=plottitles{s};
    a.Title.Position(1) = tv(1)+50;
        
end

%% save
fn = sprintf('../figures/figure_subsampling');
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=1;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

 

    
    