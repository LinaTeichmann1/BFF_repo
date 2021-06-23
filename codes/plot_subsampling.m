%% Compare results for different data sizes
%   Figure contains Bayes Factors for different sample sizes and trial
%   numbers

%   April 2021
%% Setup
load('../data_colour/subsampling_results.mat')
load('../data_colour/ds_stacked_realcolour.mat')

[v_max,i_max]=max(mean(ds_stacked_realcolour.samples))

addpath(genpath('./local_functions'))
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
plottitles=flip(prct)*n_trials;
plottitles=arrayfun( @(x) [num2str(x) ' trials'],plottitles,'UniformOutput',false);

%% plot data
f=figure(1);clf
f.Position = [f.Position(1:2) 800 800];
cols=(magma(5));
for s = 1:size(toplot,3)
ax=subplot(4,1,s);

for i = 1:size(toplot,1)
    y = toplot(i,:,s);
    p_idx = logical(p_toplot(i,:,s));
    a(i)=area(tv,y,'FaceColor',cols(i,:),'FaceAlpha',0.6,'EdgeColor',cols(i,:)); hold on
    ylim([-20,20])
    
    
    plot(tv(p_idx),-9-i*1.5-1,'.k','Color',cols(i,:),'MarkerSize',10)
    plot(tv,tv*0,'k','LineWidth',1)
 
end
    ax.YTick = [-10,0,10,20];
    ax.YLabel.String = 'BF (log)';
    ax.YLabel.Position(2) = 10;
    ax.YLabel.Position(1) = -150;
    ax.FontSize = 14;
    ax.Box='off';
    ax.Title.String = ['n=' num2str(participants(s))];
    ax.Title.Position(2) = 15;
    ax.Title.Position(1) = -50;

    t = text(0,-13,'p<0.05','FontSize',12);
    l=legend(a,{'1600 trials','1200 trials','800 trials','400 trials','160 trials'});
    
    if s == 3
        ax.XLabel.String = 'time (ms)';
    end
end

cols = viridis(3);
toplot2 = toplot(:,i_max,:);
toplot2 = squeeze(toplot2);
ax=axes('Units','normalized','Position',[0.3,0.03,0.6,0.2]);

marker_sizes = (100*5):-100:5;
a = [];
order = [3,2,1];
for i = 1:3
        a(order(i))=scatter([1,2,3,4,5],toplot2(:,i),100,'MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerFaceAlpha',1);hold on
    plot([1,2,3,4,5],toplot2(:,i),'Color',[cols(i,:),0.1],'LineWidth',2)
end
set(ax,'XLim',[0.5,5.5],'XTick',[1,2,3,4,5],'XTickLabels',{'1600 trials','1200 trials','800 trials','400 trials','160 trials'},'FontSize',14,'YLim',[0,18],'XDir','reverse')
ylabel('BF (log)')
l=legend(a,{'n = 18','n = 12','n = 6'});
l.Location='northeastoutside';

title('Peak decoding (125ms)')

annotation('textbox',[.025,.91,.4,.04],'String','A','FontSize',20,'FontWeight','bold','LineStyle','none')
annotation('textbox',[.25,.22,.4,.04],'String','B','FontSize',20,'FontWeight','bold','LineStyle','none')



% save
fn = sprintf('../figures/figure_subsampling');
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=1;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

 