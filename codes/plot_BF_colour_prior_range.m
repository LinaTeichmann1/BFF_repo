%% Compare Bayes Factors with different prior ranges (null intervals) for colour data
%   Figure 5 contains Bayes Factors for four different prior ranges
%   Uses Bayes Factor wrapper to calculate Bayes Factors (needs to be installed in R)   

%   April 2021
%% setup & load
clearvars;
% add toolboxes
addpath(genpath('./local_functions'));

% load data
load('../data_colour/ds_stacked_realcolour.mat')
load('../data_colour/significant_timepoints_realcolour_new.mat')

% colormaps
tv = ds_stacked_realcolour.a.fdim.values{1}*1000;
load('color_bwr.mat');
load('color_gray.mat');

% calculate Bayes Factors
X = ds_stacked_realcolour.samples';
null_ints = [0,0.2,0.5,0.8];
bfs = cell(1,size(null_ints,2));
for i = 1:size(bfs,2)
    bfs{i} =bayesfactor_R_wrapper(X,'args',['mu=0.5,rscale="medium",nullInterval=c(' num2str(null_ints(i)) ',Inf)'],'returnindex',1);
end

%% make the plot comparing different prior ranges (null intervals)
f=figure(1);clf
f.Position=[f.Position(1:2) 1000 400];f.Resize='off';f.PaperPositionMode='auto';f.Color='w';
top = 300; bot = 60;
a=axes('Units','Pixels','Position',[55 bot 700 top]);hold on
markers=['-d';'-s';'-o';'-^'];

null_ints_titles={'[0, Inf]','[0.2, Inf]','[0.5, Inf]','[0.8, Inf]'};
for i = 1:size(bfs,2)
    stem(a,tv,bfs{i},markers(i,:),'basevalue',1,'Color','k','LineWidth',0.01,'MarkerSize',8)
end
% plot markers on top of the stems
for i = 1:size(bfs,2)
    plot(tv,bfs{i},markers(i,:),'MarkerFaceColor',co_gray(i,1:3),'MarkerEdgeColor',co_gray(i,1:3),'LineStyle','None','MarkerSize',8);
end

set(a,'YScale','log','XLim',[tv(1),tv(end)], ...
        'YLim',[1e-10 1e10],'YTick',10.^(-10:5:10))
    
a.XLabel.String=('time (ms)');
a.YLabel.String=('BF (log)');
a.FontSize = 20;

% add the prior plots to the right
size_tiny_plots = (top-bot)/4+8;
oneway=1;
limits=5;
null_line_width=1;
x_ax_on=0;
prior_width = 1/sqrt(2);
for i = 1:4
    ax=axes('Units','Pixels','Position',[780 top-(i-1)*size_tiny_plots-9*i size_tiny_plots*2.5 size_tiny_plots]);
    null_int_start=null_ints(i);
    if i == 4
        x_ax_on=1;
    end
    plot_cauchy(ax,null_int_start,oneway,limits,null_line_width,co_bwr,x_ax_on,prior_width)
    hold on 
    plot(ax.XLim(1)+1,ax.YLim(2)-0.1,markers(i,:),'MarkerFaceColor',co_gray(i,1:3),'MarkerEdgeColor',co_gray(i,1:3),'MarkerSize',10)
    text(ax.XLim(1)+1.5,ax.YLim(2)-0.1,null_ints_titles{i},'FontSize',16)
end
set(gcf, 'InvertHardCopy', 'off'); 

%% save plot
set(gcf, 'InvertHardCopy', 'off'); 
fn = sprintf('../figures/figure5');
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im = imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=1;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

