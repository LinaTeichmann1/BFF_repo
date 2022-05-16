%% Compare Bayes Factors with different width priors for colour data
%   Figure 3 contains Bayes Factors for three different prior widths
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
%prior_widths = {2;1;0.707};
prior_widths = {sqrt(2); 1; sqrt(2)/2};
prior_widthstxt = {'ultrawide';'wide';'medium'};
bfs_o = cell(1,size(prior_widths,1));
bfs = cell(1,size(prior_widths,1));
tic;
for i = 1:size(bfs,2)
  bfs_o{i} =bayesfactor_R_wrapper(X,'args',['mu=0.5,rscale="' prior_widthstxt{i} '",nullInterval=c(0.5,Inf)'],'returnindex',1);
end
toc
tic;
for i = 1:size(bfs,2)
  m =  mean(X,2) - 0.5;
  s = std(X',0);
  d = m' ./ s;
  n = repmat(size(X,2),1, size(X,1));
  location = repmat(0,1,size(X,1));
  scale = repmat(prior_widths{i},1,size(X,1));
  ll = repmat(0.5,1,size(X,1));
  ul = repmat(Inf,1,size(X,1));
  bfs{i} = BayesFactor(d', n', location', scale', ll', ul');
end
toc
%% make the plot comparing different prior widths
f=figure(1);clf
f.Position=[f.Position(1:2) 1000 400];f.Resize='off';f.PaperPositionMode='auto';f.Color='w';
top = 300; bot = 60;
a=axes('Units','Pixels','Position',[85 bot 700 top]); hold on
markers=['-d';'-s';'-o';'-^'];
null_ints_titles={'ultrawide (r = 1.414)','wide (r = 1)','medium (r = 0.707)'};

for i = 1:size(bfs,2)
    stem(a,tv,bfs{i},markers(i,:),'basevalue',1,'Color','k','LineWidth',0.01,'MarkerSize',8)
    plot(tv,bfs{i},markers(i,:),'MarkerFaceColor',co_gray(i+1,1:3),'MarkerEdgeColor',co_gray(i,1:3),'LineStyle','None','MarkerSize',8);
end

set(a,'YScale','log','XLim',[tv(1),tv(end)], ...
        'YLim',[1e-10 1e10],'YTick',10.^(-10:5:10))
a.XLabel.String=('time (ms)');
a.YLabel.String=('BF (log scale)');
a.FontSize = 20;

% add the prior plots
size_tiny_plots = (top-bot)/3+8;
oneway=1;
limits=5;
null_line_width=1;
x_ax_on=0;
null_int_start=0.5;
widths = [sqrt(2),1,1/sqrt(2)];
for i = 1:size(bfs,2)
    ax=axes('Units','Pixels','Position',[820 top-(i-1)*size_tiny_plots-18*i size_tiny_plots*2 size_tiny_plots ]);    
    prior_width=widths(i);
    if i == 3
        x_ax_on=1;
    end
    plot_cauchy(ax,null_int_start,oneway,limits,null_line_width,co_bwr,x_ax_on,prior_width)
    hold on 
    plot(ax.XLim(1)+1,ax.YLim(2)-0.1,markers(i,:),'MarkerFaceColor',co_gray(i+1,1:3),'MarkerEdgeColor',co_gray(i,1:3),'MarkerSize',10)
    text(ax.XLim(1)+1.5,ax.YLim(2)-0.1,null_ints_titles{i},'FontSize',16)
end


%% save plot
set(gcf, 'InvertHardCopy', 'off'); 
fn = sprintf('../figures/figure3');
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im = imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=1;
%imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');



