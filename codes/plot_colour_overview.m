%% Load colour decoding results & make an overview plot (Figure 2)
%   Figure 2 contains decoding accuracy, effect size, p-values & Bayes Factors
%   Uses Bayes Factor wrapper to calculate Bayes Factors (needs to be installed in R)   

%   April 2021
%% setup & load
clearvars
% add toolboxes
addpath(genpath('./local_functions'));

% load data
load('../data_colour/ds_stacked_realcolour.mat')
load('../data_colour/significant_timepoints_realcolour_new.mat')

% colormaps
tv = ds_stacked_realcolour.a.fdim.values{1}*1000;
load('color_bwr.mat');

% calculate Bayes Factors
X = ds_stacked_realcolour.samples';
bf = bayesfactor_R_wrapper(X,'args','mu=0.5,rscale="medium",nullInterval=c(0.5,Inf)','returnindex',1);

%% make the plot 
f=figure(1);clf
f.Position = [f.Position(1:2) 900 700];

% SUBPLOT 1: decoding accuracy
a=subplot(3,1,1);hold on
mu = mean(ds_stacked_realcolour.samples*100);
se = std(ds_stacked_realcolour.samples*100)./sqrt(size(ds_stacked_realcolour.samples,1));
f=fill([tv fliplr(tv)],[mu-se fliplr(mu+se)],'k',...
    'FaceAlpha',.2,'LineStyle','none');
plot(tv,mu,'LineWidth',2,'Color','k');
a.XLim = tv([1 end]);
a.YLim = [47,60];
plot(a,tv(1:15:end),zeros(length(tv(1:15:end)))+50,'k','LineStyle','--')
xlabel('time (ms)')
ylabel({'Decoding accuracy (%)'})
a.FontSize = 14;

% SUBPLOT 2: effect size with p-values
a=subplot(3,1,2);hold on
t_stat=(mean(ds_stacked_realcolour.samples)-0.5)./(std(ds_stacked_realcolour.samples));
plot(a,tv,t_stat,'k','LineWidth',2)
plot(a,tv(significant_timepoints_realcolour),-1,'k.','MarkerSize',15,'Color',[0.5,0.5,0.5])
a.XLim = tv([1 end]);
a.YLim = [-2,4];
xlabel('time (ms)')
a.YLabel.String = ["effect size","(Cohen's d)"];
plot(a,tv(1:15:end),zeros(length(tv(1:15:end))),'k','LineStyle','--')
a.FontSize = 14;

% SUBPLOT 3: stemp plot for bayes
a=subplot(3,1,3);hold on
exponential_minmax=10;
val_col_map = logspace(-exponential_minmax,exponential_minmax,size(co_bwr,1));
for t = 1:length(tv)
    [~,idx] = min(abs(val_col_map-bf(t)));  
    st = stem(a,tv(t),bf(t),'Clipping','off','basevalue',1,'Color','k','MarkerFaceColor',co_bwr(idx,1:3),'MarkerEdgeColor','k','LineWidth',1,'MarkerSize',6);
    hold on;

end
set(a,'YScale','log','XLim',[tv(1),tv(end)], ...
        'YLim',[1e-10 1e10],'YTick',10.^(-10:5:10))
xlabel('time (ms)')
ylabel('BF (log scale)')
a.FontSize = 14;

% colorbar
colormap(co_bwr)
cbh = colorbar;
caxis([-exponential_minmax,exponential_minmax])
cbh.Units = 'normalized';
cbh.Limits = [-exponential_minmax,exponential_minmax];
cbh.Position(1) = 0.92;cbh.Position(3) = 0.01;cbh.Position(4) = a.Position(4);cbh.Position(2) = a.Position(2);
cbh.Label.String = 'Bayes Factor';
cbh.TickLabels=arrayfun(@(x) ['10^{' num2str(x) '}'], cbh.Ticks, 'UniformOutput', false);

% add the prior plot
axis_size = 0.08;
offset_box = 0.01;
ax=axes('Units','normalized','Position',[a.Position(1)+a.Position(3)-axis_size-offset_box,a.Position(2)+a.Position(4)-axis_size-offset_box, axis_size, axis_size]);
null_int_start=0.5;
oneway=1;
limits=5;
null_line_width=1.5;
x_ax_on=0;
prior_width = 1/sqrt(2);
plot_cauchy(ax,null_int_start,oneway,limits,null_line_width,co_bwr,x_ax_on,prior_width)

%% save plot
fn = sprintf('../figures/figure2');
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im = imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=1;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');



