%% Makes the different BF figures for the colour data set

%% setup & load
clearvars

% add toolboxes
addpath(genpath('./local_functions'));

% load data
load('../data_colour/ds_stacked_realcolour.mat')
load('../data_colour/significant_timepoints_realcolour_new.mat')


%% colormaps
tv = ds_stacked_realcolour.a.fdim.values{1}*1000;

num_cols = length(tv); 

co_bwr= getPyPlot_cMap('bwr_r',num_cols);



%% ***********************************************************
% ************* Figure 2: colour results *********
% ************************************************************

% make a plot comparing BFs with p-values
% run the BFs
X = ds_stacked_realcolour.samples';
bf_interval = bayesfactor_R_wrapper(X,'args','mu=0.5,rscale="medium",nullInterval=c(0.5,Inf)','returnindex',1);

% setup the plot
tv = ds_stacked_realcolour.a.fdim.values{1}*1000;
f=figure(1);clf
f.Position = [f.Position(1:2) 900 700];
num_cols = length(tv); 
co= getPyPlot_cMap('bwr_r',num_cols);

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

% SUBPLOT 2: cohen's d with p-values
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
val_col_map = logspace(-exponential_minmax,exponential_minmax,num_cols);
for t = 1:length(tv)
    [~,idx] = min(abs(val_col_map-bf_interval(t)));  
    st = stem(a,tv(t),bf_interval(t),'Clipping','off','basevalue',1,'Color','k','MarkerFaceColor',co(idx,1:3),'MarkerEdgeColor','k','LineWidth',1,'MarkerSize',6);
    hold on;

end
set(a,'YScale','log','XLim',[tv(1),tv(end)], ...
        'YLim',[1e-10 1e10],'YTick',10.^(-10:5:10))
xlabel('time (ms)')
ylabel('BF (log scale)')
a.FontSize = 14;
 % colorbar
colormap(co)
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
plot_cauchy(ax,null_int_start,oneway,limits,null_line_width,co,x_ax_on,prior_width)

% SAVE ALL 
saveas(gcf, '../figures/figure2.png')
%%
% ZOOMED IN: prior plot alone
f=figure(1);clf;
f.Position=[f.Position(1:2) 250 200];
ax = axes();
null_line_width = 4;
x_ax_on = 1;
plot_cauchy(ax,null_int_start,oneway,limits,null_line_width,co_bwr,x_ax_on,prior_width)
hold on
saveas(gcf, '../figures/figure_prior_only.png')


%% ***********************************************************
% ************************* Figure 3 *************************
% ************************************************************

% make a plot comparing different prior width
% run the BFs
tv = ds_stacked_realcolour.a.fdim.values{1}*1000;
co= getPyPlot_cMap('gray',5);

X = ds_stacked_realcolour.samples';
prior_widths = {'ultrawide';'wide';'medium'};
bfs_toplot=[];
for i = 1:3
    bfs_toplot{i} =bayesfactor_R_wrapper(X,'args',['mu=0.5,rscale="' prior_widths{i} '",nullInterval=c(0.5,Inf)'],'returnindex',1);
end

f=figure(3);clf
f.Position=[f.Position(1:2) 1000 400];f.Resize='off';f.PaperPositionMode='auto';f.Color='w';
top = 300; bot = 60;
a=axes('Units','Pixels','Position',[55 bot 800 top]);

hold on
markers=['-d';'-s';'-o';'-^'];
null_ints_titles={'ultrawide (r = 1.414)','wide (r = 1)','medium (r = 0.707)'};
for i = 1:3
    stem(a,tv,bfs_toplot{i},markers(i,:),'basevalue',1,'Color','k','LineWidth',0.01,'MarkerSize',8)
    hold on
end

for i = 1:3
    p(i)=plot(tv,bfs_toplot{i},markers(i,:),'MarkerFaceColor',co(i+1,1:3),'MarkerEdgeColor',co(i,1:3),'LineStyle','None','MarkerSize',8)
end
legend(p,null_ints_titles,'Location','NE','FontSize',24,'Box','off')


set(a,'YScale','log','XLim',[tv(1),tv(end)], ...
        'YLim',[1e-10 1e10],'YTick',10.^(-10:5:10))
a.XLabel.String=('time (ms)');
a.YLabel.String=('BF (log scale)');
a.FontSize = 20;


% add the prior plots to the right
size_tiny_plots = (top-bot)/3+8;
oneway=1;
limits=5;
null_line_width=1;
x_ax_on=0;
null_int_start=0.5;
widths = [sqrt(2),1,1/sqrt(2)];
for i = 1:3
    a=axes('Units','Pixels','Position',[900 top-(i-1)*size_tiny_plots-18*i size_tiny_plots size_tiny_plots ],'Color',[co(i+1,1:3),0.5]);
    prior_width=widths(i);
    plot_cauchy(a,null_int_start,oneway,limits,null_line_width,co_bwr,x_ax_on,prior_width)
end


set(gcf, 'InvertHardCopy', 'off'); 
saveas(gcf, '../figures/figure3.png')

%% ***********************************************************
% ************************* Figure 5 *************************
% ************************************************************

% make a plot comparing different null intervals
tv = ds_stacked_realcolour.a.fdim.values{1}*1000;
co= getPyPlot_cMap('gray',5);

% run the BFs
X = ds_stacked_realcolour.samples';
null_ints = [0,0.2,0.5,0.8];
for i = 1:4
    bfs_toplot{i} =bayesfactor_R_wrapper(X,'args',['mu=0.5,rscale="medium",nullInterval=c(' num2str(null_ints(i)) ',Inf)'],'returnindex',1);
end


f=figure(2);clf
f.Position=[f.Position(1:2) 1000 400];f.Resize='off';f.PaperPositionMode='auto';f.Color='w';
top = 300; bot = 60;
a=axes('Units','Pixels','Position',[55 bot 800 top]);

hold on
a=gca;
markers=['-d';'-s';'-o';'-^'];
null_ints_titles={'[0, Inf]','[0.2, Inf]','[0.5, Inf]','[0.8, Inf]'};
for i = 1:4
    stem(a,tv,bfs_toplot{i},markers(i,:),'basevalue',1,'Color','k','LineWidth',0.01,'MarkerSize',8)
    hold on
end

for i = 1:4
    p(i)=plot(tv,bfs_toplot{i},markers(i,:),'MarkerFaceColor',co(i,1:3),'MarkerEdgeColor',co(i,1:3),'LineStyle','None','MarkerSize',8)
end
legend(p,null_ints_titles,'Location','NE','FontSize',24,'Box','off')

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
    a=axes('Units','Pixels','Position',[900 top-(i-1)*size_tiny_plots-9*i size_tiny_plots size_tiny_plots ],'Color',[co(i,1:3),0.5]);
    null_int_start=null_ints(i);
    plot_cauchy(a,null_int_start,oneway,limits,null_line_width,co_bwr,x_ax_on,prior_width)
end
set(gcf, 'InvertHardCopy', 'off'); 
saveas(gcf, '../figures/figure5.png')



%% helper functions
function plot_cauchy(ax,null_int_start,oneway,limits,null_line_width,co,x_ax_on,prior_width)

    pd_cauchy = makedist('Stable','alpha',1,'beta',0,'gam',prior_width,'delta',0);

    x = -limits:.1:limits;
    pdf_cauchy = pdf(pd_cauchy,x);
    hold on
    if oneway==1
        area(x(x>=0),pdf_cauchy(x>=0),'FaceColor',[1,1,1])
    else
        area(x,pdf_cauchy,'FaceColor',[1,1,1])
    end
    area(x(x>=null_int_start),pdf_cauchy(x>=null_int_start),'FaceColor',co(end,1:3))
    ax.XLim=[x(1), x(end)]; ax.YLim = [0,0.5];
    if x_ax_on==0
        ax.XTick=[];
    else
        xlabel('effect size (\delta)')
        ax.XTick = -limits:limits;
        ax.FontSize=14;
    end
    ax.YTick=[];
    ax.Box='on';
    stem(0,max(pdf_cauchy),'.m','Color', co(1,1:3), 'LineWidth',null_line_width,'MarkerSize',1)


end
