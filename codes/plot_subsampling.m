%%
addpath('~/Repository/CommonFunctions/matplotlib')
load('subsampling_results.mat')
load('color_bwr.mat','co')

r = 10*[-1 1];
prct = fliplr([10 25 50 75 100]);
tv = 1000*all_res_subsampled{1}{1}.a.fdim.values{1};
f=figure(1);clf
f.Position = [f.Position(1:2) 800 800];
plotnr=0;
for s = 1:2:5
    plotnr=plotnr+1;
    a=subplot(3,2,plotnr);
    imagesc(tv,1:5,log([bfs{s}{:}])',r)
    a.YTick=1:numel(prct);
    a.YTickLabel = prct;
    a.FontSize = 14;
    title(sprintf('n=%i',3+s*3))
    a.Colormap = co;
    c=colorbar;
    c.Label.String = 'log BF';
    ylabel('trials (%)')
    if s==5
        xlabel('time (ms)')
    end
    plotnr=plotnr+1;
    a=subplot(3,2,plotnr);
    imagesc(tv,1:5,vertcat(sig{s}{:}),[0 1])
    a.YTick=1:numel(prct);
    a.YTickLabel = prct;
    a.FontSize = 14;
    title(sprintf('n=%i',3+s*3))
    a.Colormap = [1 1 1;0 0 0;];
    c=colorbar;
    c.Label.String = 'p<.05';
    ylabel('trials (%)')
    c.Ticks = [0 1];c.TickLabels = {'plos one','nature'};
    if s==5
        xlabel('time (ms)')
    end
end
    
%% save
fn = sprintf('../figures/figure_subsampling');
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=1;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

    