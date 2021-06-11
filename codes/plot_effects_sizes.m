%%
addpath('~/Repository/CommonFunctions/matplotlib/')
addpath('~/Repository/CommonFunctions/distributionPlot/')

%%
load('overview_effectsizes.mat','effects','labels','timevecs')

f = figure(1);clf
f.Position = [f.Position(1:2) 800 600];f.Resize='off';
co = tab10;
yrange = [min(cellfun(@min,effects)) max(cellfun(@max,effects))];
xrange = [min(cellfun(@min,timevecs)) max(cellfun(@max,timevecs))];

prestim=[];
poststim=[];
for e=1:numel(effects)
    prestim(e) = max(effects{e}(timevecs{e}<=0));
    poststim(e) = max(effects{e}(timevecs{e}>0));
end

c1 = co(10,:);
c2 = co(1,:);

newlabels = {' Teichmann et al., 2018 magnitude decoding',...
    ' Teichmann et al., 2020 colour congurency',...
    ' Moerel et al., 2021a attention & decision',...
    ' Moerel et al., 2021b attention & expectation',...
    ' Grootswagers et al., 2021 object attention',...
    ' Robinson et al., 2021 imagined location',...
    ' Grootswagers et al., 2019 mid level features'};
locs = [0 4 5 3 3 6 2]; %how many values for each label
labelstart = cumsum(locs)+(1:numel(locs)); %magic

a1=subplot(1,1,1);
a1.YDir='reverse';
hold on
box off
b=[];
yvals = 1:(numel(prestim)+numel(newlabels));
yvals(labelstart) = [];alpha=.9;
b(1)=barh(yvals,poststim,1,'FaceColor',c2,'FaceAlpha',alpha);
b(2)=barh(yvals,prestim,1,'FaceColor',c1,'FaceAlpha',alpha);
a1.YTick=[];
text(0*labelstart,labelstart,newlabels,'Color','k','FontSize',12)
xlabel('Estimated effect size (d)')
a1.YAxis.Visible='off';
distributionPlot(prestim','xValues',max(yvals)+3,'xyOri','flipped','addSpread',0,'distWidth',3,...
    'histOpt',1,'color',{c1},'showMM',0,'histOri','left');
distributionPlot(poststim','xValues',max(yvals)+3,'xyOri','flipped','addSpread',0,'distWidth',3,...
    'histOpt',1,'color',{c2},'showMM',0,'histOri','left');
xlabel('Maximum effect size (\delta)')
legend(fliplr(b),{'pre stimulus onset','post stimulus onset'},'Location','NE','Box','off','FontSize',16)
a1.YTick = [];
a1.FontSize=16;
ylim([0 max(yvals)+4.5])

fprintf('\nprestim avg: %.4f min: %.4f max: %.4f\n',mean(prestim),min(prestim),max(prestim))
fprintf('posttim   avg: %.4f min: %.4f max: %.4f\n',mean(poststim),min(poststim),max(poststim))

%% save
fn = sprintf('../figures/figure_effects');
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=1;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');
