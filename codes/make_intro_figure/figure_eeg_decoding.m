%% Introduction figure showing MVPA for time-series data
%   This figure is based on simulated data and has been adapted from:
%   Carlson, Grootswagers, & Robinson, 2020,
%   https://arxiv.org/pdf/1905.04820.pdf 
%   April 2021

%% setup
clearvars;

f=figure(1);clf
f.Position=[0 80 500 635];f.Resize='off';f.PaperPositionMode='auto';f.Color='w';
left = 40;

%% bottom: decoding over time
a=axes('Units','Pixels','Position',[left 30 450 140]);
m = load('cleardecodingdata.mat');
xdata = m.xdata;
xdata = [xdata 605:5:800]
ydata=m.ydata*100;
ydata=[ydata repmat(ydata(end),1,40)-(rand(40,1))']

plot(xdata,ydata,'k','LineWidth',2);hold on
plot(xdata,50+0*ydata,'k:','LineWidth',2)

a.YLim=[45 75];a.YTick=[50 60 70];

h=fill([0 100 100 0],[0 0 1 1],'k','FaceAlpha',.1,'EdgeAlpha',1);

xlabel('time (ms)')
ylabel('classification accuracy (%)')
text(0,.55,'stimulus onset','Rotation',90,'VerticalAlignment','top','HorizontalAlignment','left')
text(100,.55,'stimulus offset','Rotation',90,'VerticalAlignment','bottom','HorizontalAlignment','left')

%% middle: decoding over time
rng(50)
bottom = 200;
NexemplarsPerCat = 24;
arange = [-.5 3.5 -.5 3.5];
x_mu = [1 2]; 
x_var = [.5 .5]; 

y_mu = [2 1]; 
y_var = [.5 .5]; 

x = [x_mu(1)+x_var(1)*randn([NexemplarsPerCat 1]) x_mu(2)+x_var(2)*randn([NexemplarsPerCat 1])];
y = [y_mu(1)+y_var(1)*randn([NexemplarsPerCat 1]) y_mu(2)+y_var(2)*randn([NexemplarsPerCat 1])];

m = fitcdiscr([x;y],[zeros(length(x),1);ones(length(y),1)]);

a=axes('Units','Pixels','Position',[left bottom 200 200]);
h1 = scatter(x(:,1),x(:,2),'ok','MarkerFaceColor',[0.1725,0.6275,0.1725],'MarkerEdgeColor',[0.1725,0.6275,0.1725]);
hold on
h2 = scatter(y(:,1),y(:,2),'^k','MarkerFaceColor',[0.8392,0.1529,0.1569],'MarkerEdgeColor',[0.8392,0.1529,0.1569]);

axis(arange)

a=gca;
c = m.Coeffs(1,2).Const;
bx = m.Coeffs(1,2).Linear(1);
by=m.Coeffs(1,2).Linear(2);
h=ezplot(@(x,y) c+bx*x+by*y);h.LineWidth=2;h.LineColor='k';
a.XTick=[];a.YTick=[];a.XLabel.String='channel 1';a.YLabel.String='channel 2';a.Title=[];

legend([h1 h2],{'Green','Red'},'Location','NE','FontSize',12)

train_labels  = repmat((1:4)',[NexemplarsPerCat/4 1]);

for i = 1:4    
    ms = 4;
    a = axes('Units','Pixels','Position',[left+220+(i-1)*60 bottom+150 50 50]);
    plot(x(train_labels~=i,1),x(train_labels~=i,2),'ok','MarkerFaceColor',[0.1725,0.6275,0.1725],'MarkerEdgeColor',[0.1725,0.6275,0.1725],'MarkerSize',ms)
    hold on
    plot(y(train_labels~=i,1),y(train_labels~=i,2),'^k','MarkerFaceColor',[0.8392,0.1529,0.1569],'MarkerEdgeColor',[0.8392,0.1529,0.1569],'MarkerSize',ms)
    axis(arange)
    m = fitcdiscr([x(train_labels~=i,:);y(train_labels~=i,:)],...
        [zeros(sum(train_labels~=i),1);ones(sum(train_labels~=i),1)]);

        
    c = m.Coeffs(1,2).Const;
    bx = m.Coeffs(1,2).Linear(1);
    by=m.Coeffs(1,2).Linear(2);
    h=ezplot(@(x,y) c+bx*x+by*y);h.LineWidth=2;h.LineColor='k';
    a.XTick=[];a.YTick=[];a.XLabel=[];a.YLabel=[];a.Title=[];a.Box='off';
    
    title('train')
    
    a = axes('Units','Pixels','Position',[left+220+(i-1)*60 bottom+75 50 50]);
    plot(x(train_labels==i,1),x(train_labels==i,2),'ok','MarkerFaceColor',[0.1725,0.6275,0.1725],'MarkerEdgeColor',[0.1725,0.6275,0.1725],'MarkerSize',ms)
    hold on
    plot(y(train_labels==i,1),y(train_labels==i,2),'^k','MarkerFaceColor',[0.8392,0.1529,0.1569],'MarkerEdgeColor',[0.8392,0.1529,0.1569],'MarkerSize',ms)
    axis(arange)
    h=ezplot(@(x,y) c+bx*x+by*y);h.LineWidth=2;h.LineColor='k';
    a.XTick=[];a.YTick=[];a.XLabel=[];a.YLabel=[];a.Title=[];a.Box='off';
    
    title('test')
end

for i=1:4
    for j=1:4
        a = axes('Units','Pixels','Position',[left+220+(j-1)*60 bottom+(i-1)*12.5 50 12.5]);
        axis([0 1 0 1])
        if i==j
            fill([0 1 1 0],[0 0 1 1],'k')
            t = text(0.5,0.5,'TEST','Color','w');
        else
            t = text(0.5,0.5,'TRAIN');
        end
        t.VerticalAlignment='middle';
        t.HorizontalAlignment='center';
        a.Box='on';a.XTick=[];a.YTick=[];
        
        if i==4
            title(sprintf('iteration %i',j))
        end
    end
end
drawnow
    
%% top: scalp & eeg
bottom = 430;

% Subplot B
a=axes('Units','Pixels','Position',[left+220 bottom 230 200]);

m=load('dat.mat');dat=m.dat;dat = dat(1:10,:)
plot((repmat(50*(1:size(dat,1))'+150,1,size(dat,2))+dat)','k')
axis([1 2000 0 800])
a.XTick=[];a.YTick=[];

stim_onsets = linspace(a.XLim(1),a.XLim(end),7);
hold on
stim_n = ['01';'13';'05';'03';'15';'13'];
x_y_ratio = (a.YLim(end)-a.YLim(1))/(a.XLim(end)-a.XLim(1));
stimsize = 150;
for i  = 1:length(stim_onsets)-2
    stim = imread(['./stimuli/' stim_n(i,:) '.png']);
    image(flipud(stim), 'XData', [stim_onsets(i+1)-stimsize/x_y_ratio/2 stim_onsets(i+1)+stimsize/x_y_ratio/2], 'YData', [a.YLim(1)+10 a.YLim(1)+stimsize])
    plot([stim_onsets(i+1),stim_onsets(i+1)],[a.YLim(end)-100,a.YLim(1)+130],'Color',[0.5,0.5,0.5],'LineWidth',2)
end

title('Channel voltage over time')
a.Title.Position(2) = 720;
a.Title.FontSize = 14;

% Subplot A

ax=axes('Units','Pixels','Position',[left-20 bottom-10 200 200]);
stim = imread('brain_transparent.png');
image(flipud(stim),'XData',[0.1,0.9],'YData',[0,0.8])
axis('off')
hold on
% plot some random sensors & cables
locs_sens = [0.7,0.73;0.8,0.61;0.87,0.49];
for i  = 1:3
    plot(locs_sens(i,1),locs_sens(i,2),'ko','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize',10,'MarkerEdgeColor','k','LineWidth',2)
    plot([locs_sens(i,1),locs_sens(i,1)],[locs_sens(i,2),locs_sens(i,2)+0.12],'Color',[153, 204, 255]/255,'LineWidth',2)
    plot([locs_sens(i,1),1.2],[locs_sens(i,2)+0.12,locs_sens(i,2)+0.12],'Color',[153, 204, 255]/255,'LineWidth',2)
end

plot(0,0.92,'ko','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize',10,'MarkerEdgeColor','k','LineWidth',2)
text(0.05,0.92, 'MEG sensor / EEG channel','FontSize',12)

ax.XLim=[0,1];ax.YLim=[0,1];ax.YDir = 'normal';ax.Clipping='off';


%% annotations for all subplots
fs=20;
a=annotation('textbox','Units','Pixels','Position',[0,540,50,100],'String','A','FontSize',fs,'FontWeight','bold','LineStyle','none');
a=annotation('textbox','Units','Pixels','Position',[235,540,50,100],'String','B','FontSize',fs,'FontWeight','bold','LineStyle','none');
a=annotation('textbox','Units','Pixels','Position',[0,330,50,100],'String','C','FontSize',fs,'FontWeight','bold','LineStyle','none');
a=annotation('textbox','Units','Pixels','Position',[235,330,50,100],'String','D','FontSize',fs,'FontWeight','bold','LineStyle','none');
a=annotation('textbox','Units','Pixels','Position',[0,90,50,100],'String','E','FontSize',fs,'FontWeight','bold','LineStyle','none');

drawnow

%% save
fn = sprintf('../../figures/figure1');
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=0;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

