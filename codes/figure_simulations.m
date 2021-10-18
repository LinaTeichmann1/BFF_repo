%%
% https://www.cosmomvpa.org/
addpath('~/CoSMoMVPA/mvpa/')
% https://au.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps
addpath('~/Repository/CommonFunctions/matplotlib/')
% https://au.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m
addpath('~/Repository/CommonFunctions/distributionPlot/') 

%%
nboot=1000;
nsubvec = 2:100;
eff = flipud((0:.1:1)');
efs_real=[];efs_sim=[];
tdd = tempdir;
bfs_sim = [];
bf_args = { 'mu=0,rscale="medium",nullInterval=c(0.5,Inf)',...
            'mu=0,rscale="medium",nullInterval=c(0,Inf)'};
for a = 1:numel(bf_args)
    fprintf('Simulating %s\n',bf_args{a});
    fprintf('Distributing tasks in tempdir: %s\n',tdd);
    rng(1);
    clear F;F(1:numel(nsubvec),1)=parallel.FevalFuture;
    for b = 1:numel(nsubvec)
        efs_real(b,:) = repelem(1:numel(eff),1,nboot)';
        X = eff(efs_real(b,:)) + randn(numel(eff)*nboot,nsubvec(b));
        efs_sim(b,:) = mean(X')./std(X');
        td = sprintf('%s/b%i',tdd,b);
        if ~exist(td,'dir')
            mkdir(td);
        end
        F(b) = parfeval(@(x) bayesfactor_R_wrapper(x,...
            'args',bf_args{a},'returnindex',1,'tempdir',td),1,X);
    end
    fprintf('Collecting results\n');
    cc = clock();mm='';
    for b = 1:numel(nsubvec)
        [i,v] = fetchNext(F);
        bfs_sim(a,i,:) = v;
        mm=cosmo_show_progress(cc,b/numel(nsubvec),sprintf('%i/%i',b,numel(nsubvec)),mm);
    end
end
fprintf('\nFinished\n');
save('bf_simulations.mat','bf_args','eff','bfs_sim','efs_real','efs_sim','nsubvec');

%% plot simulations
load('bf_simulations.mat');

f=figure(2);clf
f.Position = [f.Position(1:2) 800 1000];f.Resize='off';
plotnr = 0;

bf_args = {'interval = [0.5,Inf]','','interval = [0,Inf]'};

for aa = [1 2]
    acell = {};
    for i=1:2
        plotnr=plotnr+1;
        a=subplot(4,2,plotnr);hold on
        a.FontSize=12;
        a.YScale = 'log';
        a.YAxis.MinorTick='off';
        if i==1 
            a.XLim = [2 30];
        else
            a.XLim = [30 99];
        end
        co = tab10;
        if i==1
            ylim(10.^[-5+3*(aa>2) 5])
            a.YTick = 10.^(-4:4);
        else
            ylim(10.^[-16+14*(aa>2) 16])
            a.YTick = 10.^(-15:3:15);
        end
        
        h=[];
        fill([a.XLim fliplr(a.XLim)],reshape(([1/10 10]'+0*a.XLim)',[],1),'k','FaceAlpha',.3,'LineStyle','none');
        for e=1:numel(eff)
            idx = efs_real(1,:)==e;
            x = squeeze(bfs_sim(aa,:,idx))';
            mu = median(x);
            ci = prctile(x,[5 95]);
            h(e) = plot(nsubvec,mu,'Color',co(e,:),'LineWidth',2,'DisplayName',sprintf('d=%.2f',eff(e)));
            if aa==1 || mod(e,2)
                text(nsubvec(a.XLim(2)),mu(a.XLim(2)),sprintf('d=%.1f',eff(e)),'Color',co(e,:))
            end
        end
        xlabel('number of subjects')
        if mod(plotnr,2)
            ylabel('Bayes Factor')
        end
        acell{i} = a;
    end
    for i=1:2
        plotnr=plotnr+1;
        a=subplot(4,2,plotnr);hold on
        a.FontSize=12;
        
        co = tab10;
        h=[];
        a.XLim=[0 12];
        fill([a.XLim fliplr(a.XLim)],reshape(([-1 1]'+0*a.XLim)',[],1),'k','FaceAlpha',.3,'LineStyle','none');
        y = [];
        for e=1:numel(eff)
            idx = efs_real(1,:)==e;
            y(:,e) = squeeze(bfs_sim(aa,nsubvec==(30+(i*70-70)),idx))';
        end
        for e=1:numel(eff)
            h=distributionPlot(log(y(:,e)),'addSpread',0,...
                'xValues',12-e,'color',co(e,:),'showMM',3);
            h{2}.Marker='+';
            h{2}.Color='k';
        end
        if i==1
            a.YTick = -16:4:24;
        else
            a.YTick = -40:10:60;
        end
        a.YTickLabel = sprintf('10^{%i}\n',a.YTick);
        
        a.YLim = a.YTick([1 end]);
        
        a.XTick = 1:11;
        a.XTickLabel = flipud(eff);
        xlabel('effect size (\delta)')
        if mod(plotnr,2)
            ylabel('Bayes Factor')
        end
        text(.5,.9*a.YLim(end),sprintf('n=%i',(30+(i*70-70))))
    end
end
annotation('textbox',[.065,.91,.4,.04],'String','A   With interval','FontSize',17,'FontWeight','bold','LineStyle','none')
annotation('textbox',[.065,.471,.4,.04],'String','B   Without interval','FontSize',17,'FontWeight','bold','LineStyle','none')

%% save
fn = sprintf('../figures/figure_simulations');
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=1;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

