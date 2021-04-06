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
    ax.XLim=[x(1), x(end)]; ax.YLim = [0,0.7];
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