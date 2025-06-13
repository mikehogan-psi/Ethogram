function[mfr,sfr,t,fr,frspk,trialspk] = raster_NM(tsp,evt,tpre,tpost,dt,graph)
t = [-tpre:dt:tpost];
N = length(evt);
Nt = length(t);
fr = zeros(N,Nt);
frspk = cell(1,N);
trialspk = cell(1,N);
for n = 1:N
     temp = tsp(find((tsp>evt(n)-tpre)&(tsp<(evt(n)+tpost))))-evt(n);
     frspk{n} = temp;
     fr(n,:) = hist(temp,t)/dt; 
     trialspk{n} = n*ones(length(temp),1);
end
frspk = vertcat(frspk{:});
trialspk = vertcat(trialspk{:});
t = t(2:end-1);
Nt = length(t);
fr = fr(:,2:end-1);
mfr = mean(fr);
sfr = std(fr)/sqrt(N);
if graph&length(frspk) 
    fig = figure;
    set(fig,'Position',[300 300 600 300]);
    %
    h1 = subplot(1,2,1), hold on;
    %plot(frspk,trialspk,'k.','MarkerSize',8);
    p = pcolor(t,1:size(fr,1),fr); set(p,'LineStyle','none');
    xlim([t(1) t(end)]);ylim([1 size(fr,1)]);
    %xlim([t(1) t(end)]);
    %xlabel('Time(s)','FontSize',14,'FontName','Arial');
    ylabel('#Trial','FontSize',14,'FontName','Arial');
    set(h1,'FontSize',14,'FontName','Arial');
    ylim([1 N+1]);
    %
    mfr = mean(fr);
    maxfr = 1.1*max(mfr);
    minfr = 0.9*min(mfr);
    %
    h2 = subplot(1,2,2), hold on;
    bar(t,mfr,'BarWidth',1.01,'LineStyle','none','FaceColor','k');
    set(h2,'FontSize',14,'FontName','Arial');
    xlabel('Time(s)','FontSize',14,'FontName','Arial');
    ylabel('FR(Hz)','FontSize',14,'FontName','Arial');
    xlim([t(1) t(end)]);ylim([minfr maxfr]);
end
