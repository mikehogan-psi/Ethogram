function[mfr,sfr,t,fr,frspk,trialspk,spikes_per_trial, binned_skp_counts] = raster_NM(tsp,evt,tpre,tpost,dt,graph)

% General: this function creates raster plots and peri-stimulus time histograms (PSTHs) 
%          for spike trains aligned to a set of events

% Inputs:
%        - tsp: Spike timestamps (vector of spike times, in seconds).
%        - evt: Event times to align spikes to (e.g., stimulus onset times).
%        - tpre: Time before each event to include in the window (e.g., 0.5 seconds).
%        - tpost: Time after each event to include.
%        - dt: Bin size (time resolution of histogram, in seconds).
%        - graph: Flag (0 or 1) to plot the raster and PSTH.



t = [-tpre:dt:tpost];
N = length(evt);
Nt = length(t);
fr = zeros(N,Nt);
binned_skp_counts = zeros(N,Nt);
frspk = cell(1,N);
trialspk = cell(1,N);


for n = 1:N
     temp = tsp(find((tsp>evt(n)-tpre)&(tsp<(evt(n)+tpost))))-evt(n);
     frspk{n} = temp;
     binned_skp_counts(n,:) = hist(temp,t);
     fr(n,:) = binned_skp_counts(n,:)/dt; 
     trialspk{n} = n*ones(length(temp),1);
end

spikes_per_trial = frspk;
frspk = vertcat(frspk{:});
trialspk = vertcat(trialspk{:});
t = t(2:end-1);
Nt = length(t);
fr = fr(:,2:end-1); %alter original scrip to remove -1 if needed
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

