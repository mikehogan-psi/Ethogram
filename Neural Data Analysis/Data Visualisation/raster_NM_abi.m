function[mfr,sfr,t,fr,frspk,trialspk] = raster_NM_abi(tsp,evt,tpre,tpost,dt,graph)

% General: this function creates raster plots and peri-stimulus time histograms (PSTHs) 
%          for spike trains aligned to a set of events

% Inputs:
%        - tsp: Spike timestamps (vector of spike times, in seconds).
%        - evt: Event times to align spikes to (e.g., stimulus onset times).
%        - tpre: Time before each event to include in the window (e.g., 0.5 seconds).
%        - tpost: Time after each event to include.
%        - dt: Bin size (time resolution of histogram, in seconds).
%        - graph: Flag (0 or 1) to plot the raster and PSTH.
% 
% Outputs:
%        - mfr: Mean firing rate across trials (Hz).
%        - sfr: Standard error of firing rate across trials (Hz).
%        - t: Time vector (excluding the first and last bin).
%        - fr: Firing rate matrix (trials × time bins).
%        - frspk: Vector of spike times aligned to events (all trials concatenated).
%        - trialspk: Trial numbers for each spike in frspk.

% Trial time

% initializing variables 
    t = [-tpre:dt:tpost]; % vector of time bins centered on the event time
                          % for us eg: [10:1:23] ???
    N = length(evt);      % Number of events/trials
    Nt = length(t);       % Number of time bins
    fr = zeros(N,Nt);     % Firing rate matrix
    frspk = cell(1,N);    % Cell array to hold aligned spikes per trial
    trialspk = cell(1,N); % Trial index for each spike


% STEP 1: Extracting spikes and firing rates within each trial
for n = 1:N
     temp = tsp(find((tsp>evt(n)-tpre)&(tsp<(evt(n)+tpost))))-evt(n); % extract spike event in trial time-window and align to stimulus onset time
     frspk{n} = temp;                                                 % stores aligned spike times
     fr(n,:) = hist(temp,t)/dt;                                       % bin spikes into histogram using hist and convert to firing rate (Hz)
     trialspk{n} = n*ones(length(temp),1);                            % stores spikes of this trial as trial number
end

% STEP 2: Concatenate spike times and trial indices
frspk = vertcat(frspk{:});
trialspk = vertcat(trialspk{:});

% STEP 3: Trim edges (Removes the first and last bins of the histogram to  avoid edge effects)
t = t(2:end-1);
Nt = length(t);
fr = fr(:,2:end-1);

% STEP 4: Compute PSTH statistics
mfr = mean(fr);            % Mean firing rate over trials (PSTH)
sfr = std(fr)/sqrt(N);     % Standard error of mean

% STEP 5: plotting results
if graph&length(frspk) 

    fig = figure;
    set(fig,'Position',[300 300 600 300]);
    
  % Subplot 1: Heatmap of spike counts
    h1 = subplot(1,2,1);
    hold on;
    %plot(frspk,trialspk,'k.','MarkerSize',8);
    p = pcolor(fr); set(p,'LineStyle','none');
    xlim([1 size(fr,2)]);ylim([1 size(fr,1)]);
    %xlim([t(1) t(end)]);
    %xlabel('Time(s)','FontSize',14,'FontName','Arial');
    ylabel('#Trial','FontSize',14,'FontName','Arial');
    set(h1,'FontSize',14,'FontName','Arial');
    ylim([1 N+1]);
    
    mfr = mean(fr);
    maxfr = 1.1*max(mfr);
    minfr = 0.9*min(mfr);
    
  % Subplot 2: PSTH
    h2 = subplot(1,2,2), hold on;
    bar(t,mfr,'BarWidth',1.01,'LineStyle','none','FaceColor','k');
    set(h2,'FontSize',14,'FontName','Arial');
    xlabel('Time(s)','FontSize',14,'FontName','Arial');
    ylabel('FR(Hz)','FontSize',14,'FontName','Arial');
    xlim([t(1) t(end)]);ylim([minfr maxfr]);

end
