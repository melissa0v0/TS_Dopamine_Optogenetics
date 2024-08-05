function [M_plot, Response, ste_Response] = SingleAnimal_masking(mouse_name, m)
%% initialize

close all
if m == 'm'
    d = load('0_no_masking.mat');
elseif m == 'n'
    d = load('0_no_masking.mat');
end

dat = d.signal_raw;
d_mouse = dat{mouse_name};
GCaMP = d_mouse(:,1);
opto = d_mouse(:,2);
led = d_mouse(:,3);

% plot all signals for the whole session %
% f1=figure();
% plot(GCaMP, 'g');hold on
% plot(led,'b');hold on
% plot(opto,'r');hold on
% 
% title('photometry(g), opto(r 20Hz 0.5s), LED(b)');

%% 
% get stim start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

led_on = crossing(led,[],2); %threshold(mV)
led_on_ts = (led_on(1:2:end)).';  %each opto on
led_off_ts = (led_on(2:2:end)).';


%% trial start time: opto & trial type sig 

%%%%%%%% ten 5ms lights in 0.5s %%%%%%%%%%%

% get opto start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trial_type_sig_on = crossing(opto,1:length(opto),2.5);
trial_type_sig_off = (trial_type_sig_on(2:2:end)).';
trial_type_sig_on = (trial_type_sig_on(1:2:end)).';

trial_type_sig_on(:,2) = trial_type_sig_on(:,1);
trial_type_sig_on(:,2) = 0;

for i = 1:size(trial_type_sig_on,1)
temp_B = 1;
temp_A = trial_type_sig_on(i,1);    
    for f = i+1:size(trial_type_sig_on,1)
        if trial_type_sig_on(f,1) < temp_A + 5000  %identify multiple signals within 1s
        temp_B = temp_B + 1;  
        end
    end
trial_type_sig_on(i,2) = temp_B;  %number of signals within 1s after trial_type_sig_on(i,1)
end


% delete duplicates of the same trial
for i = 1:size(trial_type_sig_on,1)
    temp_B = 1;
    temp_A = trial_type_sig_on(i,1);
    for f = i+1:size(trial_type_sig_on,1)
        if trial_type_sig_on(f,1) < temp_A + 4000
        trial_type_sig_on(f,:) = 0;
        end
    end
end


ind = find(trial_type_sig_on(:,1)>0);
trial_type_sig_on = trial_type_sig_on(ind,:);

trial_type_sig_all = [];
for iii = 1: max(trial_type_sig_on(:,2))
    trial_type = find(trial_type_sig_on(:,2)==iii);
    trial_type_ts{iii} = trial_type_sig_on(trial_type);  % use first TTL
    trial_type_sig_all = [trial_type_sig_all;trial_type_sig_on(trial_type)];
end
trial_type_sig_all = sort(trial_type_sig_all); 


% f2=figure();
% plot(opto,'m'); hold on
% for iii = 1: max(trial_type_sig_on(:,2))
%     if size(trial_type_ts{iii},1)>0
%     plot(trial_type_ts{iii},ones,'o')
%     end
% end


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean photometry signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove 60Hz noise
normG_median_divided = analyze_noise_onlyG(GCaMP,trial_type_sig_all);

% figure();
% plot(normG_median_divided);
% title('GCamp after denoise - smoothing - decay correction')

%% make matrix of GCaMP data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trigger = {trial_type_ts{10}};
% trigger = {trial_type_ts{5}, trial_type_ts{10},trial_type_ts{20}};
% trigger = {trial_type_ts{40}, trial_type_ts{80},trial_type_ts{160}};
% trigger = {trial_type_ts{30};trial_type_ts{60};trial_type_ts{120}};
triggerB = trigger;

% figure
% rectangle('Position', [0, -1, 4000, 11], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on

% plotColors = color_select(4);
% plotColors = {'r','g','b','c','m'};
% legend_name = {'10Hz','20Hz','40Hz'};
plotdata = normG_median_divided;
plotWin = [-2000:8000];
M_plot = [];
DeltaF = [];
Trial_number = [];

Response = [];
ste_Response = [];
for i = 1:length(trigger)
    ts = round(trigger{i});
    if ~exist('triggerB');
        triggerB = trigger;
    end
    tsB = round(triggerB{i});

    ind = find( tsB+ plotWin(1)>0,1,'first');
    ind2 = find( ts+ plotWin(end)< length(plotdata),1,'last');
    ts = ts(ind:ind2);
    plotind = bsxfun(@plus, repmat(plotWin,length(ts),1),ts);
    rawTrace = plotdata(plotind);

    tsB = tsB(ind:ind2);
    plotind = bsxfun(@plus, repmat(plotWin,length(ts),1),tsB);
    rawTraceB = plotdata(plotind); 

    F = mean(rawTraceB(:,1:1900),2);        %using this time window in plotwin for baseline
    deltaF = bsxfun(@minus,rawTrace,F);
    m_plot = mean(deltaF);
    s_plot = std(deltaF)/sqrt(length(ts));
    % errorbar_patch(plotWin,m_plot,s_plot,plotColors{i});
    M_plot = [M_plot; m_plot];
    
    DeltaF = [DeltaF;deltaF];
    Trial_number = [Trial_number size(deltaF,1)];
    response = deltaF(:,2000:3000); %2000 is trigger1
    response = mean(response');
    Response = [Response mean(response)];
    ste_Response = [ste_Response std(response)/sqrt(length(response))];
end

% % legend('0','1.2','2.5','5','10','20')
% xlabel('time - opto stimulation(ms)')
% title('DA sensor')
% legend(legend_name)
% 
% figure
% bar(Response)
% hold on
% errorbar(Response, ste_Response,'b.')
% hold on
% title('DA sensor')
% ylabel('response to 0.5s 12mW stimulation (0-1s)')
% h=gca;
% h.XTickLabel = {'10Hz','20Hz','40Hz'};
end