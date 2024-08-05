function [M_plot, Response_4s, ste_Response] = SingleAnimal_task_DA(mouse_name)

close all
d = load('3_task_DA_response.mat');
dat = d.signal_raw;
d_mouse = dat{mouse_name};
GCaMP = d_mouse(:,1);
opto = d_mouse(:,2);
led = d_mouse(:,3);


%% trial start time: opto & trial type sig 

%%%%%%%% ten 5ms lights in 0.5s %%%%%%%%%%%

% get opto start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opto_on = crossing(opto,[],2); %threshold(mV)
opto_on_ts = (opto_on(1:2:end)).';  %each opto on
opto_off_ts = (opto_on(2:2:end)).';
opto_trial_all = opto_on_ts(1:80:end);
opto_trial = opto_trial_all;

led_on = crossing(led,[],2); %threshold(mV)
led_on_ts = (led_on(1:2:end)).';  %each led on
led_off_ts = (led_on(2:2:end)).';
led_trial_all = led_on_ts(1:4:end);
% led_trial_all = led_on_ts(1:1:end);


% Find elements in list 1 with a difference within 100 in list 2
toRemove = [];
for i = 1:numel(led_trial_all)
    for j = 1:numel(opto_trial)
        if abs(led_trial_all(i) - opto_trial(j)) <= 100
            toRemove = [toRemove, i];
            break; % No need to check further if a match is found
        end
    end
end

% Remove elements from list 1
led_trial = led_trial_all;
led_trial(toRemove) = [];

%% clean photometry signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%

normG_median_divided = analyze_noise_onlyG(GCaMP,led_trial_all);

%% trigger & trial number

trigger = {led_trial; opto_trial};
triggerB = trigger;

Trial_idx = cell(size(trigger));
Trial_idx{2} = toRemove;
temp = 1:30;
temp(toRemove)=[];
Trial_idx{1} = temp;

%% make matrix of GCaMP data

plotdata = normG_median_divided;
plotWin = [-2000:8000];
M_plot = [];
DeltaF = cell(size(trigger));

Response_1s = cell(size(trigger));
Response_4s = zeros(2,1);
response_tbt = []; 
ste_Response = [];
Trial_number = [];
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

    F = mean(rawTraceB(:,1:2000),2);        %using this time window in plotwin for baseline
    deltaF = bsxfun(@minus,rawTrace,F);
    m_plot = mean(deltaF);
    s_plot = std(deltaF)/sqrt(length(ts));

    % average of this trial type
    if i == 1;
        plot_rectangle = true;
    else
        plot_rectangle = false;
    end
    M_plot = [M_plot; m_plot];
    
    DeltaF{i} = deltaF;
    Trial_number = [Trial_number size(deltaF,1)];
    response_1s = deltaF(:,2000:2500); % 2000 is trigger1
    response_4s = deltaF(:,2000:6000); 

    response_1s = mean(response_1s');
    response_4s = mean(response_4s');
    Response_1s{i} = response_1s;
    Response_4s(i,:) = mean(response_4s);
    ste_Response = [ste_Response std(response_1s)/sqrt(length(response_1s))];
end
