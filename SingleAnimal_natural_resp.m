function [M_plot_m, Response_4s_m, ste_Response] = SingleAnimal_natural_resp(mouse_name)
%% initialize

close all
d = load('1_natrual_novel_response.mat');
dat = d.signal_raw;
d_mouse = dat{mouse_name};
GCaMP = d_mouse(:,1);
opto = d_mouse(:,2);
led = d_mouse(:,3);


%% trial start time

led_on = crossing(led,[],2); %threshold(mV)
led_on_ts = (led_on(1:2:end)).';  %each led on
led_off_ts = (led_on(2:2:end)).';
led_trial_all = led_on_ts(1:4:end);

trigger_m = {led_trial_all};
triggerB_m = trigger_m;

%% Clean photometry 

normG_median_divided_m = analyze_noise_onlyG(GCaMP,led_trial_all);


%% GCaMP data for first 5 trials

plotdata = normG_median_divided_m;
trigger_m_last = {trigger_m{1}(1:end)};
DeltaF_m = cell(size(trigger_m_last));
M_plot_m = [];
Response_4s_m = [];
ste_Response = [];
plotWin = [-2000:8000];

for i = 1:length(trigger_m_last)
    ts = round(trigger_m_last{i});
    if ~exist('triggerB_m');
        triggerB_m = trigger_m_last;
    end
    tsB = round(triggerB_m{i});

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
    M_plot_m = [M_plot_m; m_plot];

    % figure()
    % subplot(2,3,1)
    % errorbar_patch(plotWin,m_plot,s_plot,'k',false);

    DeltaF_m{i} = deltaF;
    response_4s_m = deltaF(:,2000:6000); 
    
    % response = (response1 + response2 + response3 + response4) ./ 4;
    response_4s_m = mean(response_4s_m');
    Response_4s_m = [Response_4s_m response_4s_m];
    ste_Response = [ste_Response std(response_4s_m)/sqrt(length(response_4s_m))];
end
end