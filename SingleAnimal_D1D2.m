function [Response_1s, Response_4s, M_plot,Response_1s_m, Response_4s_m, DeltaF_m,sig] = SingleAnimal_D1D2(mouse_name,day)

if strcmp(mouse_name(1:3),'ADF')
    cd("D:\journey\Harvard\code\3_ADF_Stim");
else
    cd("D:\journey\Harvard\code\3_TAF_Stim\TAF7_13_0321")
end
addpath 'D:\journey\Harvard\code\Analysis'

file_ID = fopen([mouse_name,'_',day,'t_Analog'],'r');
CC_analog = fread(file_ID, inf, 'double', 0, 'b');

A = reshape(CC_analog, 2, [ ]);
B = reshape (A, 2, 12, [ ]); % A, 2, num of channels

file_ID2 = fopen([mouse_name,'_',day,'_Analog'],'r');
CC_analog2 = fread(file_ID2, inf, 'double', 0, 'b');

A2 = reshape(CC_analog2, 2, [ ]);
B2 = reshape (A2, 2, 12, [ ]); % A, 2, num of channels

GCaMP = B (:,1,:);
GCaMP = reshape (GCaMP,[],1)+4;
GCaMPm = B2 (:,1,:);
GCaMPm = reshape (GCaMPm,[],1)+4;

opto = B (:,5,:);
opto = reshape (opto,[],1);
optom = B2 (:,5,:);
optom = reshape (optom,[],1);

led = B (:,12,:);
led = reshape (led,[],1);
ledm = B2 (:,12,:);
ledm = reshape (ledm,[],1);


%% trial start time: opto & trial type sig 

%%%%%%%% ten 5ms lights in 0.5s %%%%%%%%%%%

% get opto start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opto_on = crossing(opto,[],2); %threshold(mV)
opto_on_ts = (opto_on(1:2:end)).';  %each opto on
opto_off_ts = (opto_on(2:2:end)).';
opto_trial_all = opto_on_ts(1:80:end);
opto_trial = opto_trial_all(1:15);

led_on = crossing(led,[],2); %threshold(mV)
led_on_ts = (led_on(1:2:end)).';  %each led on
led_off_ts = (led_on(2:2:end)).';
led_trial_all = led_on_ts(1:4:end);

led_on = crossing(ledm,[],2); %threshold(mV)
led_on_ts = (led_on(1:2:end)).';  %each led on
led_off_ts = (led_on(2:2:end)).';
led_trial_m = led_on_ts(1:4:end);
% led_trial_m = led_trial_m(1:20);

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
led_trial = led_trial(1:15);

% remove outlier
% opto_trial([11,10])=[];
% led_trial([13,12,11]) =[];

%% plot all signals for the whole session %
% f1=figure();
% plot(GCaMP, 'g');
% hold on;
% plot(led,'b');
% hold on;
% plot(opto,'r')
% 
% title('photometry(g), opto(r 20Hz 0.5s), LED(b)');

%% Clean photometry 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean photometry signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%

normG_median_divided_ = analyze_noise_onlyG(GCaMP,led_trial_all);
normG_median_divided_m_ = analyze_noise_onlyG(GCaMPm,led_trial_m);

d_ = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',12.6,'HalfPowerFrequency2',12.9, ...
               'DesignMethod','butter','SampleRate',1000);
normG_median_divided = filtfilt(d_,normG_median_divided_);
normG_median_divided_m = filtfilt(d_,normG_median_divided_m_);

%% trigger & trial number

trigger = {led_trial; opto_trial};
triggerB = trigger;
trigger_m = {led_trial_m};
triggerB_m = trigger_m;

Trial_idx = cell(size(trigger));
Trial_idx{2} = toRemove;
temp = 1:30;
temp(toRemove)=[];
Trial_idx{1} = temp;

%% make matrix of GCaMP data


traces_all = figure('Name',mouse_name);

% Create color gradient using parula colormap
numColors = 15;
plotColors_tbt = parula(numColors);
plotColors = {'k','r','y','g','c','b'};
plotdata = normG_median_divided;
plotWin = [-2000:8000];
M_plot = [];
DeltaF = cell(size(trigger));

Response_1s = cell(size(trigger));
Response_4s = cell(size(trigger));
ste_Response = [];
Trial_number = [];

figure(traces_all)
subplot(2,3,1)
rectangle('Position', [0, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
rectangle('Position', [1000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
rectangle('Position', [2000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
rectangle('Position', [3000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
ylim([-3, 7]);


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
%     deltaF = bsxfun(@minus,rawTrace,F)/2.5;
%     deltaF = bsxfun(@rdivide,rawTrace,F);
    m_plot = mean(deltaF);
%      m_plot = mean(deltaF(end-20:end,:));
%      m_plot = deltaF(6,:);
    s_plot = std(deltaF)/sqrt(length(ts));
%     errorbar_patch(plotWin,m_plot,s_plot,plotColors(i,:));

    % % single trial raw traces
    % for f = 1:size(deltaF,1)
    %     errorbar_patch(plotWin,deltaF(f,:),zeros(1,size(deltaF,2)),plotColors_tbt(f,:))
    % end
    
    % average of this trial type
    if i == 1;
        plot_rectangle = true;
    else
        plot_rectangle = false;
    end
    figure(traces_all)
    subplot(2,3,1)
    errorbar_patch(plotWin,m_plot,s_plot,plotColors{i},false);
    M_plot = [M_plot; m_plot];
    
    DeltaF{i} = deltaF;
    Trial_number = [Trial_number size(deltaF,1)];
    response_1s = deltaF(:,2000:3000); % 2000 is trigger1
    response_4s = deltaF(:,2000:6000); 
    % response2 = deltaF(:,3000:3500); % 4000 is trigger2
    % response3 = deltaF(:,4000:4500); % 6000 is trigger3
    % response4 = deltaF(:,5000:5500); % 6000 is trigger3
    
    % response = (response1 + response2 + response3 + response4) ./ 4;
    response_1s = mean(response_1s');
    response_4s = mean(response_4s');
    Response_1s{i} = response_1s;
    Response_4s{i} = response_4s;
    ste_Response = [ste_Response std(response_1s)/sqrt(length(response_1s))];
end

% legend('0','1.2','2.5','5','10','20')
xlabel('time - opto stimulation(ms)')
title('GCaMP')

%% GCaMP data for first 5 trials

plotdata = normG_median_divided_m;
trigger_m_last = {trigger_m{1}(1:end)};
DeltaF_m = cell(size(trigger_m_last));
M_plot_m = [];
Response_1s_m = cell(size(trigger_m_last));
Response_4s_m = cell(size(trigger_m_last));
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
%     deltaF = bsxfun(@minus,rawTrace,F)/2.5;
%     deltaF = bsxfun(@rdivide,rawTrace,F);
    m_plot = mean(deltaF);
%      m_plot = mean(deltaF(end-20:end,:));
%      m_plot = deltaF(6,:);
    s_plot = std(deltaF)/sqrt(length(ts));
    M_plot_m = [M_plot_m; m_plot];

    figure(traces_all)
    subplot(2,3,1)
    errorbar_patch(plotWin,m_plot,s_plot,'y',false);

    DeltaF_m{i} = deltaF;
    Trial_number = [size(deltaF,1) Trial_number];
    response_1s_m = deltaF(:,2000:3000); % 2000 is trigger1
    response_4s_m = deltaF(:,2000:6000); 
    % response2 = deltaF(:,3000:3500); % 4000 is trigger2
    % response3 = deltaF(:,4000:4500); % 6000 is trigger3
    % response4 = deltaF(:,5000:5500); % 6000 is trigger3
    
    % response = (response1 + response2 + response3 + response4) ./ 4;
    response_1s_m = mean(response_1s_m');
    response_4s_m = mean(response_4s_m');
    Response_1s_m{i} = response_1s_m;
    Response_4s_m{i} = response_4s_m;
end



%% raster plot

figure(traces_all)
subplot(2,3,4)
% scrsz = get(groot,'ScreenSize');
% figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5])
%1) bin the data
trialNum = size(DeltaF_m{1},1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(DeltaF_m{1}(:,1:length_x),trialNum, binSize,[]),2));
%imagesc(binedF,[-20 20]);
imagesc(binedF,[-7 7]);
colormap yellowblue
xlabel('time');
ylabel('trials')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
title('first twenty')


figure(traces_all)
titles = {'sensory','sensory + Opto'};
for i = 1:length(trigger)
    % scrsz = get(groot,'ScreenSize');
    % figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5])
    subplot(2,3,4+i)
    %1) bin the data
    trialNum = size(DeltaF{i},1); binSize = 100;
    length_x = plotWin(end)-plotWin(1);
    binedF = squeeze(mean(reshape(DeltaF{i}(:,1:length_x),trialNum, binSize,[]),2));
    %imagesc(binedF,[-20 20]);
    imagesc(binedF,[-7 7]);
    colormap yellowblue
    xlabel('time - water (s)');
    ylabel('trials')
    h=gca;
    h.XTick = 0:10:(length_x/binSize);
    h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
    h.YTickLabel = {};
    title(titles{i})
end



%% time course & histrgrams

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1s response %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coeff = cell(size(trigger));
% time_course=figure();
for i = 1:length(trigger)
    idx = Trial_idx{i}+Trial_number(1);
    figure(traces_all)
    subplot(2,3,2)
    scatter(idx,Response_1s{i},'MarkerFaceColor', plotColors{i}, 'MarkerEdgeColor', 'none'); hold on

    % Fit a line to the data using linear regression
    coefficients = polyfit(idx, Response_1s{i}, 1); % Linear regression coefficients

    x_values = idx;
    y_values = polyval(coefficients, x_values);
    plot(x_values, y_values, 'Color', plotColors{i}, 'LineWidth', 2);

    TSS = sum((Response_1s{i} - mean(Response_1s{i})).^2);
    RSS = sum((Response_1s{i} - y_values).^2);
    R2 = 1 - RSS/TSS;
    R2_text = ['Slope = ' num2str(coefficients(1)) ', R^2 = ' num2str(R2)];
    position = round(0.95^i,2);
    yl=ylim;
    text(15, yl(2)-0.2*i, R2_text, 'FontSize', 12,'Color', plotColors{i});

end
figure(traces_all)
subplot(2,3,2)
scatter(1:Trial_number(1),response_1s_m,'MarkerFaceColor','k', 'MarkerEdgeColor', 'none'); hold on
coefficients = polyfit(1:Trial_number(1), Response_1s_m{1}, 1); % Linear regression coefficients
x_values = 1:Trial_number(1);
y_values = polyval(coefficients, x_values);
plot(x_values, y_values, 'Color', 'k', 'LineWidth', 2);
title("Response 0.5s time course")


histogram_1s = figure('Name',mouse_name);
% Loop through each element of the cell array
for h = 1:2

    % Histo 
    figure(histogram_1s)
    subplot(2,3,1)
    histogram(Response_1s{h}, 'FaceColor', plotColors{h}, 'BinWidth', 0.25);
    xlabel('Average Z-score');
    ylabel('Trial Counts');
    hold on;
    
    % Cumulative histo
    subplot(2,3,2)
    [counts1, edges1] = histcounts(Response_1s{h},'BinWidth', 0.25);
    stairs(edges1(1:end-1), cumsum(counts1)/sum(counts1), 'LineWidth', 2,'Color',plotColors{h});
    hold on;
end
% Adjust labels and title
title('Histogram of Response 0.5s');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4s response %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coeff = cell(size(trigger));
% time_course2=figure();
for i = 1:length(trigger)
    idx = Trial_idx{i}+Trial_number(1);
    figure(traces_all)
    subplot(2,3,3)
    scatter(idx,Response_4s{i},'MarkerFaceColor', plotColors{i}, 'MarkerEdgeColor', 'none')
    hold on; % Hold the current plot

    % Fit a line to the data using linear regression
    coefficients = polyfit(idx, Response_4s{i}, 1); % Linear regression coefficients

    x_values = idx;
    y_values = polyval(coefficients, x_values);
    plot(x_values, y_values, 'Color', plotColors{i}, 'LineWidth', 2);

    TSS = sum((Response_4s{i} - mean(Response_4s{i})).^2);
    RSS = sum((Response_4s{i} - y_values).^2);
    R2 = 1 - RSS/TSS;
    R2_text = ['Slope = ' num2str(coefficients(1)) ', R^2 = ' num2str(R2)];
    position = round(0.95^i,2);
    yl=ylim;
    text(15, yl(2)-0.2*i, R2_text, 'FontSize', 12,'Color', plotColors{i});

end
figure(traces_all)
subplot(2,3,3)
scatter(1:Trial_number(1),response_4s_m,'MarkerFaceColor','k', 'MarkerEdgeColor', 'none'); hold on
coefficients = polyfit(1:Trial_number(1), Response_4s_m{1}, 1); % Linear regression coefficients
x_values = 1:Trial_number(1);
y_values = polyval(coefficients, x_values);
plot(x_values, y_values, 'Color', 'k', 'LineWidth', 2);
title("Response 4s time course")

% histogram_4s = figure('Name',mouse_name);
% Loop through each element of the cell array
for h = 1:2
    % Plot histogram for the current element with a different color
    figure(histogram_1s);
    subplot(2,3,4)
    histogram(Response_4s{h}, 'FaceColor', plotColors{h}, 'BinWidth', 0.25);
    hold on;
    subplot(2,3,5)
    [counts1, edges1] = histcounts(Response_4s{h},'BinWidth', 0.25);
    stairs(edges1(1:end-1), cumsum(counts1)/sum(counts1), 'LineWidth', 2,'Color',plotColors{h});
    hold on;
end


% Adjust labels and title
xlabel('Average Z-score');
ylabel('Trial Counts');
title('Histogram of Response 4s');

%% Box plot

data1 = Response_1s{1};
data2 = Response_1s{2};
% Box plot
figure(histogram_1s);
subplot(2,3,3);
boxplot([data1, data2], [ones(1,15), 2*ones(1,15)]);
hold on;
% scatter plot
for b = 1:2
    figure(histogram_1s);
    subplot(2,3,3);
    scatter(b*ones(size(Response_1s{b})), Response_1s{b}, plotColors{b}, 'filled');
    hold on;
end

data1 = Response_4s{1};
data2 = Response_4s{2};
% Box plot
figure(histogram_1s);
subplot(2,3,6);
boxplot([data1, data2], [ones(1,15), 2*ones(1,15)]);
hold on;
% scatter plot
for b = 1:2
    figure(histogram_1s);
    subplot(2,3,6);
    scatter(b*ones(size(Response_4s{b})), Response_4s{b}, plotColors{b}, 'filled');
    hold on;
end


%% test significance

% initial response significance
data00 = zeros(size(Response_1s_m{1}));
data01 = Response_1s_m{1};
data04 = Response_4s_m{1};
[h, pt01, ci, stats] = ttest2(data00, data01);
[h, pt04, ksstat] = kstest2(data00, data04);
% fprintf('Response t test 1s p= %.4f\n', pt01);
% fprintf('Response t test 4s p= %.4f\n', pt04);


% 1s significance
data1 = Response_1s{1};
data2 = Response_1s{2};
[h, pt1, ci, stats] = ttest2(data1, data2);
[h, pk1, ksstat] = kstest2(data1, data2);

fprintf('  95% Confidence Interval for the Difference in Means: [%.4f, %.4f]\n', ci(1), ci(2));
figure(histogram_1s);
subplot(2,3,3);
R2_text = ['p(t)=' num2str(pt1) newline 'p(KS)=' num2str(pk1)];
yl=ylim;
text(1.5, yl(2)*0.95, R2_text, 'FontSize', 12,'Color','k');

% 4s significance
data3 = Response_4s{1};
data4 = Response_4s{2};
[h, pt4, ci, stats] = ttest2(data3, data4);
[h, pk4, ksstat] = kstest2(data3, data4);

fprintf('  95% Confidence Interval for the Difference in Means: [%.4f, %.4f]\n', ci(1), ci(2));
figure(histogram_1s);
subplot(2,3,6)
R2_text = ['p(t)=' num2str(pt4) newline 'p(KS)=' num2str(pk4)];
yl=ylim;
text(1.5, yl(2)*0.95, R2_text, 'FontSize', 12,'Color','k');

sig = [pt1,pt4];

% fprintf('KS test 1s p= %.4f\n', pk1);
% fprintf('KS test 4s p= %.4f\n', pk4);

end