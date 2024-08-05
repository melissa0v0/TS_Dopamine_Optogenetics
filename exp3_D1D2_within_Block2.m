%% Response for all mice [trial type x num of trials x mice]

close all
mice_num = 12;
group = 'ADF';
day = 'd1';
if strcmp(group,'ADF')      % D2 mice
    start_num = 0;
elseif strcmp(group,'TAF')   % D1 mice
    start_num = 6;
end

%% Load the .mat file

filename = [group '_' day '.mat'];
loaded_data = load(filename);

data_all = loaded_data.data_all;
Resp1 = data_all{1};
Resp4 = data_all{2};
Trace_avg = data_all{3};
Sig = data_all{4};

%% Run through data and save
% Only need to be run at the first time
% Used to generate and save intermediate files
% including TAF_d1.mat, TAF_d2.mat, ADF_d1.mat, ADF_d2.mat
% Saving these files make further analysis faster

Resp1 = zeros(2,15,mice_num);
Resp4 = zeros(2,15,mice_num);
Trace_avg = cell(1,mice_num);
Sig = cell(1,mice_num);
for i = 1:mice_num
    mouse_name = [group, num2str(i+start_num)];
    [resp1, resp4, trace_avg, sig] = TAFADF_plot_single(mouse_name,day);
    Resp1(:,:,i) = [resp1{1}; resp1{2}];
    Resp4(:,:,i) = [resp4{1}; resp4{2}];
    Trace_avg{i} = trace_avg;
    Sig{i} = sig;
end

% saving
data_all = {Resp1, Resp4, Trace_avg, Sig};
filename = [group '_' day '.mat'];
save(filename, "data_all");


%% plot paired response comparision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1s plots;
averages = cellfun(@(x) mean(x(:,2000:3000), 2), Trace_avg, 'UniformOutput', false);
resp1_miceAVG = cat(3, averages{:});
% resp1_miceAVG = mean(Resp1,2);

RESP = figure();
set(gcf, 'Position', [50,50, 1000, 400]);
set(gcf,'color','w');
colors = [repmat('r',1,6), repmat('b',1,6)];

for i = 1:mice_num/2
    figure(RESP);
    subplot(1,4,2);
    plot(1:2,resp1_miceAVG(:,:,i),"Color",[0.8,0.8,0.8],'LineWidth',1); hold on
    scatter(1*ones(1,1), resp1_miceAVG(1,:,i), 'k', 'filled');hold on
    scatter(2*ones(1,1), resp1_miceAVG(2,:,i), 'r', 'filled');hold on
end
xlim([0.5,2.5]);
if strcmp(group,'ADF')
    ylim([-1,1]);
else
    ylim([-0.5,2.5]);
end
title([group ' ' day ' 1s Response ChrimsonR']);
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
xticks([1, 2]);
xticklabels({'no stim', 'stim'});

for i = mice_num/2+1:mice_num
    figure(RESP);
    subplot(1,4,1);
    plot(1:2,resp1_miceAVG(:,:,i),"Color",[0.8,0.8,0.8],'LineWidth',1); hold on
    scatter(1*ones(1,1), resp1_miceAVG(1,:,i), 'k', 'filled');hold on
    scatter(2*ones(1,1), resp1_miceAVG(2,:,i), 'b', 'filled');hold on
end
xlim([0.5,2.5]);
if strcmp(group,'ADF')
    ylim([-1,1]);
else
    ylim([-0.5,2.5]);
end
title('Ct');
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
xticks([1, 2]);
xticklabels({'no stim', 'stim'});

% Test signigicance for opto
data1 = resp1_miceAVG(1,1,1:6);
data2 = resp1_miceAVG(2,1,1:6);
[h, pt_ChR, ci, stats_chr] = ttest(data1, data2);
% Test signigicance for control
data1 = resp1_miceAVG(1,1,6:end);
data2 = resp1_miceAVG(2,1,6:end);
[h, pt_Ct, ci, stats_ct] = ttest(data1, data2);

R2_text = ['p(ChrimsonR)=' num2str(pt_ChR) newline 'p(Control)=' num2str(pt_Ct)];
yl=ylim;
plot_sig(pt_Ct, -stats_ct.tstat, RESP, 4,1);
plot_sig(pt_ChR,-stats_chr.tstat, RESP, 4,2);
% text(1, yl(2)*0.95, R2_text, 'FontSize', 12,'Color','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4s plots
resp4_miceAVG = mean(Resp4,2);
for i = 1:mice_num/2
    figure(RESP);
    subplot(1,4,4);
    plot(1:2,resp4_miceAVG(:,:,i),"Color",[0.8,0.8,0.8],'LineWidth',1); hold on
    scatter(1*ones(1,1), resp4_miceAVG(1,:,i), 'k', 'filled');hold on
    scatter(2*ones(1,1), resp4_miceAVG(2,:,i), 'r', 'filled');hold on
end
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
if strcmp(group,'ADF')
    ylim([-0.8,0.6]);
else
    ylim([-0.5,1.5]);
end
title([group ' ' day ' 4s Response ChrimsonR']);
xlim([0.5,2.5]);
xticks([1, 2]);
xticklabels({'no stim', 'stim'});
for i = mice_num/2+1:mice_num
    figure(RESP);
    subplot(1,4,3);
    plot(1:2,resp4_miceAVG(:,:,i),"Color",[0.8,0.8,0.8],'LineWidth',1); hold on
    scatter(1*ones(1,1), resp4_miceAVG(1,:,i), 'k', 'filled');hold on
    scatter(2*ones(1,1), resp4_miceAVG(2,:,i), 'b', 'filled');hold on
end
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
if strcmp(group,'ADF')
    ylim([-0.8,0.6]);
else
    ylim([-0.5,1.5]);
end
title('Control');
xlim([0.5,2.5]);
xticks([1, 2]);
xticklabels({'no stim', 'stim'});

% Test signigicance for opto
data1 = resp4_miceAVG(1,1,1:6);
data2 = resp4_miceAVG(2,1,1:6);
[h, pt_ChR, ci, stats_chr] = ttest(data1, data2);
% Test signigicance for control
data1 = resp4_miceAVG(1,1,6:end);
data2 = resp4_miceAVG(2,1,6:end);
[h, pt_Ct, ci, stats_ct] = ttest(data1, data2);

R2_text = ['p(ChrimsonR)=' num2str(pt_ChR) newline 'p(Control)=' num2str(pt_Ct)];
% text(1, yl(2)*0.95, R2_text, 'FontSize', 12,'Color','k');
plot_sig(pt_Ct, -stats_ct.tstat,RESP, 4,3);
plot_sig(pt_ChR, -stats_chr.tstat,RESP, 4,4);


%% plot resp difference
threshold = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1s plots
diff5 = Resp1(2,:,:) - Resp1(1,:,:);
data1 = reshape(mean(diff5(1,:,1:6),2),1,[]);
data2 = reshape(mean(diff5(1,:,7:12),2),1,[]);

diff1 = resp1_miceAVG(2,:) - resp1_miceAVG(1,:);
data1 = diff1(1:6);
data2 = diff1(7:12);

DIFF = figure();
set(gcf, 'Position', [50,50, 600, 400]);
set(gcf,'color','w');
subplot(1,2,1)
% boxplot([data1, data2], [ones(1,6), 2*ones(1,6)]);hold on
boxplot([data2, data1], [ones(1,6), 2*ones(1,6)]);hold on
for i = 1:mice_num/2
    if Sig{i}(1) <=threshold
        scatter(2*ones(size(data1)), data1(i), 'r', 'filled');hold on
    else
        scatter(2*ones(size(data1)), data1(i), 'r', 'o');hold on
    end
end
for i = mice_num/2+1:mice_num
    if Sig{i}(1) <= threshold
        scatter(1*ones(size(data2)), data2(i-6), 'b', 'filled');hold on
    else
        scatter(1*ones(size(data1)), data2(i-6), 'b', 'o');hold on
    end
end
xlim([0.5,2.5]);
if strcmp(group,'ADF')
    ylim([-0.6,0.6]);
else
    ylim([0,1]);
end
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
title('1s Response diff')
xticklabels({'Control', 'ChrimsonR'});
xticks([1,2])

[h, pt_diff05, ci, stats] = ttest2(data1, data2);
% R2_text = ['p=' num2str(pt_diff05)];
% yl=ylim;
% text(1, yl(2)*0.9, R2_text, 'FontSize', 12,'Color','k');
plot_sig(pt_diff05, stats.tstat,DIFF, 2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4s plots
diff4 = Resp4(2,:,:) - Resp4(1,:,:);
data1 = reshape(mean(diff4(1,:,1:6),2),1,[]);
data2 = reshape(mean(diff4(1,:,7:12),2),1,[]);
figure(DIFF)
subplot(1,2,2)
boxplot([data2, data1], [ones(1,6), 2*ones(1,6)]);hold on
for i = 1:mice_num/2
    if Sig{i}(2) <=threshold
        scatter(2*ones(size(data1)), data1(i), 'r', 'filled');hold on
    else
        scatter(2*ones(size(data1)), data1(i), 'r', 'o');hold on
    end
end
for i = mice_num/2+1:mice_num
    if Sig{i}(2) <= threshold
        scatter(1*ones(size(data2)), data2(i-6), 'b', 'filled');hold on
    else
        scatter(1*ones(size(data2)), data2(i-6), 'b', 'o');hold on
    end
end
xlim([0.5,2.5]);
if strcmp(group,'ADF')
    ylim([-0.8,0.6]);
else
    ylim([-0.8,0.8]);
end
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
title('4s Response diff')
xticklabels({'Control', 'ChrimsonR'});
xticks([1,2])


[h, pt_diff4, ci, stats] = ttest2(data1, data2);
% R2_text = ['p=' num2str(pt_diff4)];
% yl=ylim;
% text(1, yl(2)*0.9, R2_text, 'FontSize', 12,'Color','k');
plot_sig(pt_diff4,stats.tstat, DIFF, 2,2);


%% plot traces

Traces_all = cat(3,Trace_avg{:});
plotWin = [-2000:8000];
xrange = [-2000,8000];
xrange = [-500,1000];
% yrange = [-0.6,1.2];
yrange = [-0.5,3];

TRACEE = figure();
set(gcf, 'Position', [50,0, 800, 800]);
set(gcf,'color','w');
trace_chr_n = squeeze(Traces_all(1,:,1:6))';
m_plot = mean(trace_chr_n,1);  % normal trial in chr mice
s_plot = std(trace_chr_n)/sqrt(6);
subplot(3,2,2)
rectangle('Position', [0, -3, 4000, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
plot([-3000,8000], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
% rectangle('Position', [1000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
% rectangle('Position', [2000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
% rectangle('Position', [3000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
errorbar_patch(plotWin,m_plot,s_plot,'k',false);hold on
trace_chr_o = squeeze(Traces_all(2,:,1:6))';
m_plot = mean(trace_chr_o,1);  % opto trial in chr mice
s_plot = std(trace_chr_o)/sqrt(6);
errorbar_patch(plotWin,m_plot,s_plot,'r',false)

title("ChrimsonR")
h=gca;
h.XTickLabel = {(xrange(1)/1000):(xrange(end)/1000)};
ylim([-0.5,3]);
xlabel("time(s)");
ylabel("z-score");
xlim(xrange);
ylim(yrange);

trace_ct_n = squeeze(Traces_all(1,:,7:12))';
m_plot = mean(trace_ct_n,1);  % normal trial in ct mice
s_plot = std(trace_ct_n)/sqrt(6);
subplot(3,2,1)
rectangle('Position', [0, -3, 4000, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
plot([-3000,8000], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
% rectangle('Position', [1000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
% rectangle('Position', [2000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
% rectangle('Position', [3000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
errorbar_patch(plotWin,m_plot,s_plot,'k',false);hold on
trace_ct_o = squeeze(Traces_all(2,:,7:12))';
m_plot = mean(trace_ct_o,1);  % opto trial in ct mice
s_plot = std(trace_ct_o)/sqrt(6);
errorbar_patch(plotWin,m_plot,s_plot,'b',false)

title("Control")
h=gca;
h.XTickLabel = {(xrange(1)/1000):(xrange(end)/1000)};
ylim([-0.5,3]);
xlabel("time(s)");
ylabel("z-score");
xlim(xrange);
ylim(yrange);


%% animal-by-animal for comfirming error bar

aba = figure();
plotColors_tbt = parula(6);
for f = 1:size(trace_ct_o,1)
    figure(aba)
    subplot(1,2,1)
    errorbar_patch(plotWin,trace_ct_o(f,:),zeros(1,size(trace_ct_o,2)),plotColors_tbt(f,:))
end

for f = 1:size(trace_ct_n,1)
    figure(aba)
    subplot(1,2,2)
    errorbar_patch(plotWin,trace_ct_n(f,:),zeros(1,size(trace_ct_n,2)),plotColors_tbt(f,:))
end


%% Plot trace diff

TRACE_DIFF = figure();
subplot(1,3,1)
rectangle('Position', [0, -3, 4000, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on

trace_diff_chr = trace_chr_o - trace_chr_n;
trace_diff_ct = trace_ct_o - trace_ct_n;
m_plot = mean(trace_diff_chr,1);  % diff in chr mice
s_plot = std(trace_diff_chr)/sqrt(6);
errorbar_patch(plotWin,m_plot,s_plot,'m',false);hold on
m_plot = mean(trace_diff_ct,1);  % diff in chr mice
s_plot = std(trace_diff_ct)/sqrt(6);
errorbar_patch(plotWin,m_plot,s_plot,'k',false)
ylim([-2,2])
legend({'ChrimsonR','Control'})
title("stim trial - no stim trial")

figure(TRACE_DIFF);
subplot(1,3,2)
trace_diff = mean(trace_diff_chr,1) - mean(trace_diff_ct,1);
s_plot = zeros(size(trace_diff));
errorbar_patch(plotWin,trace_diff,s_plot,'k',false);hold on
ylim([-2,2])


subplot(1,3,3)
heatmap_range = [-1,1];
trialNum = size(trace_diff,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(trace_diff(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF',heatmap_range);
colormap yellowblue
ylabel('trials')
title("ChrimsonR diff - control diff")
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};

%% Response latency

latency_ct_n = get_latency(trace_ct_n);
latency_ct_o = get_latency(trace_ct_o);
latency_chr_n = get_latency(trace_chr_n);
latency_chr_o = get_latency(trace_chr_o);

latency = figure();
subplot(1,2,1);    % control mice
scatter(1*ones(size(latency_ct_n)), latency_ct_n, 'k', 'filled');hold on
scatter(2*ones(size(latency_ct_o)), latency_ct_o, 'r', 'filled')
xlim([0.5,2.5])
ylim([2000,2200])
title('Control')
xticks([1,2]);
xticklabels({'no stim', 'stim'});
ylabel('Response latency');

subplot(1,2,2);    % Chr mice
scatter(1*ones(size(latency_chr_n)), latency_chr_n, 'k', 'filled');hold on
scatter(2*ones(size(latency_chr_o)), latency_chr_o, 'r', 'filled')
xlim([0.5,2.5])
ylim([2000,2200])
title('ChrimsonR')
xticks([1,2]);
xticklabels({'no stim', 'stim'});
ylabel('Response latency');

%% sort heat map

% sort_tmp_ct = squeeze(resp4_miceAVG(1,:,7:12));
% sort_tmp_chr = squeeze(resp4_miceAVG(1,:,1:6));
sort_tmp_ct = squeeze(resp1_miceAVG(1,:,7:12));
sort_tmp_chr = squeeze(resp1_miceAVG(1,:,1:6));
[sort_tmp_ct,order_ct] = sort(sort_tmp_ct,'descend');
[sort_tmp_chr,order_chr] = sort(sort_tmp_chr,'descend');
data_chr_unsorted = [trace_chr_n; trace_chr_o];
data_ct_unsorted = [trace_ct_n; trace_ct_o];
data_chr = zeros(size(data_chr_unsorted));
data_ct = zeros(size(data_ct_unsorted));
for i = 1:6
    data_ct(i,:) = data_ct_unsorted(order_ct(i),:);
    data_ct(i+6,:) = data_ct_unsorted(order_ct(i)+6,:);
    data_chr(i,:) = data_chr_unsorted(order_chr(i),:);
    data_chr(i+6,:) = data_chr_unsorted(order_chr(i)+6,:);
end

%% heat map
heatmap_range = [-3,3];
plotWin = [-2000:8000];
figure(TRACEE);
xrange = [0,100];
xrange = [15,30];

% heat map chr
subplot(3,2,4)
trialNum = size(data_chr,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(data_chr(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);hold on
colormap yellowblue
xlabel('time');
ylabel('animals')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
xlim(xrange)
title("ChrimsonR orig")
line([20,20], [0, size(data_chr,1)+1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1,'LineStyle', '--');
line([0,80], [6.5,6.5], 'Color', [0.8,0.8,0.8], 'LineWidth', 1);


% heat map ct
subplot(3,2,3)
trialNum = size(data_ct,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(data_ct(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time');
ylabel('animals')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
xlim(xrange)
title("Control orig")
line([20,20], [0, size(data_chr,1)+1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1,'LineStyle', '--');
line([0,80], [6.5,6.5], 'Color', [0.8,0.8,0.8], 'LineWidth', 1);


%diff heat map
% DIFF_HEAT = figure();
subplot(3,2,5)
data_ct_diff = data_ct(7:12,:) - data_ct(1:6,:);
trialNum = size(data_ct_diff,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(data_ct_diff(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time');
ylabel('animals')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
xlim(xrange)
title("Control diff")
line([20,20], [0, size(data_ct_diff,1)+1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1,'LineStyle', '--');

subplot(3,2,6)
data_chr_diff = data_chr(7:12,:) - data_chr(1:6,:);
trialNum = size(data_chr_diff,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(data_chr_diff(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time(s)');
ylabel('animals')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
xlim(xrange)
title("ChrimsonR diff")
line([20,20], [0, size(data_chr_diff,1)+1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1,'LineStyle', '--');


%% figure settings

all_figures = findall(0, 'type', 'figure'); % Get handles of all open figures

for i = 1:length(all_figures)
    figure(all_figures(i)); % Switch to each figure
    axs = findall(gcf, 'type', 'axes'); % Get handles of all axes in the current figure
    
    for j = 1:length(axs)
        ax = axs(j); % Switch to each subplot
        % Modify settings
        box(ax, 'off'); % Turn off box
        set(ax, 'TickDir', 'out'); % Set tick direction to 'out
        set(ax, 'TickLength', 2*(get(ax, 'TickLength'))); % Increase tick length
    end
end


%%
function plot_sig(pt_ChR,t_stat, RESP,a,b)
figure(RESP);
set(gcf,'color','w')
yl=ylim;
if pt_ChR < 0.01
    subplot(1,a,b);
    text(1.4,yl(2)*0.9,['** p=',newline,num2str(pt_ChR),newline, num2str(t_stat)],'FontSize', 12,'Color','k');
elseif pt_ChR < 0.05
    subplot(1,a,b);
    text(1.5,yl(2)*0.9,['* p=',newline,num2str(pt_ChR),newline, num2str(t_stat)],'FontSize', 12,'Color','k');
else
    subplot(1,a,b);
    text(1.35,yl(2)*0.9,['n.s. p=',newline,num2str(pt_ChR),newline, num2str(t_stat)],'FontSize', 12,'Color','k');
end
end

function [latency] = get_latency(trace)
% Calculate the standard deviation for each animal
std_dev = std(trace, 0, 2); % Use '0' as the second argument to compute the standard deviation normalized by N-1
% Find the index where the value first exceeds 3 times the standard deviation
latency = zeros(size(trace, 1), 1);
for i = 1:size(trace, 1)
    exceed_indices = find(trace(i, :) > 3 * std_dev(i), 1);
    if ~isempty(exceed_indices)
        latency(i) = exceed_indices;
    end
end
end