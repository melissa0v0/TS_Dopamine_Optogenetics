%% Group name and path

close all
mice_num = 12;
group = 'TAF';
day = 'd1';
if strcmp(group,'ADF')
    start_num = 0;
else
    start_num = 6;
end

%% run through raw data and save

Resp1 = zeros(2,15,mice_num);
Resp4 = zeros(2,15,mice_num);
Trace_avg = cell(1,mice_num);
Resp1_ = cell(1,mice_num);
Resp4_ = cell(1,mice_num);
Trace_avg_ = cell(1,mice_num);
Sig = cell(1,mice_num);
for i = 1:mice_num
    mouse_name = [group, num2str(i+start_num)];
    [resp1, resp4, trace_avg, resp1_, resp4_, trace_avg_, sig] = TAFADF_plot_single(mouse_name,day);
    Resp1(:,:,i) = [resp1{1}; resp1{2}];
    Resp4(:,:,i) = [resp4{1}; resp4{2}];
    Trace_avg{i} = trace_avg;
    Sig{i} = sig;

    Resp1_{i} = [resp1_{1}];
    Resp4_{i} = [resp4_{1}];
    Trace_avg_{i} = trace_avg_;
end

% saving
data_all = {Resp1, Resp4, Trace_avg, Sig, Resp1_, Resp4_, Trace_avg_};
filename = [group '_' day '.mat'];
save(filename, "data_all");


%% Load the .mat file

filename = [group '_' 'd1' '.mat'];
loaded_data = load(filename);

data_all_d1 = loaded_data.data_all;
Resp1_d1 = data_all_d1{1};
Resp4_d1 = data_all_d1{2};
Trace_avg_d1 = data_all_d1{3};
Sig_d1 = data_all_d1{4};
Resp1_d1_ = data_all_d1{5};
Resp4_d1_ = data_all_d1{6};
Trace_all_d1_ = data_all_d1{7};

filename = [group '_' 'd2' '.mat'];
loaded_data = load(filename);

data_all_d2 = loaded_data.data_all;
Resp1_d2 = data_all_d2{1};
Resp4_d2 = data_all_d2{2};
Trace_avg_d2 = data_all_d2{3};
Sig_d2 = data_all_d2{4};
Resp1_d2_ = data_all_d2{5};
Resp4_d2_ = data_all_d2{6};
Trace_all_d2_ = data_all_d2{7};

Traces_all_d1 = cell2mat(Trace_avg_d1);
Traces_all_d1 = reshape(Traces_all_d1,2,10001,12);
Traces_all_d2 = cell2mat(Trace_avg_d2);
Traces_all_d2 = reshape(Traces_all_d2,2,10001,12);

%% Truncate resp to 20

cut_thres = 20;

Resp1_d1_cut = NaN(12,cut_thres);
Resp4_d1_cut = NaN(12,cut_thres);
for i = 1:mice_num
    tmp = Resp1_d1_{i};
    if numel(tmp) < cut_thres
        tmp(end+1:cut_thres) = NaN;
    elseif numel(tmp) > cut_thres
        tmp = tmp(1:cut_thres);
    end
    Resp1_d1_cut(i,:) = tmp;

    tmp = Resp4_d1_{i};
    if numel(tmp) < cut_thres
        tmp(end+1:cut_thres) = NaN;
    elseif numel(tmp) > cut_thres
        tmp = tmp(1:cut_thres);
    end
    Resp4_d1_cut(i,:) = tmp;
end


Resp1_d2_cut = NaN(12,cut_thres);
Resp4_d2_cut = NaN(12,cut_thres);
for i = 1:mice_num
    tmp = Resp1_d2_{i};
    if numel(tmp) < cut_thres
        tmp(end+1:cut_thres) = NaN;
    elseif numel(tmp) > cut_thres
        tmp = tmp(1:cut_thres);
    end
    Resp1_d2_cut(i,:) = tmp;

    tmp = Resp4_d2_{i};
    if numel(tmp) < cut_thres
        tmp(end+1:cut_thres) = NaN;
    elseif numel(tmp) > cut_thres
        tmp = tmp(1:cut_thres);
    end
    Resp4_d2_cut(i,:) = tmp;
end


%% Truncate resp to 15 

cut_thres = 15;

Resp1_d1_cut = NaN(12,cut_thres);
Resp4_d1_cut = NaN(12,cut_thres);
for i = 1:mice_num
    tmp = Resp1_d1_{i};
    if numel(tmp) < cut_thres
        tmp(end+1:cut_thres) = NaN;
    elseif numel(tmp) > cut_thres
        tmp = tmp(end-cut_thres+1:end);
    end
    Resp1_d1_cut(i,:) = tmp;

    tmp = Resp4_d1_{i};
    if numel(tmp) < cut_thres
        tmp(end+1:cut_thres) = NaN;
    elseif numel(tmp) > cut_thres
        tmp = tmp(end-cut_thres+1:end);
    end
    Resp4_d1_cut(i,:) = tmp;
end


Resp1_d2_cut = NaN(12,cut_thres);
Resp4_d2_cut = NaN(12,cut_thres);
for i = 1:mice_num
    tmp = Resp1_d2_{i};
    if numel(tmp) < cut_thres
        tmp(end+1:cut_thres) = NaN;
    elseif numel(tmp) > cut_thres
        tmp = tmp(end-cut_thres+1:end);
    end
    Resp1_d2_cut(i,:) = tmp;

    tmp = Resp4_d2_{i};
    if numel(tmp) < cut_thres
        tmp(end+1:cut_thres) = NaN;
    elseif numel(tmp) > cut_thres
        tmp = tmp(end-cut_thres+1:end);
    end
    Resp4_d2_cut(i,:) = tmp;
end


%% Truncate trace to 20

trace_init_d1 = NaN(12,20,10001);
for i = 1:mice_num
    tmp = Trace_all_d1_{i};
    if size(tmp,1) < 20
        tmp(end+1:20,:) = NaN;
    elseif size(tmp,1) > 20
        tmp = tmp(1:20,:);
    end
    trace_init_d1(i,:,:) = tmp;
end


trace_init_d2 = NaN(12,20,10001);
for i = 1:mice_num
    tmp = Trace_all_d2_{i};
    if size(tmp,1) < 20
        tmp(end+1:20,:) = NaN;
    elseif size(tmp,1) > 20
        tmp = tmp(1:20,:);
    end
    trace_init_d2(i,:,:) = tmp;
end

%% Truncate trace to 15

cut_thres = 15;

trace_init_d1 = NaN(12,cut_thres,10001);
for i = 1:mice_num
    tmp = Trace_all_d1_{i};
    if size(tmp,1) < cut_thres
        tmp(end+1:cut_thres,:) = NaN;
    elseif size(tmp,1) > cut_thres
        tmp = tmp(end-cut_thres+1:end,:);
    end
    trace_init_d1(i,:,:) = tmp;
end


trace_init_d2 = NaN(12,cut_thres,10001);
for i = 1:mice_num
    tmp = Trace_all_d2_{i};
    if size(tmp,1) < cut_thres
        tmp(end+1:cut_thres,:) = NaN;
    elseif size(tmp,1) > cut_thres
        tmp = tmp(end-cut_thres+1:end,:);
    end
    trace_init_d2(i,:,:) = tmp;
end

%% Get 1s response for block1 & block2n

resp1_miceAVG_d1i = squeeze(mean(trace_init_d1(:,:,2000:3000),[2,3]))';
resp1_miceAVG_d2i = squeeze(mean(trace_init_d2(:,:,2000:3000),[2,3]))';

resp1_miceAVG_d1n = squeeze(mean(Traces_all_d1(1,2000:3000,:),2))';
resp1_miceAVG_d2n = squeeze(mean(Traces_all_d2(1,2000:3000,:),2))';


%% plot paired response comparision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % use block2 n 0.5s resp
% resp1_miceAVG_d1 = squeeze(mean(Resp1_d1,2));
% resp1_miceAVG_d1n = resp1_miceAVG_d1(1,:);
% resp1_miceAVG_d2 = squeeze(mean(Resp1_d2,2));
% resp1_miceAVG_d2n = resp1_miceAVG_d2(1,:);
% resp1_miceAVG_d1d2n = [resp1_miceAVG_d1n; resp1_miceAVG_d2n];
% 
% 
% % use block1 0.5s resp
% resp1_miceAVG_d1d2n = [mean(Resp1_d1_cut,2,"omitnan")'; mean(Resp1_d2_cut,2,"omitnan")'];
% 
% % use block1 1s resp
% resp1_miceAVG_d1d2n = [resp1_miceAVG_d1i; resp1_miceAVG_d2i];

% use block2n 1s resp
resp1_miceAVG_d1d2n = [resp1_miceAVG_d1n; resp1_miceAVG_d2n];


RESP = figure();
set(gcf, 'Position', [50,50, 1000, 400]);
colors = [repmat('r',1,6), repmat('k',1,6)];

for i = 1:mice_num/2
    figure(RESP);
    subplot(1,4,2);
    plot(1:2,resp1_miceAVG_d1d2n(:,i),"Color",[0.5 0 0],'LineWidth',3); hold on
    scatter(1*ones(1,1), resp1_miceAVG_d1d2n(1,i), 'k', 'filled');hold on
    scatter(2*ones(1,1), resp1_miceAVG_d1d2n(2,i), 'r', 'filled');hold on
end
xlim([0.5,2.5]);
if strcmp(group,'ADF')
    ylim([-1,1]);
else
    ylim([-0.5,2.5]);
end
xticks([1, 2]);
xticklabels({'Day1', 'Day2'});
title('1s Response ChrimsonR');
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);

for i = mice_num/2+1:mice_num
    figure(RESP);
    subplot(1,4,1);
    plot(1:2,resp1_miceAVG_d1d2n(:,i),"Color",[0 0 0.4],'LineWidth',3); hold on
    scatter(1*ones(1,1), resp1_miceAVG_d1d2n(1,i), 'k', 'filled');hold on
    scatter(2*ones(1,1), resp1_miceAVG_d1d2n(2,i), 'b', 'filled');hold on
end
xlim([0.5,2.5]);
if strcmp(group,'ADF')
    ylim([-1,1]);
else
    ylim([-0.5,2.5]);
end
title('Ct');
xticks([1, 2]);
xticklabels({'Day1', 'Day2'});
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);

% Test signigicance for opto
data1 = resp1_miceAVG_d1d2n(1,1:6);
data2 = resp1_miceAVG_d1d2n(2,1:6);
[h, pt_ChR, ci, stats_chr] = ttest(data1, data2);
% Test signigicance for control
data1 = resp1_miceAVG_d1d2n(1,7:end);
data2 = resp1_miceAVG_d1d2n(2,7:end);
[h, pt_Ct, ci, stats_ct] = ttest(data1, data2);
% R2_text = ['p(ChrimsonR)=' num2str(pt_ChR) newline 'p(Control)=' num2str(pt_Ct)];
% yl=ylim;
% text(1, yl(2)*0.95, R2_text, 'FontSize', 12,'Color','k');
plot_sig(pt_Ct,-stats_ct.tstat, RESP, 4,1);
plot_sig(pt_ChR,-stats_chr.tstat, RESP, 4,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4s plots
resp4_miceAVG_d1 = squeeze(mean(Resp4_d1,2));
resp4_miceAVG_d1n = resp4_miceAVG_d1(1,:);
resp4_miceAVG_d2 = squeeze(mean(Resp4_d2,2));
resp4_miceAVG_d2n = resp4_miceAVG_d2(1,:);
resp4_miceAVG_d1d2n = [resp4_miceAVG_d1n;resp4_miceAVG_d2n];

% % use block 1 instead
% resp4_miceAVG_d1d2n = [mean(Resp4_d1_cut,2,"omitnan")';mean(Resp4_d2_cut,2,"omitnan")'];

for i = 1:mice_num/2
    figure(RESP);
    subplot(1,4,4);
    plot(1:2,resp4_miceAVG_d1d2n(:,i),"Color",[0.5 0 0],'LineWidth',3); hold on
    scatter(1*ones(1,1), resp4_miceAVG_d1d2n(1,i), 'k', 'filled');hold on
    scatter(2*ones(1,1), resp4_miceAVG_d1d2n(2,i), 'r', 'filled');hold on
end
xlim([0.5,2.5]);
if strcmp(group,'ADF')
    ylim([-0.8,0.6]);
else
    ylim([-0.5,2]);
end
xticks([1, 2]);
xticklabels({'Day1', 'Day2'});
title('4s Response ChrimsonR');
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1); hold on

for i = mice_num/2+1:mice_num
    figure(RESP);
    subplot(1,4,3);
    plot(1:2,resp4_miceAVG_d1d2n(:,i),"Color",[0 0 0.4],'LineWidth',3); hold on
    scatter(1*ones(1,1), resp4_miceAVG_d1d2n(1,i), 'k', 'filled');hold on
    scatter(2*ones(1,1), resp4_miceAVG_d1d2n(2,i), 'b', 'filled');hold on
end
xlim([0.5,2.5]);
if strcmp(group,'ADF')
    ylim([-0.8,0.6]);
else
    ylim([-0.5,2]);
end
title('Control');
xticks([1, 2]);
xticklabels({'Day1', 'Day2'});
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1); hold on

% Test signigicance for opto
data1 = resp4_miceAVG_d1d2n(1,1:6);
data2 = resp4_miceAVG_d1d2n(2,1:6);
[h, pt_ChR, ci, stats_chr] = ttest(data1, data2);
% Test signigicance for control
data1 = resp4_miceAVG_d1d2n(1,7:end);
data2 = resp4_miceAVG_d1d2n(2,7:end);
[h, pt_Ct, ci, stats_ct] = ttest(data1, data2);

% R2_text = ['p(ChrimsonR)=' num2str(pt_ChR) newline 'p(Control)=' num2str(pt_Ct)];
% yl=ylim;
% text(1, yl(2)*0.95, R2_text, 'FontSize', 12,'Color','k');
plot_sig(pt_Ct, -stats_ct.tstat, RESP, 4,3);
plot_sig(pt_ChR, -stats_chr.tstat, RESP, 4,4);


%% Test individual sig for cross day comparision
Sig_cross_day = zeros(2,12);
for i = 1:mice_num
    data1 = squeeze(Resp1_d1(1,:,i));
    data2 = squeeze(Resp1_d2(1,:,i));
    [h, pt_ct, ci, stats] = ttest(data1, data2);
    Sig_cross_day(1,i) = pt_ct;

    data1 = squeeze(Resp4_d1(1,:,i));
    data2 = squeeze(Resp4_d2(1,:,i));
    [h, pt_chr, ci, stats] = ttest(data1, data2);
    Sig_cross_day(2,i) = pt_chr;
end

%% plot resp difference
threshold = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0.5s plots
diff1 = resp1_miceAVG_d1d2n(2,:) - resp1_miceAVG_d1d2n(1,:);
DIFF = figure();
subplot(1,2,1)
% boxplot([data1, data2], [ones(1,6), 2*ones(1,6)]);hold on
for i = 1:mice_num/2
    if Sig_cross_day(1,i) <=threshold
        scatter(2*ones(6,1), diff1(i), 'filled','MarkerFaceColor',[0.5 0 0]);hold on
    else
        scatter(2*ones(6,1), diff1(i), 'o', 'MarkerEdgeColor',[0.5 0 0]);hold on
    end
end
for i = mice_num/2+1:mice_num
    if Sig_cross_day(1,i) <= threshold
        scatter(1*ones(6,1), diff1(i), 'filled','MarkerFaceColor',[0 0 0.4]);hold on
    else
        scatter(1*ones(6,1), diff1(i), 'o', 'MarkerEdgeColor',[0 0 0.4]);hold on
    end
end
xlim([0.5,2.5]);
if strcmp(group,'ADF')
    ylim([-1,1]);
else
    ylim([-1,1]);
end
title('1s Response diff')
xticks([1, 2]);
xticklabels({'Control', 'ChrimsonR'});
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);

[h, pt_diff05, ci, stats] = ttest2(diff1(1:6), diff1(7:12));
% R2_text = ['p=' num2str(pt_diff05)];
% yl=ylim;
% text(1, yl(2)*0.9, R2_text, 'FontSize', 12,'Color','k');
plot_sig(pt_diff05, stats.tstat, DIFF, 2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4s plots
diff4 = resp4_miceAVG_d1d2n(2,:) - resp4_miceAVG_d1d2n(1,:);
figure(DIFF);
subplot(1,2,2)
boxplot([diff4(7:12), diff4(1:6)], [ones(1,6), 2*ones(1,6)]);hold on
% boxplot([data1, data2], [ones(1,6), 2*ones(1,6)]);hold on
for i = 1:mice_num/2
    if Sig_cross_day(2,i) <=threshold
        scatter(2*ones(6,1), diff4(i), 'filled','MarkerFaceColor',[0.5 0 0]);hold on
    else
        scatter(2*ones(6,1), diff4(i), 'o', 'MarkerEdgeColor',[0.5 0 0]);hold on
    end
end
for i = mice_num/2+1:mice_num
    if Sig_cross_day(2,i) <= threshold
        scatter(1*ones(6,1), diff4(i), 'filled','MarkerFaceColor',[0 0 0.4]);hold on
    else
        scatter(1*ones(6,1), diff4(i),'o', 'MarkerEdgeColor',[0 0 0.4]);hold on
    end
end
xlim([0.5,2.5]);
if strcmp(group,'ADF')
    ylim([-0.8,0.6]);
else
    ylim([-0.8,0.8]);
end
title('4s Response diff')
xticks([1, 2]);
xticklabels({'Control', 'ChrimsonR'});
plot([-0.5, 2.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);

[h, pt_diff4, ci, stats] = ttest2(diff4(1:6), diff4(7:12));
% R2_text = ['p=' num2str(pt_diff4)];
% yl=ylim;
% text(1, yl(2)*0.9, R2_text, 'FontSize', 12,'Color','k');
plot_sig(pt_diff4, stats.tstat, DIFF, 2,2);

%% plot traces block2 no stim

plotWin = [-2000:8000];
legend_name={'Day1','Day2'};
xrange = [-500,1000];
yrange = [-0.6,3];

TRACE = figure();
set(gcf, 'Position', [50,0, 800, 800]);
set(gcf,'color','w');
data_chr_n1 = squeeze(Traces_all_d1(1,:,1:6))';
m_plot = mean(data_chr_n1,1);  % normal trial in chr mice D1
s_plot = std(data_chr_n1)/sqrt(6);
subplot(3,2,2)
rectangle('Position', [0, -3, 4000, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
plot([-3000,8000], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
% rectangle('Position', [1000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
% rectangle('Position', [2000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
% rectangle('Position', [3000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
errorbar_patch(plotWin,m_plot,s_plot,[0.3 0 0],false);hold on
data_chr_n2 = squeeze(Traces_all_d2(1,:,1:6))';
m_plot = mean(data_chr_n2,1);  % normal trial in chr mice D2
s_plot = std(data_chr_n2)/sqrt(6);
errorbar_patch(plotWin,m_plot,s_plot,[0.75 0.3 0.3],false);
title("ChrimsonR")
ylim([-0.5,3]);
xlabel("time(s)");
ylabel("z-score");
h=gca;
h.XTickLabel = {(xrange(1)/1000):(xrange(end)/1000)};
% legend(legend_name)
xlim(xrange);
ylim(yrange)

data_ct_n1 = squeeze(Traces_all_d1(1,:,7:12))';
m_plot = mean(data_ct_n1,1);  % normal trial in ct mice
s_plot = std(data_ct_n1)/sqrt(6);
subplot(3,2,1)
rectangle('Position', [0, -3, 4000, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
plot([-3000,8000], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
% rectangle('Position', [1000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
% rectangle('Position', [2000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
% rectangle('Position', [3000, -3, 500, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
errorbar_patch(plotWin,m_plot,s_plot,[0 0 0.2],false);hold on
data_ct_n2 = squeeze(Traces_all_d2(1,:,7:12))';
m_plot = mean(data_ct_n2,1);  % opto trial in ct mice
s_plot = std(data_ct_n2)/sqrt(6);
errorbar_patch(plotWin,m_plot,s_plot,[0.5 0.5 0.75],false)
title("Control")
ylim([-0.5,3]);
xlabel("time(s)");
ylabel("z-score");
h=gca;
h.XTickLabel = {(xrange(1)/1000):(xrange(end)/1000)};
xlim(xrange);
ylim(yrange);

%% sort according to 4s response d1

% %%%sort according to 4s response diff between d1&d2
% order_chr = [3,5,6,2,4,1];
% order_ct = [8,9,7,12,11,10];
sort_tmp_ct = squeeze(resp4_miceAVG_d1(1,7:12));
sort_tmp_chr = squeeze(resp4_miceAVG_d1(1,1:6));
[sort_tmp_ct,order_ct] = sort(sort_tmp_ct,'descend');
[sort_tmp_chr,order_chr] = sort(sort_tmp_chr,'descend');
data_chr_unsorted = [data_chr_n1; data_chr_n2];
data_ct_unsorted = [data_ct_n1; data_ct_n2];
data_chr = zeros(size(data_chr_unsorted));
data_ct = zeros(size(data_ct_unsorted));
for i = 1:6
    data_ct(i,:) = data_ct_unsorted(order_ct(i),:);
    data_ct(i+6,:) = data_ct_unsorted(order_ct(i)+6,:);
    data_chr(i,:) = data_chr_unsorted(order_chr(i),:);
    data_chr(i+6,:) = data_chr_unsorted(order_chr(i)+6,:);
end

%% sorted heat map block2 no stim

heatmap_range = [-2,2];
xrange = [15,30];

figure(TRACE)
subplot(3,2,4)
trialNum = size(data_chr,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(data_chr(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time');
ylabel('trials')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
xlim(xrange);
title("ChrimsonR orig")
line([20,20], [0, size(data_chr,1)+1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1,'LineStyle', '--');
line([0,80], [6.5,6.5], 'Color', [0.8,0.8,0.8], 'LineWidth', 1);


figure(TRACE)
subplot(3,2,3)
trialNum = size(data_ct,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(data_ct(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time');
ylabel('trials')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
xlim(xrange);
title("Control orig")
line([20,20], [0, size(data_chr,1)+1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1,'LineStyle', '--');
line([0,80], [6.5,6.5], 'Color', [0.8,0.8,0.8], 'LineWidth', 1);


% differencce heat map chr mice
figure(TRACE)
subplot(3,2,6)
trace_chr_diff = data_chr(7:12,:) - data_chr(1:6,:);
trialNum = size(trace_chr_diff,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(trace_chr_diff(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time');
ylabel('animals')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
xlim(xrange);
title("ChrimsonR diff")
line([20,20], [0, size(data_ct_diff,1)+1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1,'LineStyle', '--');

% differencce heat map ct mice
figure(TRACE)
subplot(3,2,5)
trace_ct_diff = data_ct(7:12,:) - data_ct(1:6,:);
trialNum = size(trace_ct_diff,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(trace_ct_diff(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time');
ylabel('animals')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
xlim(xrange);
title("Control diff")
line([20,20], [0, size(data_ct_diff,1)+1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1,'LineStyle', '--');

%% Plot traces block1

TRACE_INIT = figure();

trace_init_d1_miceavg = squeeze(mean(trace_init_d1,2,"omitnan"));
trace_init_d2_miceavg = squeeze(mean(trace_init_d2,2,"omitnan"));
plotWin = [-2000:8000];
legend_name={'Day1','Day2'};

subplot(3,2,1)
m_plot = mean(trace_init_d1_miceavg(7:12,:),1);  % init trial in ct mice d1
s_plot = std(trace_init_d1_miceavg(7:12,:))/sqrt(6);
rectangle('Position', [0, -3, 4000, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
errorbar_patch(plotWin,m_plot,s_plot,'k',false);hold on
m_plot = mean(trace_init_d2_miceavg(7:12,:),1);  % init trial in ct mice d2
s_plot = std(trace_init_d2_miceavg(7:12,:))/sqrt(6);
errorbar_patch(plotWin,m_plot,s_plot,'b',false)
title("Control")
ylim([-1,4]);
xlabel("time");
ylabel("z-score");
xlim([-1000,1000]);

subplot(3,2,2)
m_plot = mean(trace_init_d1_miceavg(1:6,:),1);  % init trial in chr mice d1
s_plot = std(trace_init_d1_miceavg(1:6,:))/sqrt(6);
rectangle('Position', [0, -3, 4000, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
errorbar_patch(plotWin,m_plot,s_plot,'k',false);hold on
m_plot = mean(trace_init_d2_miceavg(1:6,:),1);  % init trial in chr mice d1
s_plot = std(trace_init_d2_miceavg(1:6,:))/sqrt(6);
errorbar_patch(plotWin,m_plot,s_plot,'b',false);hold on
title("ChrimsonR")
ylim([-1,4]);
xlabel("time");
ylabel("z-score");
xlim([-1000,1000]);

%% Plot trace diff block1

TRACE_DIFF = figure();
subplot(1,3,1)
rectangle('Position', [0, -3, 4000, 10], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on

trace_diff_chr = trace_init_d2_miceavg(1:6,:) - trace_init_d1_miceavg(1:6,:);
trace_diff_ct = trace_init_d2_miceavg(7:12,:) - trace_init_d1_miceavg(7:12,:);
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

%% Plot heatmap block1

% control day1 & day2
subplot(3,2,5)
trace_init_ct = [trace_init_d1_miceavg(7:12,:); trace_init_d2_miceavg(7:12,:)];
trialNum = size(trace_init_ct,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(trace_init_ct(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time');
ylabel('trials')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
title("trial average day1")
xlim([10,30]);

% chr day1 & day2
subplot(3,2,6)
trace_init_chr = [trace_init_d1_miceavg(1:6,:); trace_init_d2_miceavg(1:6,:)];
trialNum = size(trace_init_chr,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(trace_init_chr(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time');
ylabel('trials')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
title("trial average day2")
xlim([10,30]);

%diff heat map
subplot(3,2,3)
trace_ct_diff = trace_init_d2_miceavg(7:12,:) - trace_init_d1_miceavg(7:12,:);
trialNum = size(trace_ct_diff,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(trace_ct_diff(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time');
ylabel('animals')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
title("Control")
xlim([10,30]);

subplot(3,2,4)
trace_chr_diff = trace_init_d2_miceavg(1:6,:) - trace_init_d1_miceavg(1:6,:);
trialNum = size(trace_chr_diff,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(trace_chr_diff(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,heatmap_range);
colormap yellowblue
xlabel('time');
ylabel('animals')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
title("ChrimsonR")
xlim([10,30]);





%% local functions

function plot_sig(pt_ChR, t_stat, RESP,a,b)
figure(RESP);
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