%% Group name and path

close all

mice_num = 6;

Trace = zeros(10001,mice_num);
Trace_m = zeros(10001,mice_num);
Resp = zeros(2,mice_num);
Resp_m = zeros(2,mice_num);
for i = 1:mice_num
    [trace, resp, ste_resp] = SingleAnimal_masking(i, 'n');
    [trace_m, resp_m, ste_resp_m] = SingleAnimal_masking(i, 'm');
    Trace(:,i) = trace;
    Trace_m(:,i) = trace_m;
    Resp(:,i) = [resp, ste_resp];
    Resp_m(:,i) = [resp_m, ste_resp_m];
end

%% Plot traces

TRACE_CT = figure();
set(gcf,'color','w')
plotWin = [-2000:8000];
xrange = [-500,4500];
yrange = [-1,3];
% legend_name = {'without masking','with masking'};
legend_name = {'Control','ChrimsonR'};

subplot(2,2,1)
rectangle('Position', [0, -1, 500, 5], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
data_ct = Trace(:,1:3)';  % control mice without masking
m_plot = mean(data_ct,1);  
s_plot = std(data_ct)/sqrt(3);
errorbar_patch(plotWin,m_plot,s_plot,'k',false);hold on
title("without masking")
ylabel('Control');
xlim(xrange);
ylim(yrange);

subplot(2,2,3)
rectangle('Position', [0, -1, 500, 5], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
data_chr = Trace(:,4:6)';  % chr mice without masking
m_plot = mean(data_chr,1);  
s_plot = std(data_chr)/sqrt(3);
errorbar_patch(plotWin,m_plot,s_plot,'r',false);hold on
ylabel('ChrimsonR');
xlabel('time(ms)');
xlim(xrange);
ylim(yrange);

subplot(2,2,2)
rectangle('Position', [0, -1, 500, 5], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
data_ct_m = Trace_m(:,1:3)';
m_plot = mean(data_ct_m,1);  % control mice with masking
s_plot = std(data_ct_m)/sqrt(3);
errorbar_patch(plotWin,m_plot,s_plot,'k',false);
title("with masking")
legend({"Control"})
xlim(xrange);
ylim(yrange);

subplot(2,2,4)
rectangle('Position', [0, -1, 500, 5], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
data_chr_m = Trace_m(:,4:6)';
m_plot = mean(data_chr_m,1);  % chr mice with masking
s_plot = std(data_chr_m)/sqrt(3);
errorbar_patch(plotWin,m_plot,s_plot,'r',false);
legend({"ChrimsonR"})
xlabel('time(ms)');
xlim(xrange);
ylim(yrange);

%% Test sig

[h,p_ct_m,ci,stat_ct_m] = ttest(mean(data_ct_m,1),0);

resp1_ct_m = mean(data_ct_m(:,2000:3000),2);
[h,p_ct_m,ci,stat_ct_m] = ttest(resp1_ct_m,0)

resp1_ct = mean(data_ct(:,2000:3000),2);
[h,p_ct,ci,stat_ct] = ttest(resp1_ct,0)
