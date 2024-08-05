%% Group name and path

close all
mice_num = 6;

Trace_init = zeros(10001,mice_num);
Trace_opto = zeros(3,10001,mice_num);
Resp_init = zeros(1,mice_num);
Resp_opto = zeros(3,mice_num);
for i = 1:mice_num
    [trace, resp, ste_resp] = SingleAnimal_natural_resp(i);
    % [trace_o, resp_o, ste_resp_o] = SingleAnimal_opto_3frequency_500ms(mouse_name);
    [trace_o, resp_o, ste_resp_o] = SingleAnimal_opto_3frequency_4s(i);
    Trace_init(:,i) = trace;
    Trace_opto(:,:,i) = trace_o;
    Resp_init(:,i) = mean(resp);
    Resp_opto(:,i) = resp_o;
end

%% Plot traces

TRACE = figure();
set(gcf,'color','w')
plotWin = [-2000:8000];
xrange = [-500,4500];
% xrange = [-1000,1000];

legend_name = {'10Hz','20Hz','40Hz','natural'};
plot_colors = [1 0 0; % Red
              0.8 0 0; % Darker red
              0.6 0 0]; % Darkest red

% subplot(2,2,1)
% rectangle('Position', [0, -1, 4000, 7], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
% data_ct = Trace_init(:,1:3)';  % control mice stimulation
% m_plot = mean(data_ct,1);  
% s_plot = std(data_ct)/sqrt(3);
% errorbar_patch(plotWin,m_plot,s_plot,'k',false);hold on
subplot(2,2,3)
rectangle('Position', [0, -1, 500, 7], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
for i = 1:3
    data_ct = squeeze(Trace_opto(i,:,1:3))';  % control mice stimulation
    m_plot = mean(data_ct,1);  
    s_plot = std(data_ct)/sqrt(3);
    figure(TRACE)
    subplot(2,2,3)
    errorbar_patch(plotWin,m_plot,s_plot,plot_colors(i,:),false);hold on
end
title("control mice")
xlabel('time(ms)');
ylabel('z-score');
xlim(xrange);
ylim([-1,6])

subplot(2,2,2)
rectangle('Position', [0, -1, 4000, 7], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
plotWin = [-2000:8000];
data_chr = Trace_init(:,4:6)';  % ChrmsonR mice stimulation
m_plot = mean(data_chr,1);  
s_plot = std(data_chr)/sqrt(3);
errorbar_patch(plotWin,m_plot,s_plot,'k',false);hold on
xlim(xrange)
ylim([-1,6])
title('first ten natural resp')

subplot(2,2,4)
rectangle('Position', [0, -1, 500, 7], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
for i = 1:3
    data_chr = squeeze(Trace_opto(i,:,4:6))';  % ChrmsonR mice stimulation
    m_plot = mean(data_chr,1);  
    s_plot = std(data_chr)/sqrt(3);
    figure(TRACE)
    subplot(2,2,4)
    errorbar_patch(plotWin,m_plot,s_plot,plot_colors(i,:),false);hold on
end

legend(legend_name)
title("ChrimsonR mice")
xlabel('time(ms)');
ylabel('z-score');
xlim(xrange);
ylim([-1,6])

%% all mice natural resp

figure(TRACE)
subplot(2,2,1)
rectangle('Position', [0, -1, 4000, 7], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
data_ctchr = Trace_init(:,:)';  % control mice stimulation
m_plot = mean(data_ctchr,1);  
s_plot = std(data_ctchr)/sqrt(6);
errorbar_patch(plotWin,m_plot,s_plot,'k',false);hold on
title('all mice first natural resp')
xlim(xrange)
ylim([-1,6])

%% plot paired resp

% 4s plots
RESP = figure();
plot_paired_resp(RESP, [Resp_init;Resp_opto], mice_num, 'stimulation',1);

%% Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%% paired response plot
function [] = plot_paired_resp(RESP, resp1_miceAVG_d1d2n, mice_num, day,plot_num)
colors = [repmat('k',1,mice_num/2) repmat('r',1,mice_num/2)];

for i = 1:mice_num/2
    figure(RESP);
    subplot(1,2,plot_num);
    plot(1:4,resp1_miceAVG_d1d2n(:,i),"Color",colors(i),'LineWidth',3); hold on
end
xlim([0.5,4.5]);
ylim([-1,6]);
xticks([1,2,3,4]);
xticklabels({'Natural', '10Hz','20Hz','40Hz'});
title(['Control']);

for i = mice_num/2+1:mice_num
    figure(RESP);
    subplot(1,2,plot_num+1);
    plot(1:4,resp1_miceAVG_d1d2n(:,i),"Color",colors(i),'LineWidth',3); hold on
end
xlim([0.5,4.5]);
ylim([-1,6]);
xticks([1,2,3,4]);
xticklabels({'Natural', '10Hz','20Hz','40Hz'});
title('ChrimsonR');

end
