%% Response for all mice [trial type x num of trials x mice]

close all
mice_num = 6;


%% acquire data

Trace = zeros(2,10001,mice_num);
Resp = zeros(2,mice_num);
for i = 1:mice_num
    mouse_name = [group, num2str(i+start_num)];
    [trace, resp, ste_resp] = SingleAnimal_task_DA(mouse_name);
    Trace(:,:,i) = trace;
    Resp(:,i) = resp;
end

%% plot traces

TRACE = figure();
set(gcf, 'Position', [50,50, 800, 300]);
set(gcf,'color','w');
plotWin = [-2000:8000];
xrange = [-1000,5000];

legend_name = {'no stimulation','stimulation'};
plot_colors = {'k','r'};

subplot(1,2,1)
rectangle('Position', [0, -1, 4000, 6], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
plot([-3000, 8000], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);hold on
for i = 1:2
    data_ct = squeeze(Trace(i,:,1:3))';  % control mice 
    m_plot = mean(data_ct,1);  
    s_plot = std(data_ct)/sqrt(3);
    figure(TRACE)
    subplot(1,2,1)
    errorbar_patch(plotWin,m_plot,s_plot,plot_colors{i},false);hold on
end
title("control mice")
h=gca;
h.XTickLabel = {(xrange(1)/1000):(xrange(end)/1000)};
xlabel('time(ms)');
ylabel('z-score');
xlim(xrange);
ylim([-1,3])

subplot(1,2,2)
rectangle('Position', [0, -1, 4000, 6], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
plot([-3000, 8000], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);hold on
for i = 1:2
    data_ct = squeeze(Trace(i,:,4:6))';  % ChrimsonR mice 
    m_plot = mean(data_ct,1);  
    s_plot = std(data_ct)/sqrt(3);
    figure(TRACE)
    subplot(1,2,2)
    errorbar_patch(plotWin,m_plot,s_plot,plot_colors{i},false);hold on
end
title("ChrimsonR mice")
h=gca;
h.XTickLabel = {(xrange(1)/1000):(xrange(end)/1000)};
XTickLabel = {(xrange(1)/1000):(xrange(end)/1000)};
xlabel('time(s)');
ylabel('z-score');
xlim(xrange);
ylim([-1,3])
% legend(legend_name)

%% plot paired resp

% 4s plots
RESP = figure();
plot_paired_resp(RESP, Resp, mice_num, group, 'task',1);

%% Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%% paired response plot
function [] = plot_paired_resp(RESP, resp1_miceAVG_d1d2n, mice_num, group, day,plot_num)
colors = [repmat('k',1,mice_num/2) repmat('r',1,mice_num/2)];

for i = 1:mice_num/2
    figure(RESP);
    subplot(1,2,plot_num);
    plot(1:2,resp1_miceAVG_d1d2n(:,i),"Color",colors(i),'LineWidth',3); hold on
end
xlim([0.5,2.5]);
ylim([-1,3]);
xticks([1, 2]);
xticklabels({'no stim', 'stim'});
title(['Task 4s Response Control']);

for i = mice_num/2+1:mice_num
    figure(RESP);
    subplot(1,2,plot_num+1);
    plot(1:2,resp1_miceAVG_d1d2n(:,i),"Color",colors(i),'LineWidth',3); hold on
end
xlim([0.5,2.5]);
ylim([-1,3]);
xticks([1, 2]);
xticklabels({'no stim', 'stim'});
title('ChrimsonR');

end


