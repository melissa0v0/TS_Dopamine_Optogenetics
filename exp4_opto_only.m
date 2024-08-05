%% Response for all mice [trial type x num of trials x mice]

close all
mice_num = 3;

%% acquire data

Trace = zeros(10001,mice_num);
Resp = zeros(1,mice_num);
for i = 1:mice_num
    mouse_name = [group, num2str(i)];
    [trace, resp, ste_resp] = SingleAnimal_opto_only(mouse_name);
    Trace(:,i) = trace;
    Resp(i) = resp;
end

%% plot traces

TRACE = figure();
set(gcf, 'Position', [50,50, 800, 300]);
set(gcf,'color','w');
plotWin = [-2000:8000];
xrange = [-500,4500];

legend_name = {'opto only'};

rectangle('Position', [0, -1, 4000, 6], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on
plot([-3000, 8000], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);hold on

m_plot = mean(Trace,2);  
s_plot = std(Trace')/sqrt(3);
errorbar_patch(plotWin,m_plot,s_plot,'r',false);hold on

title("opto only")
h=gca;
xlabel('time(ms)');
ylabel('z-score');
xlim(xrange);
ylim([-1,3])

%% test sig
[h,p_chr,ci,stat_chr] = ttest(mean(Trace(2000:6000,:)),0)


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


