function plotVOR(X1, YMatrix1, Subject, Condition, Block)

BlockNumber=mean(Block);


%CREATEFIGURE(X1, YMatrix1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 10-Jul-2019 16:11:42

% Create figure
VOR = figure('Position',[0 0 1920 1080],'visible','off');

% Create axes
axes1 = axes('Parent',VOR);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',3);
set(plot1(1),'DisplayName','Eye (Unfiltered)',...
    'Color',[0.313725501298904 0.313725501298904 0.313725501298904],...
    'LineWidth',0.5);
set(plot1(2),'DisplayName','Head','Color',[0 0 1]);
set(plot1(3),'DisplayName','Eye (Filtered)','Color',[1 0 0]);


% Create ylabel
ylabel('Velocity (�/sec)');

% Create xlabel
xlabel('Time (Seconds)');

% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-200 200]);
xlim([0 max(X1)]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',24);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Orientation','horizontal','Location','south','FontSize',24);
title(axes1,strcat('Subject: ',num2str(Subject),'  Condition: ',num2str(Condition),'  Block:',num2str(BlockNumber)))

Fig_Folder=strcat(pwd,'/Images/Subj_',num2str(Subject),'/Con_',num2str(Condition),'/VOR_fig');
PNG_Folder=strcat(pwd,'/Images/Subj_',num2str(Subject),'/Con_',num2str(Condition),'/VOR_png');
Filename=strcat('Subj_',num2str(Subject),'_Con_',num2str(Condition),'_Blk_',num2str(BlockNumber),'_VOR');

if ~exist(Fig_Folder, 'dir')
    mkdir(Fig_Folder)
end

if ~exist(PNG_Folder, 'dir')
    mkdir(PNG_Folder)
end

saveas(VOR, fullfile(Fig_Folder,strcat(Filename,'.fig')))
print(VOR,'-dpng','-r500',fullfile(PNG_Folder,strcat(Filename,'.png')))
close(VOR)