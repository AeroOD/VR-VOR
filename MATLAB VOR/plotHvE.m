function plotHVE(ValidHead, ValidFilteredEye, Subject, Condition, Block)



idx = ValidFilteredEye./ValidHead>0;
ValidHead = ValidHead(idx);
ValidFilteredEye = ValidFilteredEye(idx);

BlockNumber=mean(Block);

% Create Head vs Eye Scatter Plot
HvE = figure('Position',[0 0 1920 1080],'visible','off');

set(gca,'FontSize',20)
hold on
plot(ValidHead,ValidFilteredEye,'Color','[0 0.4470 0.7410]','LineStyle','none','marker','.')
xlabel('Head Velocity (°/sec)','FontSize',24);
ylabel('Eye Velocity (°/sec)','FontSize',24);
h=polyfit(ValidHead, ValidFilteredEye, 1);
text(-80,40,strcat('Slope: ',num2str(h(1))),'FontSize',24)
xf = [min(ValidHead), max(ValidHead)];
yf = polyval(h,xf);
plot(xf, yf, 'Color','[0.8500 0.3250 0.0980]','LineStyle','-','marker','none', 'LineWidth', 2)
title(strcat('Subject: ',num2str(Subject),'  Condition: ',num2str(Condition),'  Block:',num2str(BlockNumber)))
hold off


%Fig_Folder=strcat(pwd,'/',num2str(Subject),'/',num2str(Condition),'/HvE_fig');
%PNG_Folder=strcat(pwd,'/',num2str(Subject),'/',num2str(Condition),'/HvE_png');
%Filename=strcat(num2str(BlockNumber),'_HvE');

Fig_Folder=strcat(pwd,'/Images/Subj_',num2str(Subject),'/Con_',num2str(Condition),'/HvE_fig');
PNG_Folder=strcat(pwd,'/Images/Subj_',num2str(Subject),'/Con_',num2str(Condition),'/HvE_png');
Filename=strcat('Subj_',num2str(Subject),'_Con_',num2str(Condition),'_Blk_',num2str(BlockNumber),'_HvE');


if ~exist(Fig_Folder, 'dir')
    mkdir(Fig_Folder)
end

if ~exist(PNG_Folder, 'dir')
    mkdir(PNG_Folder)
end

set(HvE, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
saveas(HvE, fullfile(Fig_Folder,strcat(Filename,'.fig')))
print(HvE,'-dpng','-r500',fullfile(PNG_Folder,strcat(Filename,'.png')))
close(HvE)