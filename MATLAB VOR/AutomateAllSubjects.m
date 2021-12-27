%This will automate the AnalyzeVOR function through all subjects.
%This requires all subjects & condition files to be in the subfolders and
%CalibrationValues.mat to be properly filled out.

SubejctStart=2; %<-- Put the first subject you want to evaluate. The example only has Subject #2.
SubjectEnd=2; %<-- Put the last subject you want to evaluate. The example only has Subject #2.
ConditionStart=1; %<-- Put the first condition you want to evaluate. The example starts with condition 1.
ConditionEnd=3; %<-- Put the last condition you want to evaluate. The example has 3 conditions.

all_Blocks=table();
all_Cycles=table();
all_Saccades=table();

if ~exist(strcat(pwd,"/SummaryData/"))
            mkdir(strcat(pwd,"/SummaryData/"));
end

summaryfile = datestr(datetime(),'yyyy_mm_dd_HH_MM_SS');

for Subject=SubejctStart:1:SubjectEnd
    for Condition=ConditionStart:1:ConditionEnd
        
        disp(strcat("Starting:",num2str(Subject),"_",num2str(Condition)))
        AnalyzeVOR
        
        all_Blocks=[all_Blocks; BlockSummary];
        all_Cycles=[all_Cycles; CycleSummary];
        all_Saccades=[all_Saccades; SaccadeSummary];
        
       

        save(strcat(pwd,"/SummaryData/Summary_",summaryfile,".mat")) %This will save the file with each iteration, in case it crashes.

        
       
    end
end

writetable(all_Blocks,strcat(pwd,"/SummaryData/Summary_Block_",summaryfile,".xlsx"))
writetable(all_Cycles,strcat(pwd,"/SummaryData/Summary_Cycles_",summaryfile,".xlsx"))
writetable(all_Saccades,strcat(pwd,"/SummaryData/Summary_Saccades_",summaryfile,".xlsx"))