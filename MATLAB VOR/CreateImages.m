%Images can take a lot of processing time & power, therefore, they are not
%automatically created by default.

%Data must be analyzed with AutomateAllSubjects or AnalyzeVOR first.

SubejctStart=2; %<-- Put the first subject you want to create images. The example only has Subject #2.
SubjectEnd=2; %<-- Put the last subject you want to create images. The example only has Subject #2.
ConditionStart=1; %<-- Put the first condition you want to create images. The example starts with condition 1.
ConditionEnd=3; %<-- Put the last condition you want to create images. The example has 3 conditions.
BlockStart=34; %<-- Put the first block you want to create images.
BlockEnd=34; %<-- Put the last block you want to create images. If number exceeds available blocks, it will proceed until finished with the max number of blocks.

for Subject=SubejctStart:1:SubjectEnd
    for Condition=ConditionStart:1:ConditionEnd
        
        load(strcat(pwd,"/Analyzed/",num2str(Subject),"_",num2str(Condition),".mat")); %Loads analyzed VOR Data.
        
        if BlockEnd > max(RawData.Block)
            BlockEnd = max(RawData.Block);
        end
        
        for Block=BlockStart:1:BlockEnd
       
            plotVOR(RawData.TimeS(RawData.Block==Block), [RawData.CalibEye(RawData.Block==Block) RawData.Head(RawData.Block==Block) RawData.FiltCalibEye(RawData.Block==Block)], Subject, Condition, Block)
            plotHvE(RawData.Head(RawData.Block==Block), RawData.FiltCalibEye(RawData.Block==Block), Subject, Condition, Block)

        end
    end
end
