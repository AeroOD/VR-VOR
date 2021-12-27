%Original VOR Code Adapted From:
%Rey-Martinez, J., Batuecas-Caletrio, A., Matiño, E., Trinidad-Ruiz, 
%G., Altuna, X., & Perez-Fernandez, N. (2018). Mathematical methods for 
%measuring the visually enhanced vestibulo?ocular reflex and preliminary
%results from healthy subjects and patient groups. 
%Frontiers in neurology, 9, 69. https://doi.org/10.3389/fneur.2018.00069
%Original source files: https://github.com/bendermh/VVOR
%
%***Requires Signal Processing Toolbox To Function ***

%% Instructions
%In this study, Subjects were given the last name of S01, S02, etc.
% in the Otometrics vHIT software. Condition was written G1, G2, or G3
%for 3 different gain demand conditions for the first name. 
%You can replace G with 'C' or another string to suit your needs. 
%
% If executing this script manually, enter Subject and Condition in the
% comand window before proceeding.
% Example:
%
% Suject = 2;
% Condition = 3;
%
% If using AutomateAllSubjects.m to automate the analysis for all subjects
% and conditions, then Subject & Condition will automatically be loaded. 


%% Load Files

%Find and load VOR File based on subject number and condition in the
%filename. For example, find S02_G3_2019_05_21_11_48_41.txt in sub folder.
vorfile=dir(strcat('**/S',num2str(Subject,'%02.f'),'_G',num2str(Condition),"*.txt")); %finds txt file
vor = vorimport(strcat(vorfile.folder,'/',vorfile.name), [1,inf]); %table with all raw VOR data.

%Import eye calibration Values for each measurement & VOR Start/End Times
%from the XML files produced by the Otometrics software. For example, find
% S02_G3_2019_05_21_11_48_41.xml in the subfolder.
calibfile=dir(strcat('**/S',num2str(Subject,'%02.f'),'_G',num2str(Condition),"*.xml")); %finds xml file
[calibration, VORStartTime, VOREndTime] = xmlimport(strcat(calibfile.folder,'/',calibfile.name), [5,inf]);

%Load file with calibration values. The Otometrics software may crash when
%the buffer is to great. In order to aleviate the need to re-calibrate and
%stop the experiment, one may calibrate at the begining of the session and
%record the calibration value for each subject and condition. This can be
%added to CalibrationValues.mat manually or you can replace this code
%with another import function of an Excel document, etc. The sample
%CalibrationValues.mat file contains pilot Subject 2 and 3 conditions.

load 'CalibrationValues.mat'
%% VOR Table

%The term "Block" refers to the experimental measurement block. Each time
%VOR collection starts and stops is a separate block. In the text file,
%these will be imported as lines of NaN.

%Add Subject, Condition, and Block Numbers to the VOR raw data table.
vor.Subject=ones(size(vor.Time))*Subject;
vor.Condition=ones(size(vor.Time))*Condition;
vor.Block=cumsum(isnan(vor.Time));

% Replace Missing Eye Data with NaN. This is to ensure that later analysis
%does not mistake missing eye data for 0 degrees per second.
vor.Eye(vor.Eye==0)=NaN;
vor.Eye_Vert(vor.Eye_Vert==0)=NaN;

%%Create New Tables
RawData=table();
CycleSummary=table();
SaccadeSummary=table();
BlockSummary=table();


%% Extract Blocks for working block by block

for block=1:1:max(vor.Block)
    BlockEx=vor(vor.Block==block & vor.Time>1,:); %Ensures rows of all NaN are not included.
    
    %Create New Tables
    newblock=table();
    cycledetail=table();
    
    %DELETE? sz=size(BlockEx,1); %size of block
    %DELETE? sz_ones=ones(size(BlockEx,1),1); %make column of ones the same length as the block
    
    
    BlockEx.TimeS=(BlockEx.Time-min(BlockEx.Time))/10000000; % Starts this extracted block time at 0 and converts to seconds.
    BlockEx.Position=cumtrapz(BlockEx.TimeS,BlockEx.Head); % Determines horizontal head postion based on velocity data.
    BlockEx.VerticalVelocity=(BlockEx.RALP+BlockEx.LARP)/2; % Determines average vertical head movement velocity using RALP & LARP velocity.
    BlockEx.VerticalPosition=cumtrapz(BlockEx.TimeS,BlockEx.VerticalVelocity); % Determines vertical head position based on average vertical head velocity.
    BlockEx.FiltEye=medfilt1(BlockEx.Eye,30,'omitnan','truncate'); %Filters horizontal eye movement velocity with a 1-dimesional median filter, window length of 30. Ignores missing eye data.
    BlockEx.FiltEye(isnan(BlockEx.Eye))=NaN; %Makes any originally missing eye data as missing in the filtered eye data, so that filtered data is not interpolated for analysis.
    BlockEx.EyeDiff=BlockEx.Eye-BlockEx.FiltEye; % Calculates diffrence beteween raw and filtered eye veolocity data.

    % Apply calibration
    BlockEx.SubjCalibration=CalibrationValues.Calibration(CalibrationValues.Subject==Subject & CalibrationValues.Condition==Condition).*ones(size(BlockEx,1),1); % Pulls subject & condition calibration value from CalibrationValues.mat
    BlockEx.BlockCalibration=calibration(block).*ones(size(BlockEx,1),1); % Pulls calibration value actually measured this block (comes from XML file via xmlimport.m). Otometrics default is 21.
    BlockEx.Correction_Factor=BlockEx.BlockCalibration./BlockEx.SubjCalibration; %Velocity correction factor based on subject's actual calibration requirement and what was recorded by Otometrics software.
    BlockEx.CalibEye=BlockEx.Eye.*BlockEx.Correction_Factor; 
    BlockEx.FiltCalibEye=medfilt1(BlockEx.CalibEye,30,'omitnan','truncate');
    BlockEx.FiltCalibEye(isnan(BlockEx.Eye))=NaN;
    BlockEx.CalibEyeDiff=BlockEx.CalibEye-BlockEx.FiltCalibEye;
    
    %Identify Saccades. Saccades are identified as the absolute value of the
    % (raw eye velocity) - (filtered eye velocity) >= 100 degrees per
    % second. Saccades are marked in the SaccadeMarker with a value of 1.
    [SaccadePeaks, SaccadeLocs, SaccadeWidths, SaccadeProm]=findpeaks(abs(BlockEx.CalibEyeDiff),'MinPeakHeight',100,'MinPeakWidth',2,'MaxPeakWidth',40,'MinPeakDistance',10);
    BlockEx.SaccadeMarker(SaccadeLocs(SaccadePeaks<600))=1;
    
    %Identify Head Velocity Peaks >35 degrees per second.
    [HeadPeaks, HeadPeaksLocs, HeadPeaksW, HeadPeaksP] = findpeaks(abs(BlockEx.Head),'MinPeakProminence',35,'MinPeakDistance',25);
    
    NumHeadPeaks=size(HeadPeaksLocs,1); %Counts number of peaks
    
    %Excludes first and last cycle
    PeakCount=0;
    for headpks=2:1:(NumHeadPeaks-1)
        PeakCount=PeakCount+1;
        BlockEx.PeakCycle(HeadPeaksLocs(headpks):HeadPeaksLocs(headpks+1))=PeakCount;
    end
    
    %Identifies half of a sinusoidal cycles (positive & negative)
    FirstCycle=BlockEx.Head(BlockEx.PeakCycle==1);
    if FirstCycle(1)>0
        for cycles=1:2:max(BlockEx.PeakCycle)
            BlockEx.HalfCycle(BlockEx.Head<0 & (BlockEx.PeakCycle==cycles | BlockEx.PeakCycle==cycles+1))=cycles;
            BlockEx.HalfCycle(BlockEx.Head>0 & (BlockEx.PeakCycle==cycles+1 | BlockEx.PeakCycle==cycles+2))=cycles+1;
        end
    else
        for cycles=1:2:max(BlockEx.PeakCycle)
            BlockEx.HalfCycle(BlockEx.Head>0 & (BlockEx.PeakCycle==cycles | BlockEx.PeakCycle==cycles+1))=cycles;
            BlockEx.HalfCycle(BlockEx.Head<0 & (BlockEx.PeakCycle==cycles+1 | BlockEx.PeakCycle==cycles+2))=cycles+1;
        end
    end
    BlockEx.HalfCycle(BlockEx.HalfCycle==max(BlockEx.HalfCycle))=0; %gets rid of tail

    
    nanIdx=isnan(BlockEx.Head) | isnan(BlockEx.FiltEye) | BlockEx.HalfCycle==0;
    CalibnanIdx=isnan(BlockEx.Head) | isnan(BlockEx.FiltCalibEye) | BlockEx.HalfCycle==0; 
    
    for HalfCycle=1:1:max(BlockEx.HalfCycle)
      if nnz(~isnan(BlockEx.FiltEye(BlockEx.HalfCycle==HalfCycle)))>10 && nnz(~isnan(BlockEx.FiltCalibEye(BlockEx.HalfCycle==HalfCycle))) >10
            
        if mean(BlockEx.Head(~nanIdx & BlockEx.HalfCycle==HalfCycle))>1 %identify head direction. Poistive = Right; Negative = Left.
            CycleVORDirection=1;
        else 
            CycleVORDirection=-1;
        end
        
        
        CyclePctEyeMissing=sum(isnan(BlockEx.CalibEyeDiff(BlockEx.HalfCycle==HalfCycle)))/length((BlockEx.CalibEyeDiff(BlockEx.HalfCycle==HalfCycle)));
        CycleCalibEyeStd=std(BlockEx.CalibEyeDiff(BlockEx.HalfCycle==HalfCycle),'omitnan');
        
        PeakHeadVelocity=max(abs(BlockEx.Head(~nanIdx & BlockEx.HalfCycle==HalfCycle)));
        HeadAmplitude=max(BlockEx.Position(~nanIdx & BlockEx.HalfCycle==HalfCycle))-min(BlockEx.Position(~nanIdx & BlockEx.HalfCycle==HalfCycle));
        HalfCycleTime=max(BlockEx.TimeS(~nanIdx & BlockEx.HalfCycle==HalfCycle))-min(BlockEx.TimeS(~nanIdx & BlockEx.HalfCycle==HalfCycle));
        
        CycleVOR=BlockEx.Head(~nanIdx & BlockEx.HalfCycle==HalfCycle)\BlockEx.FiltEye(~nanIdx & BlockEx.HalfCycle==HalfCycle);

        
        %Determine if Valid Cycle based on expected VOR to ensure spurious
        %data are not included. In this case, we selected a minimum of 0.3
        %and maximum of 2.5, for both raw and calibrated VOR.
        if CycleVOR>0.3 && CycleVOR<2.5 
            CycleSelected=1;
        else
            CycleSelected=0;
        end
        
        CalibCycleVOR=BlockEx.Head(~CalibnanIdx & BlockEx.HalfCycle==HalfCycle)\BlockEx.FiltCalibEye(~CalibnanIdx & BlockEx.HalfCycle==HalfCycle);
        if CalibCycleVOR>0.3 && CalibCycleVOR<2.5 %Determine if Valid Cycle based on expected VOR
            CalibCycleSelected=1;
        else
            CalibCycleSelected=0;
        end
        
        
        %Omits any half cycles >200 lines (800 msec). This would typically be
        %too slow for taxing the VVOR.
        if length(BlockEx.Head(~nanIdx & BlockEx.HalfCycle==HalfCycle))>200
            CycleSelected=0;
            CalibCycleSelected=0;
        end
        
        cyclesz=size(BlockEx(BlockEx.HalfCycle==HalfCycle,1),1);
        cycleresults=table(CycleVORDirection, PeakHeadVelocity, HeadAmplitude, HalfCycleTime, CycleVOR, CycleSelected, CalibCycleVOR, CalibCycleSelected, CyclePctEyeMissing, CycleCalibEyeStd);
        
        tmpcycledetail=[table(Subject, Condition, block) table(HalfCycle) cycleresults];
        tmpnewblock=[BlockEx(BlockEx.HalfCycle==HalfCycle,1:end) repmat(cycleresults(:,:),size(BlockEx(BlockEx.HalfCycle==HalfCycle,1),1),1)]; 

        newblock=[newblock;tmpnewblock];
        cycledetail=[cycledetail;tmpcycledetail];
      else
          %do nothing
      end
    end

    %Calculate Positive, Negative, & Total VOR for both uncalibrated and
    %calibrated eye
    PosVOR=newblock.Head(~isnan(newblock.FiltEye) & newblock.CycleSelected==1 & newblock.Head>0)\newblock.FiltEye(~isnan(newblock.FiltEye) & newblock.CycleSelected==1 & newblock.Head>0);
    NegVOR=newblock.Head(~isnan(newblock.FiltEye) & newblock.CycleSelected==1 & newblock.Head<0)\newblock.FiltEye(~isnan(newblock.FiltEye) & newblock.CycleSelected==1 & newblock.Head<0);
    TotalVOR=newblock.Head(~isnan(newblock.FiltEye) & newblock.CycleSelected==1)\newblock.FiltEye(~isnan(newblock.FiltEye) & newblock.CycleSelected==1);
    
    CalibPosVOR=newblock.Head(~isnan(newblock.FiltCalibEye) & newblock.CalibCycleSelected==1 & newblock.Head>0)\newblock.FiltCalibEye(~isnan(newblock.FiltCalibEye) & newblock.CalibCycleSelected==1 & newblock.Head>0);
    CalibNegVOR=newblock.Head(~isnan(newblock.FiltCalibEye) & newblock.CalibCycleSelected==1 & newblock.Head<0)\newblock.FiltCalibEye(~isnan(newblock.FiltCalibEye) & newblock.CalibCycleSelected==1 & newblock.Head<0);
    CalibTotalVOR=newblock.Head(~isnan(newblock.FiltCalibEye) & newblock.CalibCycleSelected==1)\newblock.FiltCalibEye(~isnan(newblock.FiltCalibEye) & newblock.CalibCycleSelected==1);
    
    BlockPctEyeMissing=sum(isnan(newblock.FiltCalibEye(newblock.CalibCycleSelected==1)))/length(newblock.FiltCalibEye(newblock.CalibCycleSelected==1));
    BlockCalibEyeStd=std(newblock.CalibEyeDiff(newblock.HalfCycle==CalibCycleSelected==1),'omitnan');
    
    %Do Fourier Analysis - requires filling gaps.
    
    newblock.fillgapsHead=fillgaps(newblock.Head);
    newblock.fillgapsCalibEye=fillgaps(newblock.FiltCalibEye);
    
    [fHead,P1Head] = fourier(newblock.fillgapsHead);
    [fEye,P1Eye] = fourier(newblock.fillgapsCalibEye);
    
    
    [MaxP1Head, MaxP1HeadLoc]=max(P1Head);
    [MaxP1Eye, MaxP1EyeLoc]=max(P1Eye);
    FFT_Max_Head_Freq=fHead(MaxP1HeadLoc);
    FFT_Max_Eye_Freq=fEye(MaxP1EyeLoc);
    FFT_Gain_at_MaxHeadFreq=MaxP1Head\P1Eye(MaxP1HeadLoc);
    FFT_Gain_Range=P1Head(fHead>0.5 & fHead<1.5)\P1Eye(fHead>0.5 & fHead<1.5);

    %add to tables
    
    %Determine Number of Saccades & Characteristics
    NumValidCycles=length(unique(newblock.HalfCycle(newblock.CalibCycleSelected==1)))./2;
    NumSaccades=length(newblock.SaccadeMarker(newblock.SaccadeMarker==1 & CalibCycleSelected==1));
    SaccadesPerCycle=NumSaccades./NumValidCycles;
    
    tmpSaccadeSummary=table(newblock.TimeS(newblock.SaccadeMarker==1 & CalibCycleSelected==1),newblock.Head(newblock.SaccadeMarker==1 & CalibCycleSelected==1),newblock.CalibEye(newblock.SaccadeMarker==1 & CalibCycleSelected==1));
    tmpSaccadeSummary.Properties.VariableNames = {'TimeS','Head','CalibEye'};
    tmpSaccadeSummary.SacDirection(tmpSaccadeSummary.CalibEye./tmpSaccadeSummary.Head>0)=1;
    tmpSaccadeSummary.SacDirection(tmpSaccadeSummary.CalibEye./tmpSaccadeSummary.Head<0)=-1;
    AvgSacVel=mean(abs(tmpSaccadeSummary.CalibEye)); % Average saccade velocity in degrees per second.
    AvgSacDir=mean(tmpSaccadeSummary.SacDirection); 
    tmpSaccadeSummary.Subject(:,1)=Subject;
    tmpSaccadeSummary.Condition(:,1)=Condition;
    tmpSaccadeSummary.Block(:,1)=block;
    tmpSaccadeSummary = movevars(tmpSaccadeSummary, 'Subject', 'Before', 'TimeS');
    tmpSaccadeSummary = movevars(tmpSaccadeSummary, 'Condition', 'After', 'Subject');
    tmpSaccadeSummary = movevars(tmpSaccadeSummary, 'Block', 'After', 'Condition');
    
    
    VORresults=table(PosVOR, NegVOR, TotalVOR, CalibPosVOR, CalibNegVOR, CalibTotalVOR, NumValidCycles, NumSaccades, SaccadesPerCycle, AvgSacVel, AvgSacDir, FFT_Max_Head_Freq, MaxP1Head, FFT_Max_Eye_Freq, MaxP1Eye, FFT_Gain_at_MaxHeadFreq, FFT_Gain_Range, BlockPctEyeMissing, BlockCalibEyeStd);
    
    tmpSaccadeSummary=[tmpSaccadeSummary repmat(VORresults(:,:),size(tmpSaccadeSummary,1),1)];

    tmpBlockSummary=[table(Subject, Condition, block) VORresults];
    BlockSummary=[BlockSummary; tmpBlockSummary];
    
    tmpCycleSummary=[cycledetail repmat(VORresults(:,:),size(cycledetail,1),1)];
    CycleSummary=[CycleSummary; tmpCycleSummary];
    SaccadeSummary=[SaccadeSummary; tmpSaccadeSummary];
    RawData=[RawData; newblock];
    
end

% Save .mat files. First, check to see if a folder called Analyzed exists.
if ~exist(strcat(pwd,"/Analyzed/"))
    mkedir(strcat(pwd,"/Analyzed/"));
end

save(strcat(pwd,"/Analyzed/",num2str(Subject),"_",num2str(Condition),".mat"))