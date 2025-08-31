%% PREPROCESSING OF BEHAVIORAL AGENCY DATA FOR Neurofeedback project

% This script automatically preprocesses the behavioral agency data regarding timing for the Agency project. 
% The data-output from scanner is directly converted into a very nice
% excelsheet (called PXX_resultsTime.xlsx). This excelsheet is needed for the first level analysis. 

%v1.0 Giuseppe Zito 
%v2.0 Janine Buehler (Feb 2021)
%v3.0 Serafeim Loukas (Oct 2021)
%v4.0 Serafeim Loukas (Dec 2022)
%v5.0 Eliane Mueller (Aug 2023) - only 3 diff levels left (PV, Baseline &
%Turbulence) 

clear all; clc;

%% Define your path and import all folders (P-codes)

RootPath = '/Users/elianemueller/Documents/Neurofeedback/Online_NF/MRI_data/data/';
folder = dir(fullfile(RootPath,'P*')); %define all P-Code-Folders as variable
mySubjects = {}; %make an empty variable/list
Sequence = {'funct', 'struct'};

for i = 1:length(folder)
    mySubjects{i} = folder(i).name;
end
clear i;


%% Start the Processing of the Data


% we start a loop over all folders (j, P-Codes)
for j = 12
    disp(['***Start processing Behavioral Data of: ' (mySubjects{j})]); 
    visit = 1;
    for v = 1:2
        if visit == 1 
            %Go into the according folder - and extract the textfile called "Results"
            AgencyFolder = char(fullfile(RootPath, mySubjects(j), 'V01','2_Agency_Game','behav'));
        else 
            AgencyFolder = char(fullfile(RootPath, mySubjects(j), 'V03', '4_Agency_Game', 'behav'));
        end
    
        mytmp =  dir(fullfile(AgencyFolder, '*_Results_*.txt'));
        fileNames{j} = fullfile(AgencyFolder, mytmp.name);
         
        %Import Results Textfile
        if ~isfile(fileNames{j}); continue; end %only continue, if there is no file yet
    
        data{j} = importdata(fileNames{j}); % import the excelfile to matlab
        gamePhase{j} = data{j}.data(:,1); % exacting gamephases only (0, 11, 21, 12, 22)
    
        %Go into the according folder - and extract the textfile called "Timing"
        mytmp2 =  dir(fullfile(AgencyFolder, '*_Timing_*.txt'));
        fileNames{j} = fullfile(AgencyFolder, mytmp2.name);
            
        %Import Timing Textfile
        data{j} = importdata(fileNames{j}); % import the excelfile to matlab
    
        % Transform hours/minutes in seconds
        timeInSec{j} = data{j}.data;
    
        %convert all colums containing hours to seconds (*60*60)
        timeInSec{j}(:,1) = data{j}.data(:,1)*60*60;
        timeInSec{j}(:,4) = data{j}.data(:,4)*60*60;
        timeInSec{j}(:,7) = data{j}.data(:,7)*60*60;
        timeInSec{j}(:,10) = data{j}.data(:,10)*60*60;
        timeInSec{j}(:,13) = data{j}.data(:,13)*60*60;
        timeInSec{j}(:,16) = data{j}.data(:,16)*60*60;
    
        %convert all colums containing minutes to seconds (*60)
        timeInSec{j}(:,2) = data{j}.data(:,2)*60;
        timeInSec{j}(:,5) = data{j}.data(:,5)*60;
        timeInSec{j}(:,8) = data{j}.data(:,8)*60;
        timeInSec{j}(:,11) = data{j}.data(:,11)*60;
        timeInSec{j}(:,14) = data{j}.data(:,14)*60;
        timeInSec{j}(:,17) = data{j}.data(:,17)*60;
    
        % Sum hours+minutes+seconds
        finalStartEndTime{j} = zeros(length(data{j}.data),6); %make a cell array according to length of data and fill it with zeros
        for i = 1:length(data{j}.data)%HERE 
            finalStartEndTime{j}(i,1) = sum(timeInSec{j}(i,1:3)); %Start of Game
            finalStartEndTime{j}(i,2) = sum(timeInSec{j}(i,4:6)); %End of Game
            finalStartEndTime{j}(i,3) = sum(timeInSec{j}(i,7:9)); %Start of JOP
            finalStartEndTime{j}(i,4) = sum(timeInSec{j}(i,10:12)); %End of JOP
            finalStartEndTime{j}(i,5) = sum(timeInSec{j}(i,13:15)); %Start of JOA
            finalStartEndTime{j}(i,6) = sum(timeInSec{j}(i,16:18)); %End of JOA
        end
    
        % OLD WITH PROBLEM
        %finalStartEndTime{j} = finalStartEndTime{j}-finalStartEndTime{j}(1,1);      % Subtract the first value because that's the time "0" % Ok
        %finalStartTime{j} = horzcat(gamePhase{j},finalStartEndTime{j}(:,[1 3 5])+1.5); % Gamephase, Timeonset Game, Timeonset JOP, Timeonset JOA % Fix cross is not included in the timing
        %finalDurationTime{j} = horzcat(gamePhase{j},finalStartEndTime{j}(:,[2 4 6]) - finalStartEndTime{j}(:,[1 3 5])); % Correct
        %finalDurationTime{j}(:,2) = finalDurationTime{j}(:,2)-1.5;  % Remove the time of the fixation cross, 1.5 seconds for the Game Phase 
        
        % NEW CORRECTED
        finalStartEndTime{j} = finalStartEndTime{j}-finalStartEndTime{j}(1,1);      % Subtract the first value because that's the time "0" % Ok
        finalStartTime{j} = horzcat(gamePhase{j},finalStartEndTime{j}(:,[1 3 5])+1.5); % Add 1.5 which is a close guess of the value of tStartGame1 variable because a fixation cross was displayed between the MR trigger and the first gamephase onset. Note: 1.5 is not the true value though!
        finalDurationTime{j} = horzcat(gamePhase{j},finalStartEndTime{j}(:,[2 4 6]) - finalStartEndTime{j}(:,[1 3 5])); % Correct, duration is 15 secs
        
        
        % Group according to game phases.
        pureVisualStart{j} = finalStartTime{j}(find(gamePhase{j}==1),:);
        pureVisualDuration{j} = finalDurationTime{j}(find(gamePhase{j}==1),:);
    
        baselineStart{j} = finalStartTime{j}(find(gamePhase{j}==2),:);
        baselineDuration{j} = finalDurationTime{j}(find(gamePhase{j}==2),:);
    
%         baselineHardStart{j} = finalStartTime{j}(find(gamePhase{j}==3),:);
%         baselineHardDuration{j} = finalDurationTime{j}(find(gamePhase{j}==3),:);
%     
        turbulenceStart{j} = finalStartTime{j}(find(gamePhase{j}==4),:);
        turbulenceDuration{j} = finalDurationTime{j}(find(gamePhase{j}==4),:);
    
%         turbulenceHardStart{j} = finalStartTime{j}(find(gamePhase{j}==5),:);
%         turbulenceHardDuration{j} = finalDurationTime{j}(find(gamePhase{j}==5),:);
    
       
        % Save in an excel file
        visitNo = num2str(visit);
        filename = fullfile(AgencyFolder,[mySubjects{j}, '_resultsTime', '_V0' visitNo '.xlsx']);  
    
        header = {'Condition', 'Gamephase','JoP', 'JoA'}; %add headers to all excel-sheets
        %xlswrite(filename,[header ; num2cell(pureVisualStart{1,j})],'PureVisualStart')
        %xlswrite(filename,[header ; num2cell(pureVisualDuration{1,j})],'PureVisualDuration') 
        
        % FIX: for macbook
        writecell([header ; num2cell(zeros(4,4))], filename, 'Sheet', 'Tabelle1')% empty to match previous structure by Janine
        writecell([header ; num2cell(pureVisualStart{1,j})], filename, 'Sheet', 'PureVisualStart')
        writecell([header ; num2cell(pureVisualDuration{1,j})], filename, 'Sheet', 'PureVisualDuration')
    
        %xlswrite(filename,[header ; num2cell(baselineEasyStart{1,j})],'BaselineEasyStart')
        %xlswrite(filename,[header ; num2cell(baselineEasyDuration{1,j})],'BaselineEasyDuration')
        
        % FIX: for macbook
        writecell([header ; num2cell(baselineStart{1,j})], filename, 'Sheet', 'BaselineStart')
        writecell([header ; num2cell(baselineDuration{1,j})], filename, 'Sheet', 'BaselineDuration')       
        
        %xlswrite(filename,[header ; num2cell(baselineHardStart{1,j})],'BaselineHardStart')
        %xlswrite(filename,[header ; num2cell(baselineHardDuration{1,j})],'BaselineHardDuration')
        
        % FIX: for macbook
        % writecell([header ; num2cell(baselineHardStart{1,j})], filename, 'Sheet', 'BaselineHardStart')
        % writecell([header ; num2cell(baselineHardDuration{1,j})], filename, 'Sheet', 'BaselineHardDuration')  
    
        %xlswrite(filename,[header ; num2cell(turbulenceEasyStart{1,j})],'TurbulenceEasyStart')
        %xlswrite(filename,[header ; num2cell(turbulenceEasyDuration{1,j})],'TurbulenceEasyDuration')
        
        % FIX: for macbook
        writecell([header ; num2cell(turbulenceStart{1,j})], filename, 'Sheet', 'TurbulenceStart')
        writecell([header ; num2cell(turbulenceDuration{1,j})], filename, 'Sheet', 'TurbulenceDuration') 
        
        %xlswrite(filename,[header ; num2cell(turbulenceHardStart{1,j})],'TurbulenceHardStart')
        %xlswrite(filename,[header ; num2cell(turbulenceHardDuration{1,j})],'TurbulenceHardDuration')
        
        % FIX: for macbook
        % writecell([header ; num2cell(turbulenceHardStart{1,j})], filename, 'Sheet', 'TurbulenceHardStart')
        % writecell([header ; num2cell(turbulenceHardDuration{1,j})], filename, 'Sheet', 'TurbulenceHardDuration') 

        visit = 3; 
    end
    
    disp(['***Finished processing Behavioral Data of ' num2str(j) ' out of ' num2str(length(folder)) ' Subjects.***']);
    disp(['Saved at ' filename])

end


