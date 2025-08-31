%% Extraction of Results from Behavioral Data of the Neurofeedback Training

% This script automatically summarizes the behavioural NF training data regarding indices for the Neurofeedback project. 
% This code extracts the results from the single excel-file from each
% subject produces a summary and saves it in the Analysis folder (specified below).

% The folders have to be organized as follows: 
% Participant folder:
%  - P-XXX
%       - V01
%           - NF_TPJ_Run1
%               - behav
%                   - Pgame_P-XXX_1.mat
%                   - Pgame_P-XXX_2.mat
%                   - Pgame_P-XXX_3.mat
%                   - Pgame_P-XXX_4.mat
%                   - Pgame_P-XXX_5.mat
%                   - Pgame_P-XXX_6.mat

% The script is written for a total of 2 conditions: 
%   1 = Baseline
%   2 = Turbulence 40%



% v1.0 Eliane Mueller (July 2024)

clear all; clc;

%% Define your path and import all folders (P-codes)

RootPath = '/Users/elianemueller/Documents/Neurofeedback/Online_NF/MRI_data/data/';
folder = dir(fullfile(RootPath,'P*')); %define all P-Code-Folders as variable
subject_name = {}; %make an empty variable/list

for i = 1:length(folder)
    subject_name{i} = folder(i).name;
end
clear i;

Session = {'V01','V01', 'V01', 'V02', 'V02', 'V02','V03', 'V03', 'V03'};
run_var = {'5_FND_NF_TPJ_Run1','5_FND_NF_TPJ_Run2', '5_FND_NF_TPJ_Run3', '3_FND_NF_TPJ_Run1','3_FND_NF_TPJ_Run2','3_FND_NF_TPJ_Run3', '3_FND_NF_TPJ_Run1','3_FND_NF_TPJ_Run2','3_FND_NF_TPJ_Run3'};


%% Accuracy and Errors
% Load Pgame.mat and extract HitX and Hit O - per subject

HitX = zeros(6,9);
HitO = zeros(6,9);
HitX_O = zeros(6,9);

HitX_allSubj = zeros(length(subject_name),9);
HitO_allSubj = zeros(length(subject_name),9);
HitX_O_allSubj = zeros(length(subject_name),9);

subj_cell = cell(length(subject_name),1);


% Define data directory
for subj = 1: length(subject_name)
    for r = 1:9 
        for i = 1:6
            data_dir = fullfile([RootPath subject_name{subj} filesep Session{r} filesep run_var{r} filesep 'behav' filesep]);
            pgame_dir = dir(fullfile([data_dir, 'Pgame*'])); % get the directory of the P file and load the P.game 
            load(fullfile(pgame_dir(i).folder, pgame_dir(i).name))
            numX = size(P.game.hitX,1); %number of Hits is saved in P.game.hitX
            HitX(i,r) = (numX/48)*100; %define percentage of X that were caught (out of a total of 48)
            numO = size(P.game.hitDistractors,1); %number of Errors is saved in P.game.hitDistractors
            HitO(i,r) = (numO/48)*100;
            HitX_O(i,r) = HitX(i,r)- HitO(i,r);
        end
        HitX_allSubj(subj,r)= mean(HitX(:,r));
        HitO_allSubj(subj,r)= mean(HitO(:,r));
        HitX_O_allSubj(subj,r) = mean(HitX_O(:,r));
    end

    %save HitX and HitO as an excel file (two sheets) - per subject
    filename = fullfile(RootPath, subject_name{subj}, Session{1}, run_var{1}, 'behav', [subject_name{subj}, '_HitXO.xlsx']);

    header = {'Run1', 'Run2','Run3', 'Run4','Run5','Run6','Run7','Run8','Run9'}; %add headers to all excel-sheets
        
    writecell([header ; num2cell(HitX)], filename, 'Sheet', 'HitX')
    writecell([header ; num2cell(HitO)], filename, 'Sheet', 'HitO')
    writecell([header ; num2cell(HitX_O)], filename, 'Sheet', 'HitX - HitO')

    subj_cell(subj)=subject_name(subj);

    
end


% save for all subjects 
filename2 = fullfile(RootPath,'group_analysis','behav_results','NF', 'HitX_HitO_all_subjects.xlsx');
header = {'Subject','Run1', 'Run2','Run3', 'Run4','Run5','Run6','Run7','Run8','Run9'}; %add headers to all excel-sheets
        
writecell([header ; horzcat(subj_cell,num2cell(HitX_allSubj))], filename2, 'Sheet', 'HitX_allSubj')
writecell([header ; horzcat(subj_cell,num2cell(HitO_allSubj))], filename2, 'Sheet', 'HitO_allSubj')
writecell([header ; horzcat(subj_cell,num2cell(HitX_O_allSubj))],filename2, 'Sheet', 'HitX-HitO_allSubj')


%% Time spent in Turbulence

turb_perc = zeros(6,9);
turb_perc_all_subj = zeros(length(subject_name),9);

subj_cell = cell(length(subject_name),1);

for subj = 1: length(subject_name)
    for r = 1:9 
        for i = 1:6
            data_dir = fullfile([RootPath subject_name{subj} filesep Session{r} filesep run_var{r} filesep 'behav' filesep]);
            pgame_dir = dir(fullfile([data_dir, 'Pgame*'])); % get the directory of the P file and load the P.game 
            load(fullfile(pgame_dir(i).folder, pgame_dir(i).name))
            turb_perc(i,r)=length(find(P.game.all_diff_lvl_frames==2))/...
                length(P.game.all_diff_lvl_frames)*100;
        end
        turb_perc_all_subj(subj,r) = mean(turb_perc(:,r));
    end

     %save Turb as an excel file (two sheets) - per subject
    filename = fullfile(RootPath, subject_name{subj}, Session{1}, run_var{1}, 'behav', [subject_name{subj}, '_TimeInTurb.xlsx']);

    header = {'Run1', 'Run2','Run3', 'Run4','Run5','Run6','Run7','Run8','Run9'}; %add headers to all excel-sheets
        
    writecell([header ; num2cell(turb_perc)], filename, 'Sheet', 'Turbulence')

    subj_cell(subj)=subject_name(subj);


end

% save for all subjects (summarising the runs) 
filename2 = fullfile(RootPath,'group_analysis', 'behav_results','NF', 'TimeInTurb_all_subjects.xlsx');
header = {'Subject','Run1', 'Run2','Run3', 'Run4','Run5','Run6','Run7','Run8','Run9'}; %add headers to all excel-sheets
        
writecell([header ; horzcat(subj_cell,num2cell(turb_perc_all_subj))], filename2, 'Sheet', 'Turbulence')





