%% Subject-wise ROI activation extraction for Agency task
%
% This script extracts the mean activity within the rTPJ ROI (NF training
% mask), for each subject and for each of the two sessions of the agency task.
% It then saves this in an Excel file.
% The same code is used for extraction of ROI activity during the NF training task.
% (Adjust to extract from each run/session as required.)


clc; clear all; close all; 

% Define folders
 % Add path to utilities
addpath(genpath('./fMRI_ROI'));savepath;

% Main path to data folders
mainDir = './data/';

% Global parameters
run_settings = 1; % run directly the matlabbatch

% Find subject folders
S = dir(fullfile(mainDir, 'P*')); % get all P-Code-Folders / subjects
dirFlags = [S.isdir]; %isdir returns a 1 if S is a directory and 0 otherwise. 
subFolders = S(dirFlags);

% Create cell with the subject folder names
n_subjects= size(S,1);
subject_name= {};
for k = 1 : n_subjects
    subject_name{k} = subFolders(k).name;
end

Session = {'V01', 'V03'};


%% Get Mean activation per session within rTPJ
sess_meanValues = cell(n_subjects,2);

p_values_ttest = zeros(n_subjects,1);
DataWithinROI = zeros(672,2);

for subj = 1:n_subjects % for each subject
    for sess = 1:length(Session) % for each session / run 
        if sess == 1
            conFile = './contrastFile_sess1'; %Turbulence > Baseline contrast session 1
        elseif sess == 2
            conFile = './contrastFile_sess2'; %Turbulence > Baseline contrast session 2

        conMapVol = spm_vol(conFile);
        conMapData = spm_read_vols(conMapVol);

        roiMaskFile = './rTPJ_mask.nii';

        % Read the voxel data
        roiMaskVol = spm_vol(roiMaskFile);
        roiMaskData = spm_read_vols(roiMaskVol);
        
        
        % Extract data within the ROI (remove the background)
        DataWithinROI(:,sess) = conMapData(roiMaskData > 0);
        
       
        % Display mean contrast values within the ROI
        sess_meanValues{subj,sess} = mean(DataWithinROI(:,sess));
        fprintf('Mean BOLD signal within the ROI for Session %d (Turb > Base): %f\n',sess, sess_meanValues{subj,sess});
        end
   end

end

% Save in excel file
row_names = {'P-codes'};
diff_sess2_sess1 = cell2mat(sess_meanValues(:,2)) - cell2mat(sess_meanValues(:,1));
concat_matrix = [cell2mat(sess_meanValues), diff_sess2_sess1];
sess_table_concat = array2table(concat_matrix, 'VariableNames', column_names, 'RowNames', row_names);
writetable(sess_table_concat,'Output_directory','Sheet','Mean_Values','WriteRowNames', true);

















