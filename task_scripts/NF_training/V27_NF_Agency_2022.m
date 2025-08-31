%% NF Training task

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Code written by Serafeim Loukas (serafeim.loukas@eplf.ch) on June 16, 2022.
%
% New version of the experiment for the real time NF project and task.
% The code receives in real-time the NF value from the brain and changes the difficulty
% level of the Agency task accordingly. Only dynamically adapted 
% game phase, no judgement phase.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  VERSION LOG, by Serafeim Loukas, serafeim.loukas@epfl.ch
%  v1 June 16, 2022, Modified: initial version

%  v1 June 16, 2022, Modified: initial version
%  v2 June 17, 2022, Modified: removed unnecessary code
%  v3 June 20, 2022, Modified: modeled 2 blocks: REST/GAME/REST/GAME etc...
%  v4 June 21, 2022, Modified: modeled only 40% and 0% difficulty levels
%  v5 June 23, 2022, Modified: modeled sham real-time NF-based control of
%  the difficulty (e.g. 40% if NF is low, 0% if NF is high)
%  v6 Jully 25, 2022, Modified: added variable to save frame-wise diff
%  level
%  v7 Nov 8, 2022, Modified: auto speed adjustment & diff level set to 30%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%   BUTTON CONTROLS  %%%%%%%%%%%%%%%%%%%% 
%    -  START GAME: modeled as a 5 keypress
%
%    -  L/R MOVEMENT: keys "K" and "L" for left/right for R-handed
%                     keys "A" and "S" for left/right for L-handed
%
%    -  CONFIRMATION: "space" key
%
%    -  EMERGENCY ESCAPE: mouse click (will break the game)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%   DEBUGGING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMINDER: how to find keycode on any operation machine
% KbName: get the symbol e.g. KbName('a')
% KbName(SYMBOL): get the keycode number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Clear and setup params
close all 
clear
clc
commandwindow;

% skip sync tests
Screen('Preference', 'SkipSyncTests', 1);

%% Screen parameters


% For other setups
tmp = Screen('Resolution', 1);
resolution = [tmp.width,tmp.height]; % standard settings for other machines 


%% Request inputs
prompt = 'Operating machine (w for windows or i for i0S): ';
machineCode = input(prompt,'s');
clear prompt

prompt = 'Participant''s code: ';
participantCode = input(prompt,'s');
clear prompt
handedness = menu('Play with:', 'Right hand','Left hand');

randIndex = [2]; % 1 rest, 2 game


%% Initialization

% Initialize the manipulations 

randomCursor0 =  [0 0 0 0 0 0 0 0 0];
randomCursor33 = [-1 -1 0 0 0 0 0 0 1];
randomCursorALL = [randomCursor0; randomCursor33];


%% Lag settings (will not be used at the end)
lagDuration = 0.25; % delay in seconds (if nPhases = 5, lag will not be implemented).

%% Graphic elements
red = [255 50 0];
green = [0 240 0];
blue = [0 160 255];
lightblue = [135 206 235];
white = [255 255 255];
dark = [52 52 52];
yellow = [255 255 0];
cursorColor = white;

%% Set stimuli parameters
stimuli.nStimuli = 4*12;                      % number of stimuli (X)
distractors.nDistractors = 4*12;              % number of distractors (O)
stimuli.size = round(resolution(1)*0.009);                   % size of stimuli (pixels)
stimuli.center = [0,0];           % center of the field of stimuli (x,y)
stimuli.apertureSize = [round(resolution(1)*0.43478),4*round(resolution(2)*1.7)];  % size of rectangular aperture [w,h] (game space boundaries)
stimuli.speed = round(resolution(2)*0.002);       % pixels/frame
stimuli.duration = 4*15;       % duration of block in sec


% Define boundaries in pixels
lPix = stimuli.center(1) - stimuli.apertureSize(1)/2 + resolution(1)/2 - 10;
rPix = stimuli.center(1) + stimuli.apertureSize(1)/2 + resolution(1)/2 + 10;

performance = [];
timing = [];


%% Assign location of stimuli and distractors

% Generate positions on the grid
nSectors = 7;

xPositionGrid = zeros(nSectors,1);
for i = 1:nSectors
    xPositionGrid(i) = round((stimuli.apertureSize(1)/nSectors*i) + stimuli.center(1) - (stimuli.apertureSize(1)/(2*nSectors)) - stimuli.apertureSize(1)/2);
end

yPositionGrid = zeros(stimuli.nStimuli+distractors.nDistractors,1);
for i = 1:stimuli.nStimuli+distractors.nDistractors
    yPositionGrid(i) = round((stimuli.apertureSize(2)/(stimuli.nStimuli+distractors.nDistractors)*i) + (stimuli.center(2) - round(stimuli.apertureSize(2)*1)));
end


%% Open the screen
try
    % Open the screen
    screens = Screen('Screens');
    screenNumber = max(screens);
    Screen('Preference', 'SkipSyncTests', 1);
    [w, ~] = Screen('OpenWindow', screenNumber, [0 0 0]);
    
    % added by Serafeim
    Priority(MaxPriority(w));

    % Get the framerate
    fps = Screen('FrameRate',w);      % frames per second
    
    nFrames = secs2frames(fps, stimuli.duration);   % total number of frames before it stops
    stimuli.speed = (stimuli.apertureSize(2)+(resolution(2)/2))/nFrames;

    % Set the text parameters
    Screen('TextFont',w, 'Arial');
    Screen('TextStyle', w, 1);
        
    randTurbulence = zeros(nFrames,1);
    
    % Display fixation cross and Hide cursor
    %HideCursor;

    Screen('DrawLine', w, white, resolution(1)/2-round(resolution(1)*0.03125), resolution(2)/2, resolution(1)/2+round(resolution(1)*0.03125), resolution(2)/2, round(stimuli.size/2))
    Screen('DrawLine', w, white, resolution(1)/2, resolution(2)/2-round(resolution(1)*0.03125), resolution(1)/2, resolution(2)/2+round(resolution(1)*0.03125), round(stimuli.size/2))
    Screen('Flip', w);

    
    %% Play the game

    KbQueueCreate;
    KbQueueStart;
 
    % Wait for the start signal
    % Verify: keys_all = KbName('KeyNamesWindows'); keys_all{48}
    confirmation = 0;
    while confirmation == 0    
        % Get input
        [pressed, key] = KbQueueCheck;
        if pressed
            ch = find(key);
            if startsWith(machineCode, 'w') || startsWith(machineCode, 'W')
                keytarget_trigger = 53; % windows
            elseif startsWith(machineCode, 'i') || startsWith(machineCode, 'I')
                keytarget_trigger = 34; % iOS
            end
            if ch == keytarget_trigger % KbName('5%') == 34 for mac (iOs)
                KbQueueFlush;
                confirmation = 1;
            end
        end
    end
    KbQueueFlush;
    clear keytarget_trigger
    
    tStartGame1 = clock;
    turb_counter = 0;
    for k = 1:length(randIndex)
    
        switch randIndex(k)
            case 1      % REST
                tLag = 0;
                turbulence = 0;
                magic = 0;
                gamePhase = randIndex(k);
                %frameColor = green;
                frameColor = yellow;
            case 2      % GAME AT 40%
                tLag = 0;
                turbulence = 1;
                magic = 0;
                gamePhase = randIndex(k);
                %frameColor = green;
                frameColor = lightblue;

        end

        tLagFrames = secs2frames(fps, tLag); assert(tLagFrames == 0);

        disp(randIndex(k))
        % if REST dispay just a fixation cross

        if randIndex(k)==1 % if REST
            WaitSecs(1.5);
            tStartGame = clock;
            deltaTGame = 0;
            
            %for i = 2+tLagFrames:nFrames+tLagFrames
            for i = 2+tLagFrames:(nFrames+tLagFrames)/2
                Screen('DrawLine', w, frameColor, resolution(1)/2-round(resolution(1)*0.03125), resolution(2)/2, resolution(1)/2+round(resolution(1)*0.03125), resolution(2)/2, round(stimuli.size/2))
                Screen('DrawLine', w, frameColor, resolution(1)/2, resolution(2)/2-round(resolution(1)*0.03125), resolution(1)/2, resolution(2)/2+round(resolution(1)*0.03125), round(stimuli.size/2))
                deltaTGame = etime(clock,tStartGame);
                Screen('Flip', w);
            end
        
        tEndGame = clock;
        
        performance = horzcat(performance, [gamePhase; NaN; NaN; deltaTGame; NaN]); 
        timing = vertcat(timing, [tStartGame(4:end) tEndGame(4:end)]);  
        
        else

        
        distStimuli = zeros(2,1);
        randTurbulence = zeros(2,1);
        
        % Create stepsize
        cursorStepsize = round(stimuli.apertureSize(1)/nSectors);
        
        % Create random position for stimuli and distractors
        randPositionX = randi([1,nSectors],stimuli.nStimuli+distractors.nDistractors,1);
        xPositionGrid1 = xPositionGrid(randPositionX);

        randPositionY = randperm(stimuli.nStimuli+distractors.nDistractors);
        yPositionGrid1 = yPositionGrid(randPositionY);

        % Separate stimuli and distractors position
        stimuli.x = xPositionGrid1(1:stimuli.nStimuli);
        stimuli.y = yPositionGrid1(1:stimuli.nStimuli);

        distractors.x = xPositionGrid1(stimuli.nStimuli+1:end);
        distractors.y = yPositionGrid1(stimuli.nStimuli+1:end);
        
        clear randPositionX randPositionY i
        
        % Check inter-stimuli distance
        [stimuli.y, sortedIndex] = sort(stimuli.y);
        stimuli.x = stimuli.x(sortedIndex);
        for j = 2:length(stimuli.x)
            if pdist2(stimuli.x(j),stimuli.x(j-1)) >= pdist2(xPositionGrid(1),xPositionGrid(end-1))     % inter-stimuli distance not larger than 5 steps (with 7 sectors)
                if stimuli.x(j) > 0
                    stimuli.x(j) = xPositionGrid(end-2);
                elseif stimuli.x(j) < 0
                    stimuli.x(j) = xPositionGrid(3);
                end
            end
        end
        clear j
        

        %% Main loop
        
        % Display fixation cross
        Screen('DrawLine', w, frameColor, resolution(1)/2-round(resolution(1)*0.03125), resolution(2)/2, resolution(1)/2+round(resolution(1)*0.03125), resolution(2)/2, round(stimuli.size/2))
        Screen('DrawLine', w, frameColor, resolution(1)/2, resolution(2)/2-round(resolution(1)*0.03125), resolution(1)/2, resolution(2)/2+round(resolution(1)*0.03125), round(stimuli.size/2))
        Screen('Flip', w);
        
        Screen('TextSize',w, round(resolution(1)*0.025)); % Set the font size to draw the X
        WaitSecs(1.5);
        
        % Number of hits is initially void
        nHits = zeros(stimuli.nStimuli,1);
        nDistractors = zeros(distractors.nDistractors,1);
        hitX = 0;
        hitDistractors = 0;
        
        tStartGame = clock;
        deltaTGame = 0;
        cursorPosition = ones(nFrames+tLagFrames,1)*resolution(1)/2;

        mycounter = 0;

        all_diff_lvl_frames = zeros(1,length(2+tLagFrames:nFrames+tLagFrames));
        % Main loop for frames
        for i = 2+tLagFrames:nFrames+tLagFrames % 900 frames = 60 fps * duration (i.e., 15 secs game)

            % center the stimuli
            pixposStimuli.x = stimuli.x + resolution(1)/2;
            pixposStimuli.y = stimuli.y + resolution(2)/2;

            pixposDistractors.x = distractors.x + resolution(1)/2;
            pixposDistractors.y = distractors.y + resolution(2)/2;

            %% NF-based difficulty adjustment (within the frames loop!)
            % If NF signal is above the threshold then change the
            % difficulty accordingly
            mycounter = mycounter+1;
            NF_value = get_NF_value(mycounter, nFrames);
            % disp(NF_value)
            if NF_value <= 50 % make it difficult
                diff_level = 2;
                randomCursor = randomCursorALL(diff_level,:);
                randTurbIndex = turbulence*randi([1 length(randomCursor)],1,1);
    
            elseif NF_value > 50
                diff_level = 1;
                randomCursor = randomCursorALL(diff_level,:);
                randTurbIndex = 3; % if no turbulence, set index to 3 so that randomCursor(randTurbIndex) == 0
            end
            randTurbulence(i) = randomCursor(randTurbIndex);
            % store frame-wise difficulty level
            all_diff_lvl_frames(i-1) = diff_level;
            
            %% Get input
            [pressed, key] = KbQueueCheck;
            if pressed
                ch = find(key);
                switch handedness       % The response buttons are different if you play with right or left hand
                    
                    case 1 % right-handed
                        if startsWith(machineCode, 'w') || startsWith(machineCode, 'W')
                            keytarget_left =  75; % "k" in windows
                            keytarget_right = 76; % "l" in windows
                        elseif startsWith(machineCode, 'i') || startsWith(machineCode, 'I')
                            keytarget_left =  14; % "k" in iOS
                            keytarget_right = 15; % "l" in iOS
                        end
                        if ch == keytarget_left % go left i.e. KbName(31) = '2@'
                            KbQueueFlush;
                            cursorPosition(i) = cursorPosition(i-1)-cursorStepsize+(randTurbulence(i)*cursorStepsize);
                        elseif ch == keytarget_right % go right
                            KbQueueFlush;
                            cursorPosition(i) = cursorPosition(i-1)+cursorStepsize+(randTurbulence(i)*cursorStepsize);
                        else
                            cursorPosition(i) = cursorPosition(i-1);
                        end
                        
                    case 2 % left-handed
                        if startsWith(machineCode, 'w') || startsWith(machineCode, 'W')
                            keytarget_left =  65; % "a" in windows
                            keytarget_right = 83; % "s" in windows
                        elseif startsWith(machineCode, 'i') || startsWith(machineCode, 'I')
                            keytarget_left =  4; % "a" in iOS
                            keytarget_right = 22; % "s" in iOS
                        end
                        if ch == keytarget_left  % go left
                            KbQueueFlush;
                            cursorPosition(i) = cursorPosition(i-1)-cursorStepsize+(randTurbulence(i)*cursorStepsize);
                        elseif ch == keytarget_right  % go right
                            KbQueueFlush;
                            cursorPosition(i) = cursorPosition(i-1)+cursorStepsize+(randTurbulence(i)*cursorStepsize);
                        else
                            cursorPosition(i) = cursorPosition(i-1);
                        end
                end
            else
                KbQueueFlush;
                cursorPosition(i) = cursorPosition(i-1);
            end
            KbQueueFlush;

            % Set boundaries for the cursor
            if cursorPosition(i) > rPix-round(resolution(1)*0.01367)        % Set boundaries for cursor position
                cursorPosition(i) = max(xPositionGrid) + resolution(1)/2;         
            elseif cursorPosition(i) < lPix+round(resolution(1)*0.01367)
                cursorPosition(i) = min(xPositionGrid) + resolution(1)/2;
            end
            
            % Display the stimuli (X)
            for j = 1:stimuli.nStimuli          % A X is created with two lines crossing each other
                if find(j == hitX) > 0      % Update color for hit X
                    color = green;
                else 
                    color = white;
                end
                Screen('DrawText', w, 'X', pixposStimuli.x(j)-round(resolution(1)*0.0104), pixposStimuli.y(j)-round(resolution(1)*0.0104), color);    % A X is actually a x
            end

            % Display the distractors (O)
            for j = 1:distractors.nDistractors         
                if find(j == hitDistractors) > 0      % Update color for hit X
                    color = red;
                else 
                    color = white;
                end
                Screen('FrameOval', w, color, [pixposDistractors.x(j)-stimuli.size; pixposDistractors.y(j)-stimuli.size; pixposDistractors.x(j)+stimuli.size; pixposDistractors.y(j)+stimuli.size], round(stimuli.size/2))
            end

            Screen('DrawLine', w, frameColor, lPix, 1, lPix, resolution(2), 5)
            Screen('DrawLine', w, frameColor, rPix, 1, rPix, resolution(2), 5)
            Screen('FillRect', w, frameColor, [lPix round(resolution(2)*0.7604) rPix round(resolution(2)*0.802)])

            if randIndex(k) ~= 1    % if not PureVisual
            
                % Display cursor
                Screen('FillRect', w, cursorColor, [cursorPosition(i-tLagFrames)-round(resolution(1)*0.01367) round(resolution(2)*0.75694) cursorPosition(i-tLagFrames)+round(resolution(1)*0.01367) round(resolution(2)*0.8055)])

                % Check for hits
                for z = 1:stimuli.nStimuli
                    distStimuli = sqrt((pixposStimuli.x(z)-cursorPosition(i-tLagFrames)-randTurbulence(i)).^2 + (pixposStimuli.y(z)-round(resolution(2)*0.78125)).^2);  % Compute distances between each point
                    if magic == 0       % Implement the magic condition
                        if min(distStimuli) < round(resolution(1)*0.014)
                            nHits(z) = z;
                        end
                    else
                        if min(distStimuli) < magic
                            nHits(z) = z;
                        end
                    end
                end
                hitX = find(nHits ~= 0);

                % Check for distractors
                for z = 1:distractors.nDistractors
                    distDistractors = sqrt((pixposDistractors.x(z)-cursorPosition(i-tLagFrames)-randTurbulence(i)).^2 + (pixposDistractors.y(z)-round(resolution(2)*0.78125)).^2);  % Compute distances between each point

                    if min(distDistractors) < round(resolution(1)*0.014)
                        nDistractors(z) = z;
                    end
                end
                hitDistractors = find(nDistractors ~= 0);
            end
            
            % update the dot position
            stimuli.y = stimuli.y + stimuli.speed;
            distractors.y = distractors.y + stimuli.speed;

            deltaTGame = etime(clock,tStartGame);
            Screen('Flip', w);
        end
        tEndGame = clock;

        
        % Save the data
        performance = horzcat(performance, [gamePhase; length(hitX); length(hitDistractors); deltaTGame; diff_level]); 
        timing = vertcat(timing, [tStartGame(4:end) tEndGame(4:end)]);

        end
    end
    deltaTTotal = etime(clock,tStartGame1);
    tEndTotal = clock;
    
    KbQueueStop;
    KbQueueRelease;


    % Closing screen (end of experiment)
    Screen('TextSize',w, round(resolution(1)*0.04)); % Set the font size to draw the X
    Screen('DrawText', w, 'Game end', resolution(1)/2-round(resolution(1)*0.09), resolution(2)/2+round(resolution(1)*0.001), [255 255 255]);
    Screen('Flip', w);
    %GetClicks; 
    pause(1);
    


catch ME
    Screen('CloseAll');
    rethrow(ME)
end


Screen('CloseAll');


%% Save the results on a .txt file

dateToday = datestr(now, 'HHMMSS');

save(horzcat(participantCode,'_Results_',dateToday))

fileID = fopen(horzcat(participantCode,'_Results_',dateToday,'.txt'),'wt');
fprintf(fileID,'%s\t %s\t %s\t %s\t %s\r\n','GamePhases','HitX','HitO', 'deltaTGame', 'diff_level');  % added by Serafeim
fprintf(fileID,'%i\t %i\t %i\t %d\t %d\t %d\t %i\r\n', performance);  % added by Serafeim
fclose(fileID);

timing = timing';
fileID = fopen(horzcat(participantCode,'_Timing_',dateToday,'.txt'),'w');
fprintf(fileID,'%s\t %s\t %s\t %s\t %s\t %s\t\n', 'tGameStart(hh)', 'tGameStart(mm)', 'tGameStart(ss)', 'tGameEnd(hh)', 'tGameEnd(mm)', 'tGameEnd(ss)');
fprintf(fileID,'%i\t %i\t %d\t %i\t %i\t %d\t\n', timing);
fclose(fileID);

fileID = fopen(horzcat(participantCode,'_ResultsFrames_',dateToday,'.txt'),'wt');
fprintf(fileID,'%s\t\n','Diff_frames');  % added by Serafeim
fprintf(fileID,'%d\r\n', all_diff_lvl_frames);  % added by Serafeim
fclose(fileID);
% clear s nHits nDistractors b l r t w lPix rPix screenNumber screens randTurbulence color confirmation cursorColor cursorPosition newCursorPosition ...
%     gamePhase i hitDistractors hitX index j key pressed ch cursorStepsize k subjectiveAgency subjectivePerformance tempo magicGame z turbulence magic1 ...
%     turbulence1 xPositionGrid xPositionGrid1 yPositionGrid yPositionGrid1 tStartGame tEndGame tStartPerformance tEndPerformance tStartAgency tEndAgency resolution tStartgame1 fileID

