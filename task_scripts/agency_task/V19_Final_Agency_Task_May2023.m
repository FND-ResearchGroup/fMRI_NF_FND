%% Agency task 
% Script to the fMRI agency task used in the NF FND project. 

%  VERSION LOG
%  Initial code by Giuseppe A. Zito, 2019. 
%  v16 Oct 26th, 2021, Modified by Serafeim Loukas, 
%  v17 Sep 19th, 2022, Modified by Serafeim et al., new button numbering
%  tested on Sep 20, 2022 by Serafeim
%  v18 May 31st, 2023, Modified by Nicolas Gninenko,
%   fix of text (German version) + number of baseline vs turb trials, up to 8.6 min
%  v19 General bug fixes (resolution) and annotations, Nicolas Gninenko & Eliane Müller

% BUTTON CONTROLS
% MRI trigger: modeled as a 5 keypress
% box buttons 1/3 if R hand, 6/8 if L hand
% 2 if R hand or 7 if L hand
% escape: abort

% The script is written for a total of 3 conditions: 
% - Pure visual (PV)
% - Baseline (BH & BE), also called Non-Turbulence
% - Turbulence (TH & TE)

% Each condition consists of 3 phases: 
% 1. Gamephase
% 2. Judgement of performance
% 3. Judgement of agency

close all; tic;
clear;
clc; sca;
commandwindow;

% skip sync tests
Screen('Preference', 'SkipSyncTests', 1);

% Screen parameters
tmp = Screen('Resolution',1); % use external monitor for testing
resolution = [tmp.width,tmp.height];


%% Request input
% Format is agency_PXX_VYY.m
participantCode=mfilename; 
participantCode=participantCode(regexp(participantCode,'P*'):regexp(participantCode,'V*')-2);
%clear prompt % add path to Templates used
addpath('.\230000_Templates'); 
handedness = menu('Play with:', 'Right hand','Left hand');

%% Randomization of game phases

% Define phases 
% In this version
% BE & BH  = Basline, TE & TH = Turbulence. There is
% no difference between BE & BH, or TE or TH. 
nPhases = 5;     %(PV, BE, BH, TE, TH)       
nRepetitions = 4;   % n of repetition of the full game (nDifficulties x nPhases)


randOrder = randperm(26);
randIndex = [1 1 1 1 2*ones(1,11) 4*ones(1,11)];
randIndex = randIndex(randOrder);

clear nPhases nonRandIndex;


%% Initialization

% Initialize the manipulations 
% 40% Turbulence (in 40% of button clicks cursor movment is manipulated)
randomCursor = [-1 -1 0 0 0 0 0 0 1 1];     % positions for turbulence phase
lagDuration = 0.25;                       % delay in seconds

% Graphic elements
red = [255 50 0];
green = [0 240 0];
blue = [0 160 255];
white = [255 255 255];
lightblue = [135 206 235];
cursorColor = white;

% Set stimuli parameters
stimuli.nStimuli = 12;                % number of stimuli (X)
distractors.nDistractors = 12;              % number of distractors (O)
stimuli.size = round(resolution(1)*0.009);                   % size of stimuli (pixels)
stimuli.center = [0,0];           % center of the field of stimuli (x,y)
stimuli.apertureSize = [round(resolution(1)*0.43478),round(resolution(2)*1.7)];     % size of rectangular aperture [w,h] (game space boundaries)
stimuli.speed = 2.64;%round(resolution(2)*0.002);       % pixels/frame
stimuli.duration = 15;       % duration of block in sec

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
    %screenNumber = max(screens);
    Screen('Preference', 'SkipSyncTests', 1);
    [w, ~] = Screen('OpenWindow', 1, [0 0 0]);

    % Get the framerate
    fps = Screen('FrameRate', w);      % frames per second
    
    nFrames = secs2frames(fps, stimuli.duration);   % total number of frames before it stops
    
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
    
    % Wait for the scanner to send the trigger (corresponding to the 5 keypress)
    % Verify: keys_all = KbName('KeyNamesWindows'); keys_all{48}
    confirmation = 0;
    while confirmation == 0    
        % Get input
        [pressed, key] = KbQueueCheck;
        if pressed
            ch = find(key);
            if ch == 53 % KbName('5%') == 53 in windows
                KbQueueFlush;
                confirmation = 1;
            end
        end
    end
    KbQueueFlush;
    
    tStartGame1 = clock;
    for k = 1:length(randIndex)
    
        switch randIndex(k)
            case 1      % Pure Visual
                tLag = 0;
                turbulence = 0;
                magic = 0;
                gamePhase = randIndex(k);
                frameColor = white;
            case 2      % Baseline Easy
                tLag = 0;
                turbulence = 0;
                magic = 0;
                gamePhase = randIndex(k);
                %frameColor = green;
                frameColor = lightblue;
            case 3      % Baseline Hard
                tLag = 0;
                turbulence = 0;
                magic = 0;
                gamePhase = randIndex(k);
                frameColor = lightblue;
                %frameColor = red;
            case 4      % Turbulence Easy
                tLag = 0;
                turbulence = 1;
                magic = 0;
                gamePhase = randIndex(k);
                frameColor = lightblue;
                %frameColor = green;
            case 5      % Turbulence Hard
                tLag = 0;
                turbulence = 1;
                magic = 0;
                gamePhase = randIndex(k);
                %frameColor = red;
                frameColor = lightblue;

        end
        tLagFrames = secs2frames(fps, tLag);
        
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
        
        Screen('TextSize',w, round(resolution(1)*0.025));       % Set the font size to draw the X
        WaitSecs(1.5);
        
        % Number of hits is initially void
        nHits = zeros(stimuli.nStimuli,1);
        nDistractors = zeros(distractors.nDistractors,1);
        hitX = 0;
        hitDistractors = 0;
        
        tStartGame = clock;
        deltaTGame = 0;
        cursorPosition = ones(nFrames+tLagFrames,1)*resolution(1)/2;

        for i = 2+tLagFrames:nFrames+tLagFrames

            % center the stimuli
            pixposStimuli.x = stimuli.x + resolution(1)/2;
            pixposStimuli.y = stimuli.y + resolution(2)/2;

            pixposDistractors.x = distractors.x + resolution(1)/2;
            pixposDistractors.y = distractors.y + resolution(2)/2;

            % Implement the turbulence phase
            if turbulence == 1
                randTurbIndex = turbulence*randi([1 length(randomCursor)],1,1);
            else
                randTurbIndex = 3;
            end
            randTurbulence(i) = randomCursor(randTurbIndex);

            % Get input
            [pressed, key] = KbQueueCheck;
            if pressed
                ch = find(key);
                switch handedness       % The response buttons are different if you play with right or left hand
                    case 1
                        if ch == 49 % if 1!
                        %if ch == 31
                            KbQueueFlush;
                            cursorPosition(i) = cursorPosition(i-1)-cursorStepsize+(randTurbulence(i)*cursorStepsize);
                        elseif ch == 51 % if 3#
                        %elseif ch == 32
                            KbQueueFlush;
                            cursorPosition(i) = cursorPosition(i-1)+cursorStepsize+(randTurbulence(i)*cursorStepsize);
                        else
                            cursorPosition(i) = cursorPosition(i-1);
                        end
                    case 2
                        if ch == 56 % if 8*
                        %if ch == 37
                            KbQueueFlush;
                            cursorPosition(i) = cursorPosition(i-1)-cursorStepsize+(randTurbulence(i)*cursorStepsize);
                        elseif ch == 54 % if 6^
                        %elseif ch == 36
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


        %% Judgment of performance
        
        tStartPerformance = clock;
        if randIndex(k) == 1
            subjectivePerformance = NaN;
            hitX = 0;
            hitDistractors = 0;
            exitSignal = 0;
        else            
            newCursorPosition = resolution(1)/2;
            confirmation = 0;
            while confirmation == 0    
                % Get input
                [pressed, key] = KbQueueCheck;
                if pressed
                    ch = find(key);
                    switch handedness
                        case 1
                            if ch == 49 % keys_all = KbName('KeyNamesWindows'); keys_all{49} = '1!' in windows
                                KbQueueFlush;
                                newCursorPosition = newCursorPosition-round((rPix-lPix)/11);
                            elseif ch == 51 % if 3#
                                KbQueueFlush;
                                newCursorPosition = newCursorPosition+round((rPix-lPix)/11);
                            elseif ch == 50 % if 2@
                                KbQueueFlush;
                                confirmation = 1;
                                subjectivePerformance = newCursorPosition;
                            end
                        case 2
                            if ch == 56 % if 8*
                                KbQueueFlush;
                                newCursorPosition = newCursorPosition-round((rPix-lPix)/11);
                            elseif ch == 54 % if 6^
                                KbQueueFlush;
                                newCursorPosition = newCursorPosition+round((rPix-lPix)/11);
                            elseif ch == 55 % if 7&
                                KbQueueFlush;
                                confirmation = 1;
                                subjectivePerformance = newCursorPosition;
                            end
                    end
                end
                KbQueueFlush;

                % Set boundaries for the cursor
                if newCursorPosition > rPix-round(resolution(1)*0.0136718)        % Set boundaries for cursor position
                    newCursorPosition = resolution(1)/2+(5*round((rPix-lPix)/11));
                elseif newCursorPosition < lPix+round(resolution(1)*0.0136718)
                    newCursorPosition = resolution(1)/2-(5*round((rPix-lPix)/11));
                end

                % Update position of cursor
                Screen('TextSize',w, round(resolution(1)*0.025));               
                Screen('DrawText', w, 'PERFORMANCE', resolution(1)/2-round(resolution(1)*0.13), resolution(2)/2+round(resolution(1)*0.07), [255 255 255]);
                Screen('DrawText', w, 'Very high', rPix+30, 800, [255 255 255]);
                Screen('DrawText', w, 'Very poor', lPix-330, 800, [255 255 255]); %round(resolution(1)*0.125) ; round(resolution(2)*0.765)
                Screen('FillRect', w, [255 255 255], [lPix round(resolution(2)*0.7604) rPix round(resolution(2)*0.802)])
                Screen('FillRect', w, [130 130 130], [newCursorPosition-round(resolution(1)*0.013672) round(resolution(2)*0.75694) newCursorPosition+round(resolution(1)*0.013672) round(resolution(2)*0.8055)])
                Screen('Flip', w);

                % Get mouse stop from the operator
                exitSignal = 0;
                [~,~,buttons] = GetMouse;
                if any(buttons)
                    buttons = [];
                    exitSignal = 1;
                    subjectivePerformance = 0;
                    break;
                end
            end
        end
        deltaTPerformance = etime(clock,tStartPerformance);
        tEndPerformance = clock;

        %% Judgment of agency
        
        tStartAgency = clock;
        if randIndex(k) == 1
            subjectiveAgency = NaN;
            hitX = 0;
            hitDistractors = 0;
            exitSignal = 0;
        else            
        
            newCursorPosition = resolution(1)/2;
            confirmation = 0;
            while confirmation == 0    
                % Get input
                [pressed, key] = KbQueueCheck;
                if pressed
                    ch = find(key);
                    switch handedness
                        case 1
                            if ch == 49
                            %if ch == 30
                                KbQueueFlush;
                                newCursorPosition = newCursorPosition-round((rPix-lPix)/11);
                            elseif ch == 51
                            %elseif ch == 32
                                KbQueueFlush;
                                newCursorPosition = newCursorPosition+round((rPix-lPix)/11);
                            elseif ch == 50
                            %elseif ch == 31
                                KbQueueFlush;
                                confirmation = 1;
                                subjectiveAgency = newCursorPosition;
                            end
                        case 2
                            if ch == 56
                            %if ch == 37
                                KbQueueFlush;
                                newCursorPosition = newCursorPosition-round((rPix-lPix)/11);
                            elseif ch == 54
                            %elseif ch == 35
                                KbQueueFlush;
                                newCursorPosition = newCursorPosition+round((rPix-lPix)/11);
                            elseif ch == 55
                            %elseif ch == 36
                                KbQueueFlush;
                                confirmation = 1;
                                subjectiveAgency = newCursorPosition;
                            end
                    end
                end
                KbQueueFlush;

                % Set boundaries for the cursor
                if newCursorPosition > rPix-round(resolution(1)*0.013672)        % Set boundaries for cursor position
                    newCursorPosition = resolution(1)/2+(5*round((rPix-lPix)/11));
                elseif newCursorPosition < lPix+round(resolution(1)*0.013672)
                    newCursorPosition = resolution(1)/2-(5*round((rPix-lPix)/11));
                end

                % Update position of cursor
                Screen('DrawText', w, 'CONTROL', resolution(1)/2-round(resolution(1)*0.08), resolution(2)/2+round(resolution(1)*0.07), [255 255 255]);
                Screen('DrawText', w, 'Very high', rPix+30, 800, [255 255 255]);
                Screen('DrawText', w, 'Very low', lPix-305, 800, [255 255 255]);
                Screen('FillRect', w, [255 255 255], [lPix round(resolution(2)*0.7604) rPix round(resolution(2)*0.802)])
                Screen('FillRect', w, [130 130 130], [newCursorPosition-round(resolution(1)*0.013672) round(resolution(2)*0.75694) newCursorPosition+round(resolution(1)*0.013672) round(resolution(2)*0.8055)])
                Screen('Flip', w);

                % Get mouse stop from the operator
                [~,~,buttons] = GetMouse;
                if any(buttons)
                    buttons = [];
                    exitSignal = 1;
                    subjectiveAgency = 0;
                    break;
                end            
            end

            % Exit if mouse pressed by operator
            if exitSignal == 1
                break
            end
        end
        deltaTAgency = etime(clock,tStartAgency);
        tEndAgency = clock;
        
        % Save the data
        performance = horzcat(performance, [gamePhase; length(hitX); length(hitDistractors); subjectivePerformance-resolution(1)/2; subjectiveAgency-resolution(1)/2; deltaTGame; deltaTPerformance; deltaTAgency]);
        timing = vertcat(timing, [tStartGame(4:end) tEndGame(4:end) tStartPerformance(4:end) tEndPerformance(4:end) tStartAgency(4:end) tEndAgency(4:end)]);
    end
    deltaTTotal = etime(clock,tStartGame1);
    tEndTotal = clock;
    
    KbQueueStop;
    KbQueueRelease;

    % Closing screen
    Screen('TextSize',w, round(resolution(1)*0.04));       % Set the font size to draw the X
   % Screen('DrawText', w, 'Game end', resolution(1)/2-round(resolution(1)*0.06), resolution(2)/2+round(resolution(1)*0.078125), [255 255 255]);
    Screen('DrawText', w, 'Game end', resolution(1)/2-round(resolution(1)*0.12), resolution(2)/2-round(resolution(1)*0.02), [255 255 255]);
    Screen('Flip', w);
    GetClicks;  
    
catch ME
    Screen('CloseAll');
    rethrow(ME)
end
Screen('CloseAll');


%% Save the results on a .txt file

dateToday = datestr(now, 'HHMMSS');

save([pwd filesep horzcat(participantCode,'_Results_',dateToday,'.mat')],...
    'blue','buttons','ch','color','confirmation','cursor*','dateToday','delta*',...
    'dist*','exitSignal','fps','frameColor','gamePhase','green','handedness',...
    'hitDistractors','hitX','key','lagDuration','i','j','k','lightblue','lPix',...
    'magic','nDistractors','newCursorPosition','nFrames','nHits','nRepetitions',...
    'nSectors','participantCode','performance','pixpos*','pressed','rand*','red',...
    'resolution','rPix','screens','sortedIndex','stimuli',...
    'subjectiveAgency','subjectivePerformance','tEnd*','timing','tLag','tLagFrames',...
    'tmp','tStart*','turbulence','w','white','xPositionGrid*','yPositionGrid*','z');

fileID = fopen([pwd filesep horzcat(participantCode,'_Results_',dateToday,'.txt')],'wt');
fprintf(fileID,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\r\n','GamePhases','HitX','HitO','JoP','JoA', 'deltaTGame', 'deltaTPerformance', 'deltaTAgency');
fprintf(fileID,'%i\t %i\t %i\t %i\t %i\t %d\t %d\t %d\r\n', performance);
fclose(fileID);

timing = timing';
fileID = fopen([pwd filesep horzcat(participantCode,'_Timing_',dateToday,'.txt')],'w');
fprintf(fileID,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\r\n', 'tGameStart(hh)', 'tGameStart(mm)', 'tGameStart(ss)', 'tGameEnd(hh)', 'tGameEnd(mm)', 'tGameEnd(ss)', ...
    'tJoPStart(hh)', 'tJoPStart(mm)', 'tJoPStart(ss)', 'tJoPEnd(hh)', 'tJoPEnd(mm)', 'tJoPEnd(ss)', 'tJoAStart(hh)', 'tJoAStart(mm)', 'tJoAStart(ss)', 'tJoAEnd(hh)', 'tJoAEnd(mm)', 'tJoAEnd(ss)');
fprintf(fileID,'%i\t %i\t %d\t %i\t %i\t %d\t %i\t %i\t %d\t %i\t %i\t %d\t %i\t %i\t %d\t %i\t %i\t %d\r\n', timing);
fclose(fileID);

clear s nHits nDistractors b l r t w lPix rPix screenNumber screens randTurbulence color confirmation cursorColor cursorPosition newCursorPosition ...
     gamePhase i hitDistractors hitX index j key pressed ch cursorStepsize k subjectiveAgency subjectivePerformance tempo magicGame z turbulence magic1 ...
     turbulence1 xPositionGrid xPositionGrid1 yPositionGrid yPositionGrid1 tStartGame tEndGame tStartPerformance tEndPerformance tStartAgency tEndAgency resolution tStartgame1 fileID ...
     white timing performance participantCode ans blue dateToday distractors fps green handedness lagDuration lightblue nFrames nRepetitions nSectors ...
     randIndex randomCursor randOrder red stimuli tmp
 
 fprintf(['\nElapsed ' num2str(toc) ' sec...\n']);
