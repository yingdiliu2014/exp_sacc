% express saccade project once and for all
% Yingdi LIU, 2014,12,20, Fribourg 

%% Routines 

clear all; clc; close all;
% add the current directory (and its subfolders) to Matlab path and cd to
% the input directory
myFolderRoutine('C:\Users\Liuy\Google\Google Drive\express_saccades') 


%% basic info 

run adoptee_names.m % get the names of the participants 
load('allRT.mat'); % reaction time of each of the 12 participants. 


%% load eye data 

eyedataAll = cell(12,1);

for subjNum = 1:12 
    subjname = ['at',num2str(subjNum)];
    suffix = '_APES_eye.txt';
    eyedata = load([subjname,suffix]);
    
    for cc = 1:2 % 1 for gap, 2 for overlap
        conditionData = eyedata(eyedata(:,2)==cc,:);
        
        trials = conditionData(:,3);
        directions = conditionData(:,4); 
        [trialStartIdx, trialEndIdx] = separateTrials(trials,3);
    
    
end
    
    
    
    