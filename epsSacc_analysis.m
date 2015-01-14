%% analyze express saccade data

%% data structure
% eyedatatowrite=[blk blkType tr direction period 0 echantillon eyex1 eyey1 eyex2 eyey2 trialstart trialstart currentTime]; % 0 stands for 'data before trial'
% eyedatatowrite=[blk blkType tr direction period 1 echantillon eyex1 eyey1 eyex2 eyey2 trialstart oldTime currentTime]; % 1 stands for fixation period;
% eyedatatowrite=[blk blkType tr direction period 2 echantillon eyex1 eyey1 eyex2 eyey2 trialstart gapTime currentTime]; % 2 stands for gap (or the last 0.2s for overlap);
% eyedatatowrite=[blk blkType tr direction period 3 echantillon eyex1 eyey1 eyex2 eyey2 trialstart targetTime currentTime]; % 3 for target period;


%%

clear all; clc; close all;
cd '/Users/yingdiliu/Drive/express_saccades'
addpath(genpath('/Users/yingdiliu/Drive/express_saccades'))


%% 1.
load('allRT.mat');
allRTnew = cell(12,2);
allRTnew(1:10,:)=allRT;
allRT = allRTnew; 

for subjNum = 11:12
    
    subjname = ['at',num2str(subjNum)];
    suffix = '_APES_eye.txt';
%     suffix = '_APES_eye_new.txt';
    eyedata = load([subjname,suffix]);
    
    %     eyedatas{subjNum}=eyedata; % need to increase java heap space
    
    %     %% change the time resolution if needed.
    %     times = eyedata(:,14);
    %     if median(diff(times))==0
    %         newtimes = times;
    %         % poor time resolution
    %         trials = lap(:,3);
    %         [trialStartIdx, trialEndIdx] = separateTrials(trials,3);
    %         for tr = 1:length(trialStartIdx)
    %             samples = trialStartIdx(tr):trialEndIdx(tr);
    %             t = times(samples);
    %             newt = linspace(t(1),t(end),length(t));
    %             newtimes(samples)=newt;
    %         end
    %         eyedata(:,14)=newtimes;
    %     end
    
    
    for cc = 1:2 % 1 for gap, 2 for overlap
        conditionData = eyedata(eyedata(:,2)==cc,:);
        
        trials = conditionData(:,3);
        directions = conditionData(:,4); 
        [trialStartIdx, trialEndIdx] = separateTrials(trials,3);
        
        rts = zeros(1,length(trialStartIdx));
        abnormals = []; % abnormal trials (subject made sacc before target appear)
        tocheck = []; % those trials that gradient is all under 10.
        for tr = 1:length(trialStartIdx)
            samples = trialStartIdx(tr):trialEndIdx(tr);
            trialData = conditionData(samples,:);
            xs = trialData(:,8);
            ts = trialData(:,14);
            sr = 1/((ts(end)-ts(1))/length(ts));
            
            s0 = find(trialData(:,6)==0);
            s1 = find(trialData(:,6)==1);
            s2 = find(trialData(:,6)==2);
            s3 = find(trialData(:,6)==3);
            
            %% check where blinks occur
            blinks = getBlinks2(xs,3000,sr,0);
            
            if length(unique(blinks))>1 && visual == 2
                figure;
                plot(blinks); ylim([-1 2]); hold on
                height = ylim;
                lines = height(1)+[0.2;0.4;0.6;0.8]*(height(2)-height(1));
                % mark different stages
                plot(s0,repmat(lines(1),1,length(s0)),'y','LineWidth',2); hold on;
                plot(s1,repmat(lines(2),1,length(s1)),'g','LineWidth',2); hold on;
                plot(s2,repmat(lines(3),1,length(s2)),'k','LineWidth',2); hold on;
                plot(s3,repmat(lines(4),1,length(s3)),'r','LineWidth',2); hold on;
                title('where do the blinks occur?')
                pause
                close all
            end
            
            %% visualize the trial
            if visual > 1;
                xs(blinks==1)=0;
                figure; plot(xs); hold on
                % mark different stages
                height = ylim;
                lines = height(1)+[0.2;0.4;0.6;0.8]*(height(2)-height(1));
                plot(s0,repmat(lines(1),1,length(s0)),'y','LineWidth',2); hold on;
                plot(s1,repmat(lines(2),1,length(s1)),'g','LineWidth',2); hold on;
                plot(s2,repmat(lines(3),1,length(s2)),'k','LineWidth',2); hold on;
                title('x position of the whole trial','FontSize',14)
                pause
                close all
            end
            
            %% detect saccade
            % all index relative to (within) the trial.
            
            afterFixSamples = [s2; s3];
            afterFix = trialData(afterFixSamples,:);
            x = afterFix(:,8);
            t = afterFix(:,14);
            newt = linspace(t(1),t(end),length(t))';
            v = computeVelEngbert(x,newt);
            g = gradient(x);
            saccStart = find(abs(g)>10,1);
            if isempty(find(abs(g)>10,1))==1
                beep
                tocheck = [tocheck; tr];
                saccStart = find(abs(g)==max(abs(g)),1);
            end
            
            % reaction time:
            targetOnset = length(s2)+1;
            rt = round((t(saccStart)-t(targetOnset))*1000); % ms
            
            if saccStart <= targetOnset
                % abnormal trial
                abnormals = [abnormals; tr];
                rts(tr) = NaN;
            else
                rts(tr) = rt;
            end
            
            
            
            %% visualise saccade detection
            if visual > 0
                
                gapIdx = 1:length(s2);
                
                % ignore blinks
                blink = blinks(afterFixSamples);
                x(blink==1)=0;
                figure('Position',get(0,'ScreenSize'));
                
                subplot(3,1,1); plot(x,'r'); hold on;
                height = ylim;
                width = xlim;
                gapLine = height(1)+0.5*(height(2)-height(1));
                plot(gapIdx,repmat(gapLine,1,length(s2)),'k','LineWidth',2); hold on;
                drawLine = linspace(round(height(1)),round(height(2)),10);
                plot(repmat(saccStart,1,length(drawLine)),drawLine,':','LineWidth',2); hold on;
                note = ['reactionTime = ', num2str(rt),' ms'];
                text(width(2)/2,gapLine,note,'FontSize',14)
                title('x position','FontSize',16)
                
                subplot(3,1,2); plot(g,'r'); hold on
                height = ylim;
                gapLine = height(1)+0.5*(height(2)-height(1));
                plot(gapIdx,repmat(gapLine,1,length(s2)),'k','LineWidth',2); hold on;
                drawLine = linspace(round(height(1)),round(height(2)),10);
                plot(repmat(saccStart,1,length(drawLine)),drawLine,':','LineWidth',2); hold on;
                title('gradient(x)')
                
                subplot(3,1,3); plot(v,'r'); hold on;
                height = ylim;
                gapLine = height(1)+0.5*(height(2)-height(1));
                plot(gapIdx,repmat(gapLine,1,length(s2)),'k','LineWidth',2); hold on;
                drawLine = linspace(round(height(1)),round(height(2)),10);
                plot(repmat(saccStart,1,length(drawLine)),drawLine,':','LineWidth',2); hold on;
                title('velocity')
                
                pause
                close all
            end
        end
%         rts = rts(88:end); % remove the first 87 trials for subj2 overlap
        allRT{subjNum,cc} = rts; 
    end
end

save('allRT','allRT')


%% analysis 
cd '/Users/yingdiliu/Drive/express_saccades'
clear all; 
load allRT

%% 1. average latency for each subj in two conditions 

for subjnum = 1:size(allRT,1)
    
    for cc = 1:2
        
        rt = allRT{subjnum,cc};
        
        allmeans(subjnum,cc)=nanmean(rt); 
        allmedians(subjnum,cc)=nanmedian(rt);
        allstds(subjnum,cc)=nanstd(rt);
        
    end

end


figure; 
plot(allmeans(:,2))
hold on
plot(allmedians(:,2),'--o')
plot(allmeans(:,1),'r')
plot(allmedians(:,1),'r--o')
legend('ovlpMean','ovlpMedian','gapMean','gapMedian')



%% Visualise distribution 

toolboxDir = '/Applications/MATLAB_R2012b.app/toolbox/distrib'; 
addpath(toolboxDir)



figure('Position',get(0,'ScreenSize')); 
for subjNum = 1:size(allRT,1)/2
    
    for cc = 1:2
        
        subplot(size(allRT,1)/2,2,(subjNum-1)*2+cc)
        
        rt = allRT{subjNum,cc}; 
        rt = rt(rt<500); 
        rt = rt(rt>50); 
        
        [n,xout] = hist(rt,50); 
        step = median(diff(xout));
        bar(xout,(n./(length(rt)*step))); 
        xlim([0 500])
        ylim([0 0.04])
        hold on
        
        line120 = repmat(120,100,1);
        plot(line120,linspace(0,0.04,100)','r')
        
        title(['subject ',num2str(subjNum)])
        
        es = length(find(rt<120)); 
        esPct = es/length(rt)*100;
        message1 = ['mean: ', num2str(round(mean(rt))),'ms'];
        message2 = ['ES: ', num2str(round(esPct)),'%'];
        text(350,0.03,message1,'FontSize',15)
        text(350,0.01,message2,'FontSize',15)
    end
end

figure('Position',get(0,'ScreenSize')); 
for subjNum = 7:12
    
    for cc = 1:2
        
        subplot(6,2,(subjNum-7)*2+cc)
        
        rt = allRT{subjNum,cc}; 
        rt = rt(rt<500); 
        rt = rt(rt>50); 
        
        [n,xout] = hist(rt,50); 
        step = median(diff(xout));
        bar(xout,(n./(length(rt)*step))); 
        xlim([0 500])
        ylim([0 0.04])
        hold on
        
        line120 = repmat(120,100,1);
        plot(line120,linspace(0,0.04,100)','r')
        
        title(['subject ',num2str(subjNum)])
        
        es = length(find(rt<120)); 
        esPct = es/length(rt)*100;
        message1 = ['mean: ', num2str(round(mean(rt))),'ms'];
        message2 = ['ES: ', num2str(round(esPct)),'%'];
        text(350,0.03,message1,'FontSize',15)
        text(350,0.01,message2,'FontSize',15)
    end
end



%% fitting ex-gaussian function to saccade latency distribution

% pd = fitdist(rt','Kernel','Kernel','epanechnikov','Width',10);
% x_values = 50:1:500;
% thepdf = pdf(pd,x_values);
% plot(x_values,thepdf)
% 
% param=simple_egfit(rt);
% model = exgausspdf(param(1),param(2),param(3),0:500);
% plot(model,'LineWidth',3)
% 
% area under curve
% auc = trapz(50:500,model);


