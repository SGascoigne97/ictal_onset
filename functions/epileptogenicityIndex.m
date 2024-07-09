function [EI,Nd,t,ERmaster,Na] = epileptogenicityIndex(timeSeries)


% Epileptogenicity of brain structures in human temporal lobe epilepsy: a
% quantified study from intracerebral EEG
% Fabrice Bartolomei, Patrick Chauvel, Fabrice Wendling
% Brain, Volume 131, Issue 7, July 2008, Pages 1818-1830 

%%%%%%%%%
% Usage 
%%%%%%%%%

% Input: 

% timeSeries 
% This is a n x m matrix, containing n samples recorded from m channels. 

% Output: 

% EI: epileptigenicity index. 
% This is an n x 1 matrix of epileptogenicity index values. 
% [~, ind] = max(EI). ind here corresponds to the single lead with maximal
% EI. 

% Nd 
% This is the detection time. This is the time at which each lead is
% "flagged" as first having a rise in high-frequency activity. 
% Note that Nd is not an absolute time, but rather serves as an index into
% t. 
% NdAbsolute = t(Nd) 

% t 
% This is the timescale used by the output variables. 

% ERmaster
% This is an n x length(t) matrix. 
% It corresponds to the energy ratio of each lead over time (where time is
% given by t). 
% The energy ratio is higher when there is greater high-frequency activity.
% It's given by the ratio of beta and gamma frequency power over theta and
% alpha frequency power. 
% The EI is calculated by finding leads that have early and large rises in
% ER. 

% Na 
% Alarm time. 
% Like Nd, this is an n x 1 matrix. 
% Somewhat similar to detection time, but occurs later. 
% Detection time is when high-frequency activity begins. 
% Alarm time on the other hand is when high-frequency activity becomes
% great enough to rise above threshold and "alarm". 
% Note that Na is not an absolute time, but rather serves as an index into
% t. 
% NaAbsolute = t(Na) 


%% Preamble 

freqWindows = [3.5 7.4; 7.4 12.4; 12.4 24; 24 140];
% theta; alpha; beta; gamma 


Fs = 256;            % Sampling frequency

%% Definition of a statistic that abruptly increases as fast oscillations appear in the signal  

[~,f,t] = spectrogram(timeSeries(:,1),1000,800,[],Fs);
ERmaster = zeros(size(timeSeries,2),length(t));

[~,inds] = sort(abs(f - 60),'ascend');
inds = inds(1:20); % Remove line noise 
f(inds) = []; 


for jj = 1:size(timeSeries,2)
    
    X = timeSeries(:,jj);
    
    s = spectrogram(X,1000,800,[],Fs);
    s(inds,:) = [];
    
    G = s .* conj(s) / (2 * pi);
    ER = zeros(size(freqWindows,1),length(t));
    
    for ii = 1:size(freqWindows,1)
        
        [~, range1] = min(abs(f - freqWindows(ii,1)));
        [~, range2] = min(abs(f - freqWindows(ii,2)));
        
        ER(ii,:) = sum(G(range1:range2,:));
    end
    ER = (ER(3,:) + ER(4,:)) ./ (ER(1,:) + ER(2,:));
    
    ERmaster(jj,:) = ER;
    
end

%% Optimal detection of rapid discharges  

ERn = cumsum(ERmaster,2) ./ (1:size(ERmaster,2));

v = .5;
% This number sort of serves as a built-in "time cost". High-frequency
% activity needs to be robust enough to surpass this v number and still
% overcome the threshold. 
% It's arbitrary but I think 0.5 is what Bartolomei used. 

v = ones(size(ERmaster)) * -v;

Un = ERmaster - ERn + v; 
Un = cumsum(Un,2);

% UnTemp = Un; 

lambda = 15; 
% Lambda also serves as a somewhat arbitrary threshold value. 

Nd = nan(size(Un,1),1);
Na = ones(size(Un,1),1) * length(t);

for ii = 1:size(Un,2)
    [un,ind] = min(Un(:,1:ii),[],2);
    
    UnDiff = Un(:,ii) - un; 
    pastThreshold = UnDiff > lambda; 
    
    pastThreshold(~isnan(Nd)) = false; 
    
    Nd(pastThreshold) = ind(pastThreshold);
    Na(pastThreshold) = ii;
    
%     Un(pastThreshold,ii:end) = Un(pastThreshold,ii:end) - Un(pastThreshold,ii);
%         
%     UnTemp(pastThreshold,:) = Un(pastThreshold,:);
%     UnTemp(pastThreshold,1:ii - 1) = nan; 
   
end



%% Definition and computation of the EI 

tau = 1; 
H = 5; 

hSpan = find(t >= H,1);

Nd(isnan(Nd)) = length(t) - hSpan;
Nd(Nd > length(t) - hSpan) = length(t) - hSpan;

EI = nan(size(Nd));

N0 = min(Nd); 


for ii = 1:length(Nd) 
    EI(ii) = 1 / (t(Nd(ii)) - t(N0) + tau) * mean(ERmaster(ii,Nd(ii):Nd(ii) + hSpan));
end

EI = EI / max(EI(:));


%% Plotting

return 
[~,inds] = sort(EI,'descend'); 
for ii = 1:5
    
    % if isnan(Nd(inds(ii))); continue; end
    clf
    plot(t,Un(inds(ii),:))
    hold on
    plot(t(Nd(inds(ii))),Un(inds(ii),Nd(inds(ii))),'ro')
    pause
end


end
    

