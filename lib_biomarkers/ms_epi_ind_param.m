function [EI_tbl,Nd_tbl,ER_tbl,Onset_chan_tbl] = ms_epi_ind_param(data_tbl,opts)
% calculates epilpetogenicity index. Implementation based on Bartolomei et
% al., 2008
% 
% inputs:
%   - data_tbl: full data table
%   - opts: for more info on parameters look below (HAVEN'T PUT THIS IN
%   YET)
%
% outputs:
%   - EI: epileptigenicity index - This is an n x 1 matrix of 
%         epileptogenicity index values. 
%   - Onset_chan ([~, ind] = max(EI)) - corresponds to the single lead 
%         with maximal EI. 
%   - ERmaster - the energy ratio of each lead over time.
%   - Nd - detection time (in seconds)
% 
% outputs are in order of meta_tbl, as we assume that input is in order of
% meta_tbl

arguments
        data_tbl
%         opts.rec_thresh (1,1) double {mustBeNumeric, mustBePositive} = 0.1
%         opts.ict_buffer (1,1) double {mustBeInteger, mustBePositive} = 10 %in seconds
%         opts.mad_thresh (1,1) double {mustBeNumeric, mustBePositive} = 5
%         opts.hsc       (1,1) double {mustBeNumeric, mustBePositive} = 0.1 %percentage of seizure to be ignored at the start for calculation of max spread
        opts.v (1,1) double {mustBeNumeric, mustBePositive}= 0.5 % -1: no plotting, 0:new figure, k>0: figure with index k
        opts.lambda (1,1) double {mustBeNumeric, mustBePositive} = 15
end

    v=opts.v;
    lambda=opts.lambda;

%% Preamble 

freqWindows = [3.5 7.4; 7.4 12.4; 12.4 24; 24 140];
% freqWindows = [4 8; 8 13; 13 30; 30 100];
% theta; alpha; beta; gamma 

% Define tables to store output
EI_tbl = cell(size(data_tbl,1),1);
Nd_tbl = cell(size(data_tbl,1),1);
ER_tbl = cell(size(data_tbl,1),1);
Onset_chan_tbl = cell(size(data_tbl,1),1);

% Iterate through seizures
for sz = 1:size(data_tbl,1)
    %% Definition of a statistic that abruptly increases as fast oscillations appear in the signal  
    Fs = data_tbl.segment_fs(sz);            % Sampling frequency
    sz_data = data_tbl.segment_data{sz,:};
    timeSeries = sz_data(:,Fs*data_tbl.segment_pre(sz):size(sz_data,2)-(Fs*data_tbl.segment_post(sz)));
    
    [~,f,t] = spectrogram(timeSeries(1,:),Fs,0,[],Fs);
    ERmaster = zeros(size(timeSeries',2),length(t));
    
    %[~,inds] = sort(abs(f - 60),'ascend');
    %inds = inds(1:20); % Remove line noise 
    %f(inds) = []; 
    
    
    for jj = 1:size(timeSeries,1)
        
        X = timeSeries(jj,:);
        
        s = spectrogram(X,Fs,0,[],Fs);
        %s(inds,:) = [];
        
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
    
    %v = .5;
    % This number sort of serves as a built-in "time cost". High-frequency
    % activity needs to be robust enough to surpass this v number and still
    % overcome the threshold. 
    % It's arbitrary but I think 0.5 is what Bartolomei used. 
    
    v = ones(size(ERmaster)) * -v;
    
    Un = ERmaster - ERn + v; 
    Un = cumsum(Un,2);
    
    % UnTemp = Un; 
    
    %lambda = 15; 
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
    %return 
    % [~,inds] = sort(EI,'descend'); 
    % for ii = 1:5
    %     
    %     % if isnan(Nd(inds(ii))); continue; end
    %     clf
    %     plot(t,Un(inds(ii),:))
    %     hold on
    %     plot(t(Nd(inds(ii))),Un(inds(ii),Nd(inds(ii))),'ro')
    %     pause
    % end


    EI_tbl{sz} = EI; 
    Nd_tbl{sz} = Nd;
    ER_tbl{sz} = ERmaster;
    [~, Onset_chan_tbl{sz}] = max(EI);
 end


end
    

