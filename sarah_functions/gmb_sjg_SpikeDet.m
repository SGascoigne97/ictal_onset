function Spk = gmb_sjg_SpikeDet(EEG, cand_sp, Fs, t0, Tw)

% NEED TO ADD PARAMATER VALIDATION STEPS

%% Spike detection Parameters [https://doi.org/10.1073/pnas.1911461117]
% First Step
% Bandpass filtering
% 42Hz in this case to remove 50Hz line noise
% Similar ratio as 50Hz to remove 60 Hz line noise
bp1_f = [20,42];
% Threshold as n-times the STD
tr_n = 3;

% Second Step
% Bandpass filtering
bp2_f = [1,35];

% Median rectifing parameter
scl_Rect = 70; % Diving by the median & multipling by this moves the median to "70uV"

% Morphology Criteria
Mrp_crt = [600, 7, 10];

% Resampling Frequency
rFs = 512;

% Resample to 256Hz
rEEG = resample(EEG,rFs,Fs);

% Number of channels
[nWn,nCh] = size(rEEG);

% Valey search window to each side
vsW = round(0.05*rFs);

% Deal with NaN gaps
% Store NaN positions
% Place there random values
MN = nanmedian(rEEG); MN(isnan(MN)) = 0.5;
ST = nanstd(rEEG); ST(isnan(ST)) = 0.5;
nan_dx = isnan(rEEG);

% Random vlues betweeen med +- std
RND = rand(size(rEEG)).*(ST*2)-ST+MN;
RND(~nan_dx) = 1;
riEEG = rEEG;
riEEG(nan_dx) = 1;

% Combination of both
riEEG = riEEG.*RND;

% Tails @ end and beginng to avoid loss of data
Tl = Tw*rFs;%


% % Spikde detection - This is done separately then added to this algorithm
% as an argument
% % step 1 ----------------------------------------------------------
% % Filtering
% [d_a,d_b] = butter(2,bp1_f/(rFs/2),"bandpass");
% fEEG_1=filtfilt(d_a,d_b,riEEG);
% % fEEG_1 = bandpass(riEEG,bp1_f,rFs);
% fEEG_1(nan_dx) = nan;
% 
% % Thresholding
% cand_sp = zeros(size(fEEG_1));
% for k=1:nCh
%     cand_sp(:,k) = fEEG_1(:,k)>nanstd(fEEG_1(Tl:(nWn-Tl),k))*tr_n;
% end

% step 2 ----------------------------------------------------------
% Filtering

[d_a,d_b] = butter(2,bp2_f/(rFs/2),"bandpass");
fEEG_2=filtfilt(d_a,d_b,riEEG);

%fEEG_2 = bandpass(riEEG,bp2_f,rFs);

fEEG_2(nan_dx) = nan;

% Scaling across all channels; median to 70 uV
aux = fEEG_2(Tl:(nWn-Tl),:);
mEEG = fEEG_2/nanmedian(abs(aux(:)))*scl_Rect;

Spk = cell(1,nCh);

for k=1:nCh
    % Find peaks in this channle
    s = mEEG(:,k);

    if isempty(find(~isnan(s),1))
        continue;
    end

    [~,lck] = findpeaks(s,'MinPeakDistance',vsW);

    % Generate a template and match it to the Thresholding
    temp = zeros(size(s));
    temp(lck) = 1;

    temp = repmat(temp, size(cand_sp(:,k),1)/size(temp,1),1);

    % All peaks
    Locs = find(temp & cand_sp(:,k));
    % Remove the ones on the filler segments
    % Only on segmets 2 (30 to 60s) & 3 (60 t0 90s)
    Locs(Locs<Tl | Locs>(nWn-Tl)) = [];

    if isempty(Locs)
        continue;
    end

    for l=1:length(Locs)
        % Locating Valeys
        % 10 samples window
        ldx = (-vsW:vsW)+Locs(l);
        % Remove outlier datapoints
        ldx(ldx<1 | ldx>length(s)) = [];

        % Valey ~= -Peak
        [~,aux] = findpeaks(-s(ldx));
        v_lck = [max(aux(aux<(vsW+1))),min(aux(aux>(vsW+1)))]-(vsW+1);

        if length(v_lck)<=1
            % What if Just 1 valey?
            % Do not considere the event
            Locs(l) = 0;
            continue;
            
        end

        % step 3 --------------------------------------------------
        % Morphology Criteria

        % Criteria of Amplitude
        PK = s(Locs(l));
        VL = s(v_lck+Locs(l));
        Amp = PK-VL;

        % Criteria of duration
        Dur = abs(v_lck')/rFs*1000;

        % Criteria of slope
        Slope = Amp./Dur;

        % Determine is Spike
        Crit = [Amp,Slope,Dur]>Mrp_crt;

        % All criteria on both sides matched
        Locs(l) = Locs(l)*double(sum(Crit(:))==6);
    end

    % Remove not Spikes
    Locs(~Locs) = [];
    if isempty(Locs)
        Spk{k} = [];
        continue;
    end

    % Storing locations
    Spk{k} = t0+(Locs);%/rFs;%/(60*60*24); % Spk gives time points since the start of recording

end




