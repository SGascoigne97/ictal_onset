function [data_tbl, cell_imprint,  sz_count_pat] = ...
    calc_imprint_mahal(data_tbl, opts)% Compute seizure imprint based on EEG recordings

% input:
%   - data_tbl: full data table
%   - optional inputs
%       - window_size: window size for which imprint will be computed
%       - min_sz_count: minimum number of seizures to have been recorded per patient
%       - folder: folder to store markers in

% output
%   - data_tbl: full data table (seizures with no imprint have been
%   removed)
%   - cell_imprint: table including segment_id, imprint (channels x timepoints)
%                   MAD scores, preictal MAD scores, preictal Mahalanobis
%                   distances
%   - sz_count_pat: count of focal seizures for all patients meeting
%   inclusion criteria for the minimum number of seizures recorded

    arguments
        data_tbl
        opts.window_size (1,1) double {mustBeNumeric} = 1; % Decide window size (in seconds) for which markers and imprint are computed
        opts.min_sz_count (1,1) double {mustBeNumeric} = 5; % Same criteria as earlier but some patients will have fewer seizures following seizures with nio activity in imprint
        opts.folder = 'onset_calcs'; % folder to store markers in
        opts.window_overlap (1,1) double {mustBeNumeric} = 0;
        opts.rec_thresh (1,1) double = 0.1 % Proportion of seizure used to determine location of activity in imprint
        opts.mad_thresh (1,1) double {mustBeNumeric, mustBePositive} = 5
        opts.ict_buffer (1,1) double {mustBeNumeric} = 10 % Number of seconds to shift back clinically labelled onset time to ensure 'true' onset is captured
    end
    
    %fill in optional arguments
    window_size = opts.window_size;
    min_sz_count = opts.min_sz_count;
    folder = opts.folder;
    window_overlap = opts.window_overlap;
    rec_thresh = opts.rec_thresh;
    mad_thresh = opts.mad_thresh;
    ict_buffer = opts.ict_buffer;
    
    % Set basefolder to store markers
    patient = data_tbl.patient_id{1};
    basefolder = sprintf('%s/%s', folder, string(patient));
    sampling_rate=data_tbl.segment_fs(1) ;%just using the first, as all sampling was the same in our data after preproc.
    
    linelength_db = LL_db([basefolder '/LL_db/']); %setup folder for all Line Length measures
    linelength_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap);
    linelength_db.paramset_tbl                       % display all currently tracked paramsets
    linelength_db.calc(data_tbl,[]);                       % calculate all parametersets for all segments in data
    
    energy_db = Energy_db([basefolder '/Energy_db']);
    energy_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap);
    energy_db.paramset_tbl                       % display all currently tracked paramsets
    energy_db.calc(data_tbl,[]);                       % calculate all parametersets for all segments in data
    
    %bands: [1 4; 4 8; 8 13; 13 30; 30 60, 60 100]; 
    bandpower_db = BP_db([basefolder '/BP_db']);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[1 4]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[4 8]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[8 13]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[13 30]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[30 60]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[60 100]);
    bandpower_db.paramset_tbl                       % display all currently tracked paramsets
    bandpower_db.calc(data_tbl,[]);   
    
    
    %% Import markers (line length, energy, bandpowers)
    calcs_ll = linelength_db.get(1);    % get all calculation outputs in variable. 
    calcs_energy = energy_db.get(1);
    calcs_bp_delta = bandpower_db.get(1);
    calcs_bp_theta = bandpower_db.get(2);
    calcs_bp_alpha = bandpower_db.get(3);
    calcs_bp_beta = bandpower_db.get(4);
    calcs_bp_gamma = bandpower_db.get(5);
    calcs_bp_hgamma = bandpower_db.get(6);
    
    %% calculate imprint and recruitment markers
    val_tbl=[calcs_ll.LL_ms calcs_energy.energy_ms calcs_bp_delta.bp ...
        calcs_bp_theta.bp calcs_bp_alpha.bp calcs_bp_beta.bp calcs_bp_gamma.bp calcs_bp_hgamma.bp];
    %% Create a table of seizure data
    sz_mat_tab = table(data_tbl.segment_id, cell(size(val_tbl,1),1), calcs_ll.t_wndw,...
        'VariableNames',["segment_id", "feat_mat", "tw"]);
    for sz = 1:size(val_tbl,1)
        sz_mat=log(cat(3,val_tbl{sz,:})); % Here we log-transform markers 
        sz_mat_tab.feat_mat{sz} = sz_mat;
    end
    [imprint_out,cell_imprint,~,cell_madscores,cell_pre_features_mad, cell_pre_mahal_mat] = mahal_imprint(data_tbl,sz_mat_tab,...
        calcs_ll.t_wndw, "rec_thresh", rec_thresh, 'mad_thresh', mad_thresh, "ict_buffer", ict_buffer);  % using the same t_wndw for all features as using same window length or overlap
    
    % Only keep seizures with seizure activity detected (activity in at least
    % one channel in imprint) 
%     incl_sz = data_tbl.segment_id;
    rm_sz_no_onset = false(size(imprint_out.num_chan_imprint==0)); % hack
    
    rm_sz_mad =  cellfun(@max, cellfun(@max,cell_madscores.mahal_MAD,'UniformOutput',false)) < 5;
    rm_sz = (rm_sz_no_onset + rm_sz_mad) > 0;

    %% Uncomment this section if we choose to remove seizures with low maximum MAD
    %rm_sz_no_imprint = false(size(imprint_out.num_chan_imprint==0)); % hack
    % Remove seizures with low maximum MAD score as activity is not
    % sufficiently distinct from preictal activity
    %rm_sz_low_MAD = cellfun(@max,cellfun(@max, cell_madscores.("MAD scores") , 'UniformOutput', false)) <5;
    %rm_sz = (rm_sz_no_imprint + rm_sz_low_MAD) > 0;

    data_tbl = data_tbl(~rm_sz,:);
    cell_imprint = cell_imprint(~rm_sz,:);
    % Count the number of focal seizures per patient 
    sz_count_pat = tabulate(data_tbl.patient_id);
    sz_count_pat = sz_count_pat(cell2mat(sz_count_pat(:,2))>=min_sz_count,:);
    
    % Only include patients meeting inclusion criteria
    data_tbl = data_tbl(matches(string(data_tbl.patient_id), string(sz_count_pat(:,1))),:);
    cell_imprint.cell_madscores = cell_madscores;  
    cell_imprint.cell_pre_features_mad = cell_pre_features_mad;
    cell_imprint.cell_pre_mahal_mat = cell_pre_mahal_mat;
    cell_imprint = cell_imprint(matches(string(cell_imprint.segment_id), string(data_tbl.segment_id)),:);
end
