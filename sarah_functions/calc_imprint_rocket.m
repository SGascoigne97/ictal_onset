function [data_tbl, metadata_tbl, cell_imprint,  sz_count_pat] = calc_imprint(data_tbl, metadata_tbl, opts)
% Compute seizure imprint based on EEG recordings

% input:
%   - data_tbl: full data table
%   - metadata_tbl: full metadata table
%   - optional inputs
%       - window_size: window size for which imprint will be computed
%       - min_sz_count: minimum number of seizures to have been recorded per patient
%       - folder: folder to store markers in

% output
%   - data_tbl: full data table (seizures with no imprint have been
%   removed)
%   - metadata_tbl:
%   - cell_imprint
%   - sz_count_pat: count of focal seizures for all patients meeting
%   inclusion criteria for the minimum number of seizures recorded

    arguments
        data_tbl
        metadata_tbl
        opts.window_size (1,1) double {mustBeNumeric} = 1; % Decide window size (in seconds) for which markers and imprint are computed
        opts.min_sz_count (1,1) double {mustBeNumeric} = 5; % Same criteria as earlier but some patients will have fewer seizures following seizures with nio activity in imprint
        opts.folder = 'onset_calcs'; % folder to store markers in
        opts.window_overlap (1,1) double {mustBeNumeric} = 0;
        opts.rec_type (1,1) {mustBeMember(opts.rec_type, ["sec", "prop"])} = "prop" % do we require a consistent window across seizures (sec) or a threshold that varies with seizure durations (prop)
        opts.rec_thresh (1,1) double = 0.1 % Proportion of seizure used to determine location of activity in imprint
        opts.mad_thresh (1,1) double {mustBeNumeric, mustBePositive} = 5
    end
    
    %fill in optional arguments
    window_size = opts.window_size;
    min_sz_count = opts.min_sz_count;
    folder = opts.folder;
    window_overlap = opts.window_overlap;
    rec_type = opts.rec_type;
    rec_thresh = opts.rec_thresh;
    mad_thresh = opts.mad_thresh;
    
    % Set basefolder to store markers
    patient = data_tbl.patient_id{1};
    basefolder = sprintf('%s/%s', folder, patient);
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
    [imprint_out,cell_imprint,~,~] = ms_imprint(metadata_tbl,val_tbl,...
        calcs_ll.t_wndw, "rec_type", rec_type, "rec_thresh",rec_thresh, 'mad_thresh', mad_thresh);  % using the same t_wndw for all features as using same window length or overlap
    
    % Count the number of focal seizures per patient 
    sz_count_pat = tabulate(data_tbl.patient_id);
    sz_count_pat = sz_count_pat(cell2mat(sz_count_pat(:,2))>=min_sz_count,:);

    % Add in segment IDs 
    incl_sz = data_tbl.segment_id;
    cell_imprint = table(cell_imprint, incl_sz, 'VariableNames', {'Imprint'; 'segment_id'});
end
