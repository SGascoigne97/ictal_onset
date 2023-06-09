function [comp_table] = comp_auto_clo(clo, auto, opts)
% Compare automatically detected onset with clincially labelled onset on an
% individual basis
% Sarah J Gascoigne 05/06/2023

% input:
%   - clo:  binary array of clinically labelled onset
%   - auto: binary array or matrix (same number of rows as clo) showing
%           automatically detected onset for all seizures
%   
%   - optional inputs
%       - per_sz_or_all_sz: string indicating if output should be per seizure
%           or with respect to onset across seizures
%       - min_sz: if using all_sz, this gives the threshold for the number
%           of seizures needed for the region to be included in comaprison
%       - n_perm: number of permutations used for computing z-score

% output
%   - comp_table: table containing raw and normalised Jaccard's index and 
%     percentage of CLO captured

    arguments
        clo (:,1) logical % binary array of clinically labelled onset
        auto logical {mustHaveSameNrows(clo, auto)} % onset_output table for one specific patient 
        
        opts.per_sz_or_all_sz (1,1) string {mustBeMember(opts.per_sz_or_all_sz, ['all_sz'; 'per_sz'])} = 'all_sz'
        opts.sz_prop_thresh (1,1) {mustBeNumeric(opts.sz_prop_thresh),...
            mustBeLessThan(opts.sz_prop_thresh, 1),...
            mustBeGreaterThan(opts.sz_prop_thresh, 0)} = 0.25
        opts.n_perm (1,1) {mustBeNumeric(opts.n_perm)} = 1000 
       
    end
    
    %fill in optional arguments
    per_sz_or_all_sz = opts.per_sz_or_all_sz;
    if matches(per_sz_or_all_sz, 'all_sz')
        sz_prop_thresh = opts.sz_prop_thresh;
        min_sz = ceil(sz_prop_thresh*size(auto,2));
        auto = summarise_across_szs(auto, min_sz);
    end
    n_perm = opts.n_perm;

    % Compute compatison measues (Percentage of CLO captured and Jaccard's
    % index between CLO and automatically detected onset)
    % In both cases we will also use a permutation test to determine if the
    % value is greater or less than expected based on chance

    if sum(clo) == 0
        comp_table = [];
        return
    end
    
    for sz = 1:size(auto,2)
        auto_sz = auto(:,sz);
        % Compute Jaccard's index
        jacc = jaccard(clo,auto_sz);
        % Compute the percentage of CLO captured by automatic detection
        % algorithm
        perc_clo_capt = sum(auto_sz + clo == 2)/sum(clo);

        % Perform permutation test 
        jacc_perm = zeros(1,n_perm);
        perc_perm = jacc_perm; 
        for perm = 1:n_perm
            rng(perm)
            perm_clo = clo(randperm(length(clo)));
            rng(perm+1)
            perm_auto = auto_sz(randperm(length(auto_sz)));
            jacc_perm(perm) = jaccard(perm_clo,perm_auto);
            perc_perm(perm) = sum(perm_auto + perm_clo == 2)/sum(perm_clo);
        end

        jacc_z = (jacc - mean(jacc_perm))/std(jacc_perm);
        perc_z = (perc_clo_capt - mean(perc_perm))/std(perc_perm);

        % Compose a table as output
        sz_table = [perc_clo_capt, perc_z, jacc, jacc_z];

        % Add seizure to comparison table
        if size(sz_table,1) > 0
            if exist('comp_table', 'var')
                comp_table = [comp_table; sz_table];
            else
                comp_table = sz_table;
            end
        end
    end

    comp_table = array2table(comp_table);
    comp_table.Properties.VariableNames = {'Percentage of CLO captured'; ...
        'Perc_z'; 'Jaccard'; 'Jaccard_z'};

    for col = 1:4
        col_name = comp_table.Properties.VariableNames(col);
        comp_table.(col_name{:}) = comp_table.(col_name{:});
    end
end

% Custom validation function
function mustHaveSameNrows(a,b)
    % Test for equal size
    if ~isequal(size(a,1),size(b,1))
        eid = 'Size:notEqual';
        msg = 'Number of rows in CLO and AUTO must be the same';
        throwAsCaller(MException(eid,msg))
    end
end


