function json_data = db_extract_VI_status(json_data)
% Short function to extract VI status for seizure data so we can remove
% non-visually confirmed seizures befroe running the pipeline
%
% Sarah J. Gascoigne
% CNNP Lab, Newcastle University
% March 2023

arguments
    json_data struct
end
    for seg = 1:length(json_data)
        sz_details = json_data(seg).seizure_details;
        json_data(seg).VI_status = string(sz_details.sz_VI_status);
    end

end