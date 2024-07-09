function save_plot(f,dir,name,formats)
% save_plot(f,dir,name,formats)
% Saves the plot figure f in the desired directory (must be subfolder in 
% current working direcory), under the specified name and file format. If 
% the directory does not exist, the function creates the directory.
%
% If name is a scalar, it is turned into a three-character string (e.g., 2
% --> '002', 10 --> '010')
%
% See saveas for possible file formats. To specify multiple file formats,
% use a cell array (e.g., formats = {'jpeg','epsc'}
%
% Gabrielle M. Schroeder
% Newcastle University School of Computing 
% 31 January 2018

% make directory
[status,~,~]=mkdir(dir);
if ~status
    error('Directory could not be created')
end

% make sure name is a string
if ~ischar(name)
    name=num2str(name,'%03.f');
end

% save as each of the requested formats
if iscell(formats)
    n_formats=length(formats);
    for i=1:n_formats
        saveas(f,[dir '/' name],formats{i});
    end
else
    saveas(f,[dir  '/'  name],formats);
end

    