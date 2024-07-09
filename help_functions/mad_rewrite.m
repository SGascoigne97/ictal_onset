function y = mad_rewrite(x,flag,dim) %uses nanmean/nan median which is causing an error - rewriting here
%MATLAB Code Generation Library Function

%   Limitations:
%   Does not permit empty dim input.

%   Copyright 1993-2018 The MathWorks, Inc.
%#codegen

DOMEDIAN = logical(flag);
if isempty(x)
    y = sum(x,dim,'native','omitnan'); % Should get size and type right.
    y(:) = coder.internal.nan(class(x));
elseif DOMEDIAN
    % Compute the median of the absolute deviations from the median.
    xm = median(x,dim, 'omitnan');
    xx = abs(bsxfun(@minus,x,xm));
    y = median(xx,dim, 'omitnan');
else
    % Compute the mean of the absolute deviations from the mean.
    xm = mean(x,dim, 'omitnan');
    xx = abs(bsxfun(@minus,x,xm));
    y = mean(xx,dim, 'omitnan');
end