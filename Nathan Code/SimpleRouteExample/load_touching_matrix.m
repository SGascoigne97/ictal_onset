function [names, touching, roiIndex] = load_touching_matrix(scale)
if scale == 36
    roiIndex = 1;
    load('touching.mat', 'scale36');
    lhnames = scale36.lhnames;
    rhnames = scale36.rhnames;
    lh = scale36.lh;
    rh = scale36.rh;
elseif scale == 60
    roiIndex = 2;
    load('touching.mat', 'scale60');
    lhnames = scale60.lhnames;
    rhnames = scale60.rhnames;
    lh = scale60.lh;
    rh = scale60.rh;
elseif scale == 125
    roiIndex = 3;
    load('touching.mat', 'scale125');
    lhnames = scale125.lhnames;
    rhnames = scale125.rhnames;
    lh = scale125.lh;
    rh = scale125.rh;
elseif scale == 250
    roiIndex = 4;
    load('touching.mat', 'scale250');
    lhnames = scale250.lhnames;
    rhnames = scale250.rhnames;
    lh = scale250.lh;
    rh = scale250.rh;
end
lhmask = logical(~strcmp(lhnames,"unknown").*~strcmp(lhnames,"corpuscallosum"));
rhmask = logical(~strcmp(rhnames,"unknown").*~strcmp(rhnames,"corpuscallosum"));
lhnames = strcat('l.', lhnames);
rhnames = strcat('r.', rhnames);
lhnames = lhnames(lhmask);
rhnames = rhnames(rhmask);

names = cat(1, lhnames,rhnames);
touching = zeros(length(names));
touching(1:length(lhnames),1:length(lhnames)) = lh(lhmask,lhmask);
touching(length(lhnames)+1:end,length(lhnames)+1:end) = rh(rhmask,rhmask);
end
