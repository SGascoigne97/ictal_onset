function indexes = find_preserved(names_1,names_2)
[tf,loc] = ismember(names_1,names_2);
[~,p] = sort(loc(tf));
indexes = find(tf);
indexes = indexes(p);
end