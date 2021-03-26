function result = invdepfun(allFunctions, lookFor)
% Return all functions that depend on a given function
%
% Example: invdepfun({'myfun1', 'myfun2', 'myfun3'}, 'myfun4') returns all of
% 'myfun1', 'myfun2', 'myfun3' that use 'myfun4'.

    filename = which(lookFor);

    result = {};
    for i = 1 : numel(allFunctions)
        deps = depfun(allFunctions{i}, '-quiet');
        if any(strcmpi(deps, filename))
            result{end + 1} = allFunctions{i};
        end
    end
end