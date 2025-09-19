function [reshaped, seqLen] = reshapeAllPaths(all_paths)
    % Validate input
    if ~isnumeric(all_paths)
        error('all_paths must be numeric');
    end

    if isempty(all_paths)
        reshaped = reshape(all_paths, [0,0,2]);
        seqLen = 0;
        return;
    end

    [numTrials, ncols] = size(all_paths);

    if mod(ncols,2) ~= 0
        error('Number of columns in all_paths must be even (interleaved x,y pairs)');
    end

    seqLen = ncols / 2;

    % Preallocate output: (numTrials, seqLen, 2)
    reshaped = zeros(numTrials, seqLen, 2, class(all_paths));

    % Deinterleave: columns 1:2:end -> x, 2:2:end -> y
    x_cols = all_paths(:, 1:2:end);
    y_cols = all_paths(:, 2:2:end);

    % Place into 3D array
    reshaped(:,:,1) = x_cols;
    reshaped(:,:,2) = y_cols;
end