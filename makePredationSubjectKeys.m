function keys = makePredationSubjectKeys(all_id_numbers, all_bg_file)
% makePredationSubjectKeys  Build stable subject keys for paired stats.
%
% keys = makePredationSubjectKeys(all_id_numbers, all_bg_file)
%
% Inputs:
%   all_id_numbers : numeric vector (n x 1) prey/object id per sample
%   all_bg_file    : cellstr/string (n x 1) background filename per sample
%
% Output:
%   keys           : string (n x 1) unique-ish key "id|bgfile"
%
% Rationale:
%   In this project, the same sample across conditions is typically defined
%   by the prey/object id and the background/movie identity. Using a key
%   allows robust matching even if file order differs across conditions.

if nargin < 2
    error('Requires all_id_numbers and all_bg_file.');
end

if isempty(all_id_numbers)
    keys = strings(0,1);
    return;
end

ids = all_id_numbers(:);

if iscell(all_bg_file)
    bg = string(all_bg_file(:));
elseif isstring(all_bg_file)
    bg = all_bg_file(:);
elseif ischar(all_bg_file)
    bg = string({all_bg_file});
    bg = bg(:);
else
    error('Unsupported all_bg_file type: %s', class(all_bg_file));
end

n = min(numel(ids), numel(bg));
ids = ids(1:n);
bg = bg(1:n);

keys = string(ids) + "|" + bg;

% Normalize missing
keys(ismissing(keys)) = "";

end
