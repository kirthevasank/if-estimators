function ifSetup(path)
  if ~exist('path', 'var') | isempty(path)
    fprintf('Path not given. Using current working directory.\n');
    path = pwd;
  end
  addpath(genpath(path));
end

