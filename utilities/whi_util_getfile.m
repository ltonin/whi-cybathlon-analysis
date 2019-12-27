function files = whi_util_getfile(root, type, varargin)
% WHI_UTIL_GETFILE Get file with inclusion/exclusion criteria
%
%   files = WHI_UTIL_GETFILE(root, type) returns all files in the directory 
%   root and in the subdirectories with the extension that match with 'type'.
%
%   files = WHI_UTIL_GETFILE(root, type, 'include', {PATTERN1, PATTERN2, ...})
%   returns all files that match the provided patterns in the filename
%
%   files = WHI_UTIL_GETFILE(root, type, 'exclude', {PATTERN1, PATTERN2, ...})
%   returns all files excluding the one matching the provided patterns in
%   the filename



    %% Input parser
    
    % Default patterns
    defaultInclude = {};
    defaultExclude = {};
    defaultLevel = 1;
    
    % Validation functions
    isextension = @(x) assert(ischar(x) && isequal(regexp(x, '\.\w+'), 1), 'type must be an extension (e.g., ''.gdf'')') ;
    isrootdir   = @(x) assert(isfolder(x), 'root must be a valid directory');
    ispattern   = @(x) assert(iscellstr(x) || ischar(x), 'inclusion/exclusion patterns must be a char array or a cell of chars');
    isnumber    = @(x) assert(isnumeric(x) && isscalar(x), 'level must be a scalar');
    
    % Add parameters
    p = inputParser;
    p.addRequired('root', isrootdir);
    p.addRequired('type', isextension);
    p.addParameter('include', defaultInclude, ispattern);
    p.addParameter('exclude', defaultExclude, ispattern);
    p.addParameter('level', defaultLevel, isnumber);
    
    % Parse input
    parse(p, root, type, varargin{:});
    
    % Convert patterns to cell array
    rule_include    = char2cell(p.Results.include);
    rule_exclude    = char2cell(p.Results.exclude);
    level           = p.Results.level;
    
    %% Get files (recurserverly)
    list = get_files(root, type, level);
    
    %% Include/Exclude patterns
    idx_include = include(list, rule_include);
    idx_exclude = exclude(list, rule_exclude);
    
    idx = idx_include & ~idx_exclude;
  
    %% Select valid files
    files = list(idx);
    nfiles = length(files);
    
    if(nfiles == 0)
        error(['[io] - No files found with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);
    end
end

function files = get_files(root, type, level)

    % list entries in the given directory directory
    entries = dir(root);
    
    % remove '.' and '..' entries
    entries = entries(~ismember({entries.name}, {'.', '..'}));
    
    nentries = length(entries);
    
    files = {};
    
    for eId = 1:nentries
        centry = entries(eId);
        
        % If it is a dir, go recurseverly, otherwise concat filename 
        if (centry.isdir == true) && level > 1
            files = cat(1, files, get_files(fullfile(centry.folder, centry.name), type, level-1));
        else
            if matchextension(centry.name, type)
                files = cat(1, files, fullfile(centry.folder, centry.name));
            end
        end
        
    end
end

function idx = include(list, patterns)

    npatts = length(patterns);
    nfiles = length(list);
    idx  = true(nfiles, 1);
    
    for pId = 1:npatts
        idx = idx & contains(list, patterns{pId});
    end
    
end

function idx = exclude(list, patterns)
    npatts = length(patterns);
    nfiles = length(list);
    idx  = false(nfiles, 1);
    
    for pId = 1:npatts
        idx = idx | contains(list, patterns{pId});
    end
end


function matched = matchextension(filename, extension)
    matched = false;
    [~, ~, fileext] = fileparts(filename);
    if(strcmp(fileext, extension) == true)
        matched = true;
    end
end

function out = char2cell(in)
    
    out = in;
    if ischar(in)
        out{1} = in;
    end
end