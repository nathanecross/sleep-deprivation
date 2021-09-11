function name = find_preproc_file(path, keyword1, keyword2)
    
    files = struct2table(dir(char(path)));
    if nargin <3
        index = find(contains(files.name, keyword1) & startsWith(files.name,'.')==0);
    else
        index = find(contains(files.name, keyword1) & contains(files.name, keyword2) & startsWith(files.name,'.')==0);
    end
    name = files.name(index);
    
end