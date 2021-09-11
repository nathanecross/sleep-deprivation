function [SubjectList, NumberOfSubjects] = cycle_directory_contents(path)

clc
cd(path); 
di = dir; di = struct2table(di);
di = di(:,1); di = table2cell(di);
IDZ = strfind(di(:,1),'.'); TID = cellfun('isempty', IDZ);
SubjectList = di(TID);
subjonly = ~contains(SubjectList,'derivatives');
SubjectList = SubjectList(subjonly);

NumberOfSubjects = length(SubjectList);

end

