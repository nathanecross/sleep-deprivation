function Tasktimes = concatenate_taskonsets(file, times, keyword1)
       
    load('times/Tasktimes.mat');
    fields = fieldnames(Tasktimes);
    fields = fields(contains(fieldnames(Tasktimes),'_conv'));
    
    for t = 1:length(times)
        for f = 1:length(file)
            filename = char(file(f));
            task = extractAfter(filename,'-');
            index = fields(logical(contains(fields,task).*contains(fields,times(t,3))));
            if f == 1
                Tasktimes.(sprintf('All_%s_conv',string(times(t,3)))) = Tasktimes.(sprintf('%s',string(index)));
            elseif f > 1
                Tasktimes.(sprintf('All_%s_conv',string(times(t,3)))) = [Tasktimes.(sprintf('All_%s_conv',string(times(t,3)))) Tasktimes.(sprintf('%s',string(index)))];
            end
        end
    end

end