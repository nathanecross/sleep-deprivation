function [vector_ordered, L] = reorder_matrices(vector, nparc, nnetw, lut)
    NumberOfParcels = str2num(nparc);
    lut_header_x = lut(2, :); 
    netremove = sprintf('%sNetworks_',nnetw); istart = length(netremove); %remove first part of string
    % 
    % rename each parcel for easier sorting:
    %(new naming convention: network, hemisphere, parcel#)
    for i = 1:length(lut_header_x)
        lutcell = lut_header_x(1,i);
        lut_header_x(1,i) =  {lutcell{1,1}(istart+1:end)};
        lutcell = lut_header_x(1,i);
        hemi = char({lutcell{1,1}(1:2)});
        netw = char(extractBetween({lutcell{1,1}},"H_","_"));
        lut_header_x(1,i) = cellstr([netw '_' hemi '_' num2str(i)]);
    end
    lut = cell2table(lut); % make lut a table
    Session_mean_tab = array2table(vector');
    Session_mean_cell = num2cell(vector');
    % Create a matrix of session mean with y-axis header of parcel labels:
    net = num2cell(zeros(2,NumberOfParcels));
    for i=1:NumberOfParcels
        a = char(lut_header_x(1,i));
        lut.Properties.VariableNames{i} = a;
        Session_mean_tab.Properties.VariableNames{i} = lut.Properties.VariableNames{i};
        net{1,i} = Session_mean_tab.Properties.VariableNames{i};
        net(2,i) = Session_mean_cell(1,i);
    end
    
    nettab = cell2table(net'); % Make a table in order to sort parcels
    nettab = sortrows(nettab, 'Var1');
    vector_ordered = nettab.Var2;
    L = nettab.Var1;

end