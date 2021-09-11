function [matrix_ordered, matrix_bordered, L] = reorder_matrices(matrix, nparc, nnetw)
    NumberOfParcels = str2num(nparc);
    % create header of parcel names based off lookup table
    lut = transpose(table2cell(readtable(sprintf('labels/Schaefer2018_%sParcels_%sNetworks_order.txt', nparc,nnetw))));
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
    lut_header_y = lut_header_x'; %create y-axis header
    lut = cell2table(lut); % make lut a table
    Session_mean_tab = array2table(matrix);
    Session_mean_cell = num2cell(matrix);
    % Create a matrix of session mean with y-axis header of parcel labels:
    net = num2cell(zeros(NumberOfParcels,NumberOfParcels+1));
    for i=1:NumberOfParcels
        a = char(lut_header_x(1,i));
        lut.Properties.VariableNames{i} = a;
        Session_mean_tab.Properties.VariableNames{i} = lut.Properties.VariableNames{i};
        net{i,1} = Session_mean_tab.Properties.VariableNames{i};
        net(i,2:end) = Session_mean_cell(i,:);
    end
    top = net(1:end,1)'; filler = num2cell(zeros(1)); top = [filler top]; % create x-axis header of parcel labels
    net = [top; net]; % Create a matrix of session mean with x- and y-axis headers of parcel labels
    nettab = cell2table(net); % Make a table in order to sort parcels
    nettab.net1{1} = 'A'; % Dummy code the first row
    B = sortrows(nettab, 'net1'); 
    B = table2cell(B);
    B = B';
    B = cell2table(B);
    %B.B2 = nettab.net2;
    %B.B1{2} = 'A';
    B = sortrows(B, 'B1');
    B = table2cell(B);
    matrix_ordered = cell2mat(net(2:end,2:end));
    L = B(:,1);
    % Create a bordered matrix for visualisation
    Session_mean_bordered = B(2:end,2:end);
    Session_mean_bordered = cell2mat(Session_mean_bordered);
    Session_mean_bordered = [zeros(1,NumberOfParcels); Session_mean_bordered];
    Session_mean_bordered = [zeros(NumberOfParcels+1,1) Session_mean_bordered];
    Session_mean_bordered = [Session_mean_bordered; zeros(1,NumberOfParcels+1)];
    Session_mean_bordered = [Session_mean_bordered zeros(NumberOfParcels+2,1)];
    matrix_bordered = Session_mean_bordered;


    
end
