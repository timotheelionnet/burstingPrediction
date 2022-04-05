function bgData = readBioGridFile(fname)
% from a tab3 format BIOGRID download file (fname is the path to that file)
% this function collects all *physical* interactions between *human proteins*, 
% outputs the result as a 2 column list where the first column is the gene
% symbol of the bait, the second column is the gene symbol of the hit.

% read first line to get the header names
fid = fopen(fname);
h = fgetl(fid);

%% find out the indices of the columns for "Official Symbol Interactor A" and
% "Official Symbol Interactor B" 
hcell = strsplit(h,'\t');

c1 = find(ismember(hcell,'Official Symbol Interactor A')); % bait col
c2 = find(ismember(hcell,'Official Symbol Interactor B')); % interactor col
c3 = find(ismember(hcell,'Experimental System Type')); % physical or genetic interaction
c4 = find(ismember(hcell,'Organism Name Interactor A')); % bait organism
c5 = find(ismember(hcell,'Organism Name Interactor B')); % interactor organism

% it is anticipated that c1<c2<c3<c4

fmt = [ repmat('%*s',1,c1-1), ...
    '%s ', repmat('%*s',1,c2-c1-1),...
    '%s ', repmat('%*s',1,c3-c2-1),...
    '%s ', repmat('%*s',1,c4-c3-1),...
    '%s ', repmat('%*s',1,c5-c4-1),...
    '%s ',repmat('%*s',1,numel(hcell)-c5)];

% read the columns for baits and interactors

bgData = textscan(fid, fmt,'Delimiter','\t');

fclose(fid);

% reformat as one cell array
bgData = [bgData{:,1},bgData{:,2},bgData{:,3},bgData{:,4},bgData{:,5}];

%% clear non-human interactions
nStart = size(bgData,1);
bgData = bgData(ismember(bgData(:,3),'physical') ...
      & ismember(bgData(:,4),'Homo sapiens') ...
      & ismember(bgData(:,5),'Homo sapiens'),1:2);

disp(['Excluded ',num2str(nStart - size(bgData,1)),...
    ' genetic + non-human interactions out of ',num2str(nStart),...
    '. There remain ', num2str(size(bgData,1)),...
    ' human physical interactions.']);

