SNR_th = 10;

topLevelFolder = pwd; % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'
% Get a list of all files and folders in this folder.
files = dir(topLevelFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..


start_date = '20230718';
start_time = '233801_234335';
end_date = '20230723';
end_time = '101022_093059';
for i = 1:length(subFolderNames)
    if subFolderNames{i} == start_date
        start_folder = i;
    end
    if subFolderNames{i} == end_date
        end_folder = i;
    end
end


for j=start_folder:end_folder
topLevelFolder2 = ['C:\Users\zivaqua.WISMAIN\Documents\Research\Q-router analysis\',subFolderNames{j}] ; % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'
% Get a list of all files and folders in this folder.
files2 = dir(topLevelFolder2);
% Get a logical vector that tells which is a directory.
dirFlags2 = [files2.isdir];
% Extract only those that are directories.
subFolders2 = files2(dirFlags2); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames2 = {subFolders2(3:end).name}; % Start at 3 to skip . and ..

if start_folder==end_folder
    xx = find(strcmp(subFolderNames2, start_time));
    yy = find(strcmp(subFolderNames2, end_time));
elseif j==start_folder
    xx = find(strcmp(subFolderNames2, start_time));
    yy = (length(subFolderNames2)-2);
elseif j==end_folder
    xx = 1;
    yy = find(strcmp(subFolderNames2, end_time));
else
    xx = 1;
    yy = (length(subFolderNames2)-2);
end

    for k=xx:yy
        load([topLevelFolder2,'\',subFolderNames2{k},'\workspace.mat'],'SNR','num_transits')
        load([topLevelFolder2,'\',subFolderNames2{k},'\analysis.mat'])
%         load([topLevelFolder2,'\',subFolderNames2{k},'\workspace.mat'])

        if k==xx & j==start_folder
            cR = num_transits*0;
            cT = num_transits*0;
            tBcR = num_transits*0;
            tDcR = num_transits*0;
            tRcR = num_transits*0;
            tTcR = num_transits*0;
            num_transits_tot = num_transits*0;
        end

        num_transits_tot = num_transits_tot + num_transits;
        cR = cR + control_R.*(SNR>SNR_th);
        cT = cT + control_T.*(SNR>SNR_th);
        tBcR = tBcR + B_cR.*(SNR>SNR_th);
        tDcR = tDcR + D_cR.*(SNR>SNR_th);
        tRcR = tRcR + target_R_cR.*(SNR>SNR_th);
        tTcR = tTcR + target_T_cR.*(SNR>SNR_th);
    
    end
end


c_routing_eff = cR./(cR+cT);
c_routing_pts = (cR+cT);
t_routing_fid = tRcR./(tRcR+tTcR);
t_routing_pts = (tRcR+tTcR);
t_coherence = tBcR./(tBcR+tDcR);
t_coherence_pts = (tBcR+tDcR);

cond = (c_routing_eff>0.74).*(t_routing_fid>0.74).*(t_coherence>0.85).*(t_coherence_pts>20);



