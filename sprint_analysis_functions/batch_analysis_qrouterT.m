SNR_th = 1;


topLevelFolder = pwd; % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'
% Get a list of all files and folders in this folder.
files = dir(topLevelFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:6).name}; % Start at 3 to skip . and ..
% % 
% figure('Name','North')
% hold on
% bar(edges,histogram_N)


% nS
start_date = '20240118';
start_time = '163621_161746';
end_date = '20240118';
end_time = '181945_180203';
% 
% % nN
% start_date = '20240118';
% start_time = '142932_140759';
% end_date = '20240118';
% end_time = '145424_143246';
% %  
% % sN
% start_date = '20240118';
% start_time = '170100_164118';
% end_date = '20240118';
% end_time = '173303_171442'; 
% % % sS
% start_date = '20240118';
% start_time = '163621_161746';
% end_date = '20240118';
% end_time = '163621_161746';


for i = 1:length(subFolderNames)
    if subFolderNames{i} == start_date
        start_folder = i;
    end
    if subFolderNames{i} == end_date
        end_folder = i;
    end
end


for j=start_folder:end_folder
topLevelFolder2 = ['U:\Lab_2023\Experiment_results\QRAM\Analysis\',subFolderNames{j}] ; % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'
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
            tR = num_transits*0;
            tT = num_transits*0;
            tB = num_transits*0;
            tD = num_transits*0;
            num_transits_tot = num_transits*0;
        end

        num_transits_tot = num_transits_tot + num_transits;
        tB = tB + target_B.*(SNR>SNR_th);
        tD = tD + target_D.*(SNR>SNR_th);
        tR = tR + target_R.*(SNR>SNR_th);
        tT = tT + target_T.*(SNR>SNR_th);
    
    end
end


% t_routing_fid = tT./(tR+tT); % qrouter and SPRINT-T
t_routing_fid = tR./(tR+tT); % SPRINT-R

t_routing_pts = (tR+tT);
t_coherence = tB./(tB+tD);
t_coherence_pts = (tB+tD);

cond = (t_routing_fid>0.75).*(t_coherence>0.85).*(t_coherence_pts>20);



