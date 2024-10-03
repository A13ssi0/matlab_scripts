close all; clearvars; clc;

load("dataset_no_processing_v3.mat")

data = nan(size(data));
events = nan(size(events));

subject = 'filippo';        % filippo francesco

%%
if strcmpi(subject,'francesco')     % standard folder concatenation
    calPath =  '/media/alessio/Elements/whi2/whi2_filippo_data/';
    evalPath = '/home/alessio/francesco_training_2024/evaluation_rhlh/';
    
    calFold = dir([calPath,'*f1']);
    evalFold = dir([evalPath,'*f1']);
    
    folders = struct2table([calFold; evalFold]);
    folders = table2struct(sortrows(folders, 'name'));
elseif strcmpi(subject,'filippo')   % all files in the same folder
    path  = '/media/alessio/Elements/whi2/';
    fold = dir([path]);
    folders = struct2table(fold);
    folders = table2struct(sortrows(folders, 'name'));
end

%%

%events
%     'FILE_NO'
%     'RUN    '
%     'MODE   '
%     'TYP    '
%     'POS    '
%     'DUR    '
%     'EOG    '

%mode_names
%     'UNKNOWN:-1'
%     'OFFLINE:0 '
%     'ONLINE:1  '
%     'CONTROL:2 '


%%
idx_data_last_file = 0;
idx_events_last_file = 0;
file_no = 1;
n_hit = [];
n_miss = [];
for actual_folder = folders'
    files = dir([actual_folder.folder, '/', actual_folder.name,'/*.gdf']);
    for fl = files'
        file_path = [fl.folder,'/',fl.name];
        disp(fl.name)
        [s,h, warning_flag] = sload(file_path);
        if ~strcmpi(h.Label(1),'eeg:1')
            disp('Skipperd')
            continue
        end
%         if warning_flag
%             disp('DELETING -----------------')
%             delete(file_path)
%         end
        
        data_length = size(s,1);
        data(idx_data_last_file+(1:data_length),:) = s(:,1:16);

        filenames(file_no,:) = pad(file_path,size(filenames,2));
        
        t_events = h.EVENT;

        idx_781 = find(t_events.TYP==781);
        if isempty(idx_781)
            disp(file_path)
        end
        fl_events.POS = t_events.POS(idx_781);
        fl_events.DUR = t_events.DUR(idx_781);
        fl_events.TYP = t_events.TYP(idx_781-1);

        events_length = length(fl_events.POS);
        events(idx_events_last_file+(1:events_length), 1) = file_no-1; %python start with 0
        events(idx_events_last_file+(1:events_length), 2) = file_no-1;

        if contains(file_path,'calibration')
            events(idx_events_last_file+(1:events_length), 3) = 0;
            n_hit = cat(1,n_hit,nan);
            n_miss = cat(1,n_miss,nan);

        elseif contains(file_path,'evaluation')
            events(idx_events_last_file+(1:events_length), 3) = 1;
            n_hit = cat(1,n_hit,sum(t_events.TYP==897));
            n_miss = cat(1,n_miss,sum(t_events.TYP==898));

        elseif contains(file_path,'control')
            events(idx_events_last_file+(1:events_length), 3) = 2;
            n_hit = cat(1,n_hit,nan);
            n_miss = cat(1,n_miss,nan);

        else
            events(idx_events_last_file+(1:events_length), 3) = -1;
            n_hit = cat(1,n_hit,nan);
            n_miss = cat(1,n_miss,nan);
        end

        events(idx_events_last_file+(1:events_length), 4) = fl_events.TYP;
        events(idx_events_last_file+(1:events_length), 5) = fl_events.POS+idx_data_last_file;
        events(idx_events_last_file+(1:events_length), 6) = fl_events.DUR;
        events(idx_events_last_file+(1:events_length), 7) = 0;

        idx_data_last_file = idx_data_last_file + data_length;
        file_no = file_no+1;
        idx_events_last_file = idx_events_last_file + events_length;
    end
end

[row,~] = find(isnan(data),1,"first");
data(row:end,:) = [];

filenames(file_no:end,:) = [];

[row,~] = find(isnan(events),1,"first");
events(row:end,:) = [];

%%
%save('filippo_dataset.mat','data',"events","filenames","column_names","mode_names","channel_labels","n_miss","n_hit")
