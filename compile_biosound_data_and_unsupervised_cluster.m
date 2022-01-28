% Script summary: this script performs supervised classification on the Saline vs. Deaf
% acoustic dataset, which has already been processed by the BioSound
% pipeline.
%
% This script has two main parts
%   1) Data compilation. The data need to be loaded from the server and
%   organized into a single struct
%   2) Classification. Once the data are organized we perform
%   classification
%
% Script author: Alvince Pongos, alvince_pongos@berkeley.edu

%% Data Compilation

data_root_path = 'Y:\users\JulieE\DeafSalineGroup151\';

% Select only those subdirectories which have been processed and noted in 
processed_log = read_log_file('Y:\users\JulieE\DeafSalineGroup151\RecOnlyLogDeafSalWhat.txt');

% As described by the "processed_log", grab the *200.mat files in each of
% the subdirectories
compiled_data = get_compiled_data(processed_log, data_root_path);

%% Perform unsupervised classification

% Reduce space using PCA so that we can visualize the data
[coeff,score,latent,tsquared,explained] = plot_pca_variance_explained(compiled_data.piezo_matrix);

% Get deaf or control labels
deaf_or_control = get_deaf_or_control(compiled_data.batIDs);

% Project matrix to lower-dim PCA space. We choose first 8 coefficients
% since this explains 99% of data variance
do_pca_and_project_data_onto_top_PCs(compiled_data.piezo_matrix, coeff, deaf_or_control)

%% Perform supervised classification


%% Unsupervised classification

%% General helper functions
function BatLabels = get_bat_labels()
    BatLabels.name = {'11648', '14461', '14463', '14464', '65696', '71043', '71047', '71351', '71353', '71354'};
    BatLabels.sex = {'F', 'M', 'F', 'M', 'F', 'M', 'F', 'M', 'F', 'F'};
    BatLabels.deaf = [0 0 1 0 0 1 1 1 0 1];
    BatLabels.name_deaf_dict = containers.Map(BatLabels.name,BatLabels.deaf);
end

function deaf_or_control = get_deaf_or_control(batIDs)
    BatLabels = get_bat_labels()
    deaf_or_control = false(length(batIDs),1);
    
    for i=1:length(batIDs)     
        %keyboard()
        deaf_or_control(i) = BatLabels.name_deaf_dict(batIDs{i}{1});
    end
end

%% Helper functions for unsupervised clustering and plotting
function [coeff,score,latent,tsquared,explained] = plot_pca_variance_explained(X)
    % Plotting reference https://www.mathworks.com/matlabcentral/answers/713038-variance-explained-pca
    [coeff,score,latent,tsquared,explained] = pca(X);%, 'NumComponents', 2);
    figure
    hold on
    bar(explained)
    plot(1:numel(explained), cumsum(explained), 'o-', 'MarkerFaceColor', 'r')
    yyaxis right
    h = gca;
    h.YAxis(2).Limits = [0 100];
    h.YAxis(2).Color = h.YAxis(1).Color;
    h.YAxis(2).TickLabel = strcat(h.YAxis(2).TickLabel, '%');
    xlabel('Number of Principle Components')
    ylabel('Variance Explained')
    title('Number of Principle Components and Variance Explained of Predefomed Acoustic Features')
    hold off
end

function do_pca_and_project_data_onto_top_PCs(X, coeff, deaf_or_control)
    data_projected_onto_PCA = X*coeff(:,1:2);

    % Fit Gaussian mixture model
    GMModel = fitgmdist(data_projected_onto_PCA, 2);

    % Show gaussian contour scatterplot
    figure
    group = deaf_or_control; % [ones(height(data_projected_onto_PCA),1)]; % data labels
    h = gscatter(data_projected_onto_PCA(:,1),data_projected_onto_PCA(:,2),group, 'rk');
    hold on
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
    g = gca;
    fcontour(gmPDF,[g.XLim g.YLim])
    title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
    legend(h,'Bat Group 1','Bat Group 2')
    xlabel('PC 1')
    ylabel('PC 2')
    hold off
end

%% Helper functions for Data Compilation
function struct_out  =  read_log_file(path2txt)
    fid = fopen(path2txt);
    processed_log = textscan(fid,'%s%s%s');
    struct_out.subject = processed_log{1,1}(2:end);
    struct_out.date = fix_dates( processed_log{1,2}(2:end) );
    struct_out.time = processed_log{1,3}(2:end);

    fclose(fid);
end

function boolean_out = check_if_valid_date(date_name)
    skip_these = {'20191123','20191124','20200122'};
    boolean_out = false;
    if any(strcmp(skip_these,date_name))
       boolean_out = false;
    else
       boolean_out = true;
    end
end

function fixed_dates = fix_dates(dates)
    % The names in the logger do not match the names in the directory,
    % prefix the '20' on each date to match
    fixed_dates = {};
    for i=1:length(dates)
        fixed_dates{i, 1} = ['20' dates{i,1}];
    end
end

% TODO: Ask how to normalize table sizes for features with variable sized
% arrays
function [table1_raw, table2_piezo] = get_bioSound_data_as_table(loaded_200_file)
    % Initialize tables
    table1_raw = table;
    table2_piezo = table;
    
    % Append each BioSound output to output tables
    for i = 1:length(loaded_200_file.BioSoundCalls)   
        % Grab the file names
        [filepath1,name1,ext1] = fileparts(loaded_200_file.BioSoundFilenames{i, 1});
        [filepath2,name2,ext2] = fileparts(loaded_200_file.BioSoundFilenames{i, 2});
        
        % Grab respective data
        table1_raw_tmp = struct2table(loaded_200_file.BioSoundCalls{i, 1},'AsArray',true);
        table2_piezo_tmp = struct2table(loaded_200_file.BioSoundCalls{i, 2},'AsArray',true);
        
        % Add batID as a field. BatID needs to be grabbed from the file
        % name.
        table1_raw_tmp.batID = extractBetween(name1,"Bat","_");
        table2_piezo_tmp.batID = extractBetween(name2,"Bat","_");
        
        % Add filename just in case
        table1_raw_tmp.filename = name1;
        table2_piezo_tmp.filename = name2;
        
        % Append data to end of table
        table1_raw = [table1_raw; table1_raw_tmp];
        table2_piezo = [table2_piezo; table2_piezo_tmp];
    end
end

function PAF = list_of_predefined_acoustical_features()
    % List of predefined acoustical features (PAF)
    PAF = {'meantime','stdtime','skewtime','kurtosistime','entropytime', ...
                                       'maxAmp','rms', ... % Temporal envelope features
                                       'q1','q2','q3',  ...
                                       'meanspect','stdspect','skewspect', ...
                                       'kurtosisspect','entropyspect',  ...% Spectral envelope features
                                       'meansal','minfund','maxfund','cvfund', ...% Fundamental features
                                        };
end

function matrix_of_PAF = get_predefined_acoustical_features(BioSoundCalls)

    % List of predefined acoustical features (PAF)
    predefined_acoustical_features = list_of_predefined_acoustical_features();
    
   % Reduce struct to only PAF
   all_field_names = fieldnames(BioSoundCalls);
   toRemove = all_field_names(~ismember(all_field_names,predefined_acoustical_features));
   
   only_PAF_struct = rmfield(BioSoundCalls,[toRemove]);
   
   % Almost done, but some of the PAF need to be calculated
   try
       only_PAF_struct.Pk2 = nanmean(BioSoundCalls.f0_2);
       only_PAF_struct.second_V = sum(~isnan(BioSoundCalls.f0_2)) / length(BioSoundCalls.f0_2);
       only_PAF_struct.meanF0 = nanmean(BioSoundCalls.f0);
   catch
        keyboard()
   end
   
   % Convert struct to cell to matrix
   cell_of_PAF = struct2cell(only_PAF_struct);
   cell_of_PAF(cellfun(@isempty,cell_of_PAF)) = {NaN};
   matrix_of_PAF = cell2mat(cell_of_PAF)';
   
end

function [matrix_raw, matrix_piezo, cell_of_IDs] = get_bioSound_data_as_matrix(loaded_200_file)
    % This function grabs only the 1 x 1 numerical PAF data like in Julie's
    % 2016 Zebra Finch paper. 
    
    matrix_raw = [];
    matrix_piezo = [];
    
    cell_of_IDs = {};
    
    recording_count = 1;
    % Append each BioSound output to output tables
    for i = 1:height(loaded_200_file.BioSoundCalls)   
        % Check that there are data in these files, otherwise continue
       % try
        if length(fieldnames(loaded_200_file.BioSoundCalls{i,1})) < 18
            continue
        end
        % catch
        %    keyboard()
        %end
        
        % Grab the file names
        try
            [filepath1,name1,ext1] = fileparts(loaded_200_file.BioSoundFilenames{i, 1});
            % Append to cell of IDs
            cell_of_IDs{i,1} = extractBetween(name1,"Bat","_");
        catch
            warning = 'Error no BioSoundFilename';
            keyboard()
            %cell_of_IDs{i,1} = 'Error no BioSoundFilename';
        end
        
        
        
        % Grab respective data
        matrix_raw_tmp = get_predefined_acoustical_features(loaded_200_file.BioSoundCalls{i, 1});
        matrix_piezo_tmp = get_predefined_acoustical_features(loaded_200_file.BioSoundCalls{i, 2});
        
        if recording_count == 1
            matrix_raw = matrix_raw_tmp;
            matrix_piezo = matrix_piezo_tmp;
        else
            matrix_raw = [matrix_raw; matrix_raw_tmp];
            matrix_piezo = [matrix_piezo; matrix_piezo_tmp];
        end
        
        if height(matrix_raw) ~= height(cell_of_IDs)
          keyboard() 
       end
        
        recording_count = recording_count + 1;
    end
   
end

function boolean_if_enough_data = check_if_data_has_all_PAF(loaded_200_file)
% This function checks if all PAF are in struct
    PAF = list_of_predefined_acoustical_features();
    
    boolean_if_enough_data = false(height(loaded_200_file.BioSoundCalls),1);
    
    for i = 1:height(loaded_200_file.BioSoundCalls)
        %keyboard()
        if isstruct(loaded_200_file.BioSoundCalls{i,1})
            if all(ismember(PAF, fieldnames(loaded_200_file.BioSoundCalls{i,1}) ) ) 
                boolean_if_enough_data(i) = true;
            end
        end
    end
end

function compiled_data = get_compiled_data(processed_log, data_root_path)
    % Iterate through the log to select which subdirectories have ready
    % data
    compiled_data = struct;
    file_count = 1;
    for i = 1:length(processed_log.date)
        
        % Find all *200.mat files in subdirectory
        subdir_name = processed_log.date{i};
        path_to_200_mat_files = fullfile(data_root_path, subdir_name, 'audiologgers'); % 200.mat files are where the extracted features are.
        all_200_mat_files = dir(fullfile(path_to_200_mat_files,'*200.mat'));
        
        % Now iterate through all of the valid *200.mat files
        if ~isempty(all_200_mat_files) && check_if_valid_date(subdir_name)
            % Each *200.mat file contains multiple recordings, iterate
            % through those
            for j = 1:length(all_200_mat_files)
                % Turn struct data into table for matlab classification
                % methods
                file_200_name = all_200_mat_files(j).name
                loaded_200_file = load(fullfile(path_to_200_mat_files,file_200_name));
                
                if ~isfield(loaded_200_file,'BioSoundCalls')
                    warning = 'This file does not have BioSoundCalls'
                    %keyboard()
                    continue
                end
                % Remove empty rows of BioSoundCalls and BioSoundFiles
                non_empty_bool = ~cellfun('isempty',loaded_200_file.BioSoundCalls);
                enough_data_bool = check_if_data_has_all_PAF(loaded_200_file);
                non_empty_bool = non_empty_bool & enough_data_bool; % Take the logical intersection so you keep only good data
                
                some_empty_data1 = loaded_200_file.BioSoundCalls(:,1);
                some_empty_data2 = loaded_200_file.BioSoundCalls(:,2);
                loaded_200_file.BioSoundCalls = [some_empty_data1(non_empty_bool(:,1)), some_empty_data2(non_empty_bool(:,2)) ];
                
                some_empty_filenames1 = loaded_200_file.BioSoundFilenames(:,1);
                some_empty_filenames2 = loaded_200_file.BioSoundFilenames(:,2);
                loaded_200_file.BioSoundFilenames = [some_empty_filenames1(non_empty_bool(:,1)), some_empty_filenames2(non_empty_bool(:,2)) ];
                
                % Only pass if struct not empty
                if isstruct(loaded_200_file) 
                    [raw_matrix_tmp, piezo_matrix_tmp, batIDs] = get_bioSound_data_as_matrix(loaded_200_file);

                    % Append table to compiled table
                    if file_count == 1
                        compiled_data.raw_matrix = raw_matrix_tmp;
                        compiled_data.piezo_matrix = piezo_matrix_tmp;
                        compiled_data.batIDs = batIDs;
                    else
                        compiled_data.raw_matrix = [compiled_data.raw_matrix; raw_matrix_tmp];
                        compiled_data.piezo_matrix = [compiled_data.piezo_matrix; piezo_matrix_tmp];    
                        compiled_data.batIDs = [compiled_data.batIDs; batIDs];
                    end
                    
                    % Debug if matrix and cells don't line up!
                    if height(compiled_data.raw_matrix) ~= height(compiled_data.batIDs)
                        keyboard()
                    end
                    file_count = file_count + 1;
                end
            end
        end
    end
end