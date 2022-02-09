% Notebook summary: this notebook is divided into parts.
% 1) The first part grabs the data and organizes bioSound struct into a
% matrix.
% 2) The second part performs PCA on the bioSound Matrix to find components
% that explain the most variance.
% 3) The third part projects the bioSound matrix onto the principle components
% and then performs gaussian mixture modeling to find unsupervised
% clusters. We then color-code by Deaf vs Control group to see if the
% clusters correspond to treatment. 
% 4) The fourth part tries to find the 
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

% Get deaf, and sex labels
deaf_or_control = get_deaf_or_control(compiled_data.batIDs);

% Change boolean to do analysis on all or subsets of data
do_male_and_female = 1;
do_male_only = 0;
do_female_only = 0;

if do_male_and_female
    % Reduce space using PCA so that we can visualize the data
    [coeff,score,latent,tsquared,explained] = plot_pca_variance_explained(compiled_data.piezo_matrix, "Piezo, Male and Female");
    [coeff2,score2,latent2,tsquared2,explained2] = plot_pca_variance_explained(compiled_data.raw_matrix, "Raw, Male and Female");

    % Project matrix to lower-dim PCA space. 
    do_pca_and_project_data_onto_top_PCs(compiled_data.piezo_matrix, coeff, deaf_or_control.male_and_female, "Piezo, Male and Female")
    do_pca_and_project_data_onto_top_PCs(compiled_data.raw_matrix, coeff2, deaf_or_control.male_and_female, "Raw, Male and Female")
end

if do_male_only
    % Reduce space using PCA so that we can visualize the data
    data_piezo_mo = compiled_data.piezo_matrix(deaf_or_control.male_bool,:);
    data_raw_mo = compiled_data.raw_matrix(deaf_or_control.male_bool,:);
    %keyboard()
    [coeff_mo,score_mo,latent_mo,tsquared_mo,explained_mo] = plot_pca_variance_explained(data_piezo_mo, "Piezo, Male Only");
    [coeff2_mo,score2_mo,latent2_mo,tsquared2_mo,explained2_mo] = plot_pca_variance_explained(data_raw_mo, "Raw, Male Only");

    % Project matrix to lower-dim PCA space. 
    PCA_GMM_struct_mo = do_pca_and_project_data_onto_top_PCs(data_piezo_mo, coeff_mo, deaf_or_control.male_only, "Piezo, Male Only");
    PCA_GMM_struct2_mo = do_pca_and_project_data_onto_top_PCs(data_raw_mo, coeff2_mo, deaf_or_control.male_only, "Raw, Male Only");
    
    % This showed interesting clusters, Look into the weights
    [sorted_coeffs_mo, sorted_idx_matrix_mo] = sort(coeff_mo, 'descend');
    bioSound_labels = get_bioSound_labels_from_indices(sorted_idx_matrix_mo);    
end

if do_female_only
    % Reduce space using PCA so that we can visualize the data
    data_piezo_fe = compiled_data.piezo_matrix(deaf_or_control.female_bool,:);
    data_raw_fe = compiled_data.raw_matrix(deaf_or_control.female_bool,:);
    %keyboard()
    [coeff_fe,score_fe,latent_fe,tsquared_fe,explained_fe] = plot_pca_variance_explained(data_piezo_fe, "Piezo, Female Only");
    [coeff2_fe,score2_fe,latent2_fe,tsquared2_fe,explained2_fe] = plot_pca_variance_explained(data_raw_fe, "Raw, Female Only");

    % Project matrix to lower-dim PCA space. 
    do_pca_and_project_data_onto_top_PCs(data_piezo_fe, coeff_fe, deaf_or_control.female_only, "Piezo, Female Only")
    do_pca_and_project_data_onto_top_PCs(data_raw_fe, coeff2_fe, deaf_or_control.female_only, "Raw, Female Only")
end

%% Investigate spectrograms, and histogram of features near mean of clusters
% Find points closest to cluster means
[sorted_distances, sorted_idx, sorted_element_IDs] = find_points_closest_to_cluster_means(PCA_GMM_struct_mo, compiled_data, deaf_or_control, 'Males_only');

% Using the sorted distances and sorted idx, return the spectrogram. This
% function saves the spectrograms too
grab_spectrograms(sorted_element_IDs)

%% Fit and compare Linear Mixed Effect Models
% Make table from data sources
tbl = make_table_from_matrices(compiled_data.piezo_matrix, compiled_data.batIDs, deaf_or_control.male_and_female, score);

% Do LMEMs
lme_sex_deaf_interaction = fitlme(tbl, "PCScore~sex*Deaf_or_Control+(1|BatID)");
lme_sex_deaf_no_interaction = fitlme(tbl, "PCScore~sex+Deaf_or_Control+(1|BatID)");
lme_sex = fitlme(tbl, "PCScore~sex+(1|BatID)");
lme_deaf = fitlme(tbl, "PCScore~Deaf_or_Control+(1|BatID)");

% Compare LMEMs
test_deaf_to_sexPlusdeaf = compare(lme_sex, lme_sex_deaf_no_interaction);
test_sex_to_sexPlusdeaf = compare(lme_deaf, lme_sex_deaf_no_interaction);
test_interaction_to_noInteraction = compare(lme_sex_deaf_no_interaction, lme_sex_deaf_interaction);

% Do additional tests
lme_sex_deaf_entropytime = fitlme(tbl, "PCScore~sex+Deaf_or_Control+entropytime+(1|BatID)");
lme_sex_deaf_entropytime_stdtime = fitlme(tbl, "PCScore~sex+Deaf_or_Control+entropytime+stdtime+(1|BatID)");
lme_sex_entropytime = fitlme(tbl, "PCScore~sex+entropytime+(1|BatID)");
lme_deaf_entropytime = fitlme(tbl, "PCScore~Deaf_or_Control+entropytime+(1|BatID)");

% Compare LMEMs
test_deaf_to_sexPlusdeafPlusentropytime  = compare(lme_sex_entropytime, lme_sex_deaf_entropytime);
test_sex_to_sexPlusdeafPlusentropytime  = compare(lme_deaf_entropytime, lme_sex_deaf_entropytime);
test_entropytime_to_baseline = compare(lme_sex_deaf_no_interaction, lme_sex_deaf_entropytime);
test_stdtime_to_entropytime_and_baseline = compare(lme_sex_deaf_no_interaction, lme_sex_deaf_entropytime);


%% Make histograms of all of the biosound features in the male-only, female-only 
% run linear mixed effects model and control for identity
make_histogram_bioSound_features(data_piezo_mo, deaf_or_control.male_only, "Piezo, Male Only")
make_histogram_bioSound_features(data_piezo_fe, deaf_or_control.female_only, "Piezo, Female Only")


%% Helper functions for doing Linear Mixed Effect Models
function tbl = make_table_from_matrices(main_data, BatIDs, Deaf_or_not, PCScores)
    % Turn the matrix into a table
    tbl = array2table(main_data);
    
    % Name the table variables
    PAF = list_of_predefined_acoustical_features();
    PAF_plus_these = horzcat(PAF,{'Pk2', 'second_V', 'meanF0'});    % Recall that some of the variables were post-calculated, include these
    tbl.Properties.VariableNames = PAF_plus_these;
    
    % Add the identifiers such as BatID, Sex, and Deaf
    %keyboard()
    tbl.BatID = categorical(vertcat(BatIDs{:}));
    tbl.Deaf_or_Control = categorical(Deaf_or_not);
    tbl.sex = categorical(get_bat_sex_cells(BatIDs)');
    tbl.PCScore = PCScores(:,1);
    
end

%% Helper functions for the historgrams
function make_histogram_bioSound_features(data, deaf_or_control, title_affix)
    % Make histograms for each feature, compare deaf and control in same
    % plot
    
    % Get bioSound dictionary for feature labels
    bioSound_dict = get_bioSound_dict();
    
    % Plot each feature
    for i=1:width(data)
        deaf_data = data(deaf_or_control,i);
        control_data = data(~deaf_or_control,i);
        
        feature = bioSound_dict(i);
        figure
        hold on
        h1 = histogram(deaf_data);
        h2 = histogram(control_data);
        [p, h, stats] = ranksum(deaf_data, control_data);
        f_title = {"Histogram of feature " + feature, 'Wilcox-test pvalue: ' + string(p), title_affix};
        title(f_title);
        
        legend({'deaf','control'})
        
        % Save image
        save_as = strcat("Histogram_of_" + feature, '_pvalue_' + string(p), title_affix, '.jpg');
        file_save_path = fullfile( 'C:\Users\tobias\Documents\Alvince\DeafSalineGroup151_subset', save_as );
        saveas(gcf, file_save_path)
        hold off
        %keyboard()
    end
end

%% Helper functions to find features near mean of clusters and plotting histograms
function [sorted_distances, sorted_idx, sorted_element_IDs] = find_points_closest_to_cluster_means(PCA_GMM_struct,compiled_data, deaf_or_control, which_group)
% This function iterates through each cluster mean (usually only the top 
% two principle components) and finds which data points are closest to those
% means. Then, this function returns the sorted distances, sorted indices of those
% points, and the sorted element_IDs which are strings that can be used to
% grab spectrogram data.

    % Grab cluster means
    cluster_means = PCA_GMM_struct.GMModel.mu;
    data_projected_onto_PCs = PCA_GMM_struct.data_projected_onto_PCs;
    % For each cluster mean
    for i=1:height(cluster_means)
        % Find the euclidean distance of each point.
        distances(:,i) = sqrt((cluster_means(i,1)-data_projected_onto_PCs(:,1)).^2 + (cluster_means(i,2)-data_projected_onto_PCs(:,2)).^2);
    end
    
    % Sort descending
    [sorted_distances, sorted_idx] = sort(distances);    
    
    % Sort element IDs, these IDs are the string patterns to grab already
    % processed spectogram
    if strcmp(which_group , 'Males_only')
        %keyboard()
        males_only_ids = compiled_data.elementIDs(deaf_or_control.male_bool);
        sorted_element_IDs = males_only_ids(sorted_idx);
    end
end

function images = PDFtoImg(pdfFile, save_as)
    import org.apache.pdfbox.*
    import java.io.*
    filename = fullfile(pdfFile);
    jFile = File(filename);
    document = pdmodel.PDDocument.load(jFile);
    pdfRenderer = rendering.PDFRenderer(document);
    count = document.getNumberOfPages();
    images = [];
    for ii = 1:count
        bim = pdfRenderer.renderImageWithDPI(ii-1, 300, rendering.ImageType.RGB);
        images = [images (filename + "-" +"Page" + ii + ".png")];
        tools.imageio.ImageIOUtil.writeImage(bim, save_as, 300);
    end
    document.close()
end

function grab_spectrograms(sorted_element_IDs) 
    path_to_data = 'Y:\users\JulieE\DeafSalineGroup151\';
    vocExtracts_affix = 'audiologgers\VocExtracts';
   
    % Iterate over the two PCs
    for i=1:width(sorted_element_IDs)
        % Grab only the top 5 spectrograms closest to the cluster means
        for j = 1:5 %width(sorted_element_IDs)
            %keyboard
            subdir = strcat('20', extractBetween(sorted_element_IDs(j,i), 'DeSa_', '_') );
            subdir = strrep(subdir, '_Raw', '_Piezo');
            spectrogram_file_name = strcat(sorted_element_IDs(j,i), '.pdf');
            spectrogram_file_name = strrep(spectrogram_file_name, '_Raw', '_Piezo');
            spectrogram_file_path = fullfile(path_to_data,subdir,vocExtracts_affix,spectrogram_file_name);
            
            % Also grab the wav file
            wav_file_path = strrep(spectrogram_file_path,'.pdf','.wav');
            %open(spectrogram_file_path{1});
            
            % Save file
            save_as = strcat('Top ', string(j) ,' datapoint closest to cluster ', string(i), '.jpg');
            file_save_path = fullfile( 'C:\Users\tobias\Documents\Alvince\DeafSalineGroup151_subset', save_as );
            wave_wav_as = strrep(file_save_path, '.jpg','.wav');
            %keyboard
            try
                PDFtoImg(spectrogram_file_path{1}, file_save_path);
                %keyboard
                % Copy over wav file
                copyfile(wav_file_path{1}, wave_wav_as)
            catch
                keyboard
            end
            
        end
    end
    
end

%% General helper functions
function BatLabels = get_bat_labels()
    % This function returns a data structures that will help in
    % classification labels and controlling for variables
    BatLabels.name = {'11648', '14461', '14463', '14464', '65696', '71043', '71047', '71351', '71353', '71354'};
    BatLabels.sex = {'F', 'M', 'F', 'M', 'F', 'M', 'F', 'M', 'F', 'F'};
    BatLabels.deaf = [0 0 1 0 0 1 1 1 0 1];
    BatLabels.name_deaf_dict = containers.Map(BatLabels.name,BatLabels.deaf);
    BatLabels.name_F_or_M_dict = containers.Map(BatLabels.name, BatLabels.sex);
end

function bat_sex_cells = get_bat_sex_cells(Bat_IDs)
    BatLabels = get_bat_labels();
    bat_sex_cells = {};
    
    for i=1:length(Bat_IDs)
        %keyboard()
        bat_sex_cells{i} = BatLabels.name_F_or_M_dict(Bat_IDs{i}{1});
    end
    
end

function deaf_or_control_struct = get_deaf_or_control(batIDs)
    % This function returns boolean arrays so that when doing clustering we
    % can group the data and control for variables like gender
    % Output labels
    %   male_and_female: boolean 1=deaf, 0=control
    %   male_only: boolean 1=deaf, 0=control
    %   female_only: boolean 1=deaf, 0=control
    %   male_bool: boolean 1=male, 0=female
    %   female_bool: boolean 1=female, 0=male
    BatLabels = get_bat_labels();
    deaf_or_control_struct.male_and_female = false(length(batIDs),1);
    deaf_or_control_struct.male_and_female_bool = true(length(batIDs),1);
    deaf_or_control_struct.male_only = false(length(batIDs),1);
    deaf_or_control_struct.male_bool = false(length(batIDs),1);
    deaf_or_control_struct.female_only = false(length(batIDs),1);
    deaf_or_control_struct.female_bool = false(length(batIDs),1);
    
    for i=1:length(batIDs)     
        %keyboard()
        deaf_or_control_struct.male_and_female(i) = BatLabels.name_deaf_dict(batIDs{i}{1});
        if BatLabels.name_F_or_M_dict(batIDs{i}{1}) == 'M'
            deaf_or_control_struct.male_bool(i) = 1;
            if BatLabels.name_deaf_dict(batIDs{i}{1}) == 1 
                deaf_or_control_struct.male_only(i) = 1;
            end
        elseif BatLabels.name_F_or_M_dict(batIDs{i}{1}) == 'F'
            deaf_or_control_struct.female_bool(i) = 1;
            if BatLabels.name_deaf_dict(batIDs{i}{1}) == 1 
                deaf_or_control_struct.female_only(i) = 1;
            end
        end
        
    end
    
    % Check that female and male sums equal the total population
    if (sum(deaf_or_control_struct.male_bool) + sum(deaf_or_control_struct.female_bool)) ~= length(batIDs)
        warning = 'Something is wrong. The counts dont add up.'
        keyboard()
    end
    deaf_or_control_struct.male_only = deaf_or_control_struct.male_only(deaf_or_control_struct.male_bool);
    deaf_or_control_struct.female_only = deaf_or_control_struct.female_only(deaf_or_control_struct.female_bool);
    %keyboard()
end

function bioSound_dict = get_bioSound_dict()
% This function returns a dictionary of bioSound indices mapped to their
% feature names
PAF = list_of_predefined_acoustical_features();
PAF_plus_these = horzcat(PAF,{'Pk2', 'second_V', 'meanF0'}); % This order was gotten from get_predefined_acoustical_features(BioSoundCalls)

%keyboard()
bioSound_dict = containers.Map(1:length(PAF_plus_these), PAF_plus_these);
end

function bioSound_labels = get_bioSound_labels_from_indices(idx_matrix)
    % This function accepts idx matrix (Integers) of bioSound and returns the appropriate
    % bioSound labels (strings)
    bioSound_dict = get_bioSound_dict();
    bioSound_labels = cell(size(idx_matrix));

    % Relabel each element in cell
    for i=1:length(idx_matrix)
        idx_bool = idx_matrix == i;
        [bioSound_labels{idx_bool}] = deal(bioSound_dict(i));
    end

end

%% Helper functions for unsupervised clustering and plotting
function [coeff,score,latent,tsquared,explained] = plot_pca_variance_explained(X, title_affix)
    % Plotting reference https://www.mathworks.com/matlabcentral/answers/713038-variance-explained-pca
    [coeff,score,latent,tsquared,explained] = pca(X);%, 'NumComponents', 2);
    figure
    hold on
    bar(explained)
    %keyboard()
    plot(1:numel(explained), cumsum(explained), 'o-', 'MarkerFaceColor', 'r')
    yyaxis right
    h = gca;
    h.YAxis(2).Limits = [0 100];
    h.YAxis(2).Color = h.YAxis(1).Color;
    h.YAxis(2).TickLabel = strcat(h.YAxis(2).TickLabel, '%');
    xlabel('Number of Principle Components')
    ylabel('Variance Explained')
    title(['Number of Principle Components and Variance Explained of Predefomed Acoustic Features: ' title_affix])
    
    % Save photo
    %keyboard()
    file_save_path = fullfile( 'C:\Users\tobias\Documents\Alvince\DeafSalineGroup151_subset', strcat('PCA_exp_', title_affix, '.jpg') );
    saveas(gcf, file_save_path);
    hold off
end

function PCA_GMM_struct = do_pca_and_project_data_onto_top_PCs(X, coeff, deaf_or_control, title_affix)
    data_projected_onto_PCA = X*coeff(:,1:2); % Project data onto PCA coeffs

    % Fit Gaussian mixture model
    GMModel = fitgmdist(data_projected_onto_PCA, 2, 'Replicates',10);

    % Show gaussian contour scatterplot
    figure
    group = deaf_or_control; % [ones(height(data_projected_onto_PCA),1)]; % data labels
    h = gscatter(data_projected_onto_PCA(:,1),data_projected_onto_PCA(:,2),group, 'rk');
    hold on
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
    g = gca;
    fcontour(gmPDF,[g.XLim g.YLim])
    title(['{\bf Scatter Plot and Fitted Gaussian Mixture Contours} ' title_affix])
    legend(h,'Bat Group 1','Bat Group 2')
    xlabel('PC 1')
    ylabel('PC 2')
    
    % Save photo
    file_save_path = fullfile( 'C:\Users\tobias\Documents\Alvince\DeafSalineGroup151_subset', strcat('Scatter', title_affix, '.jpg') );
    saveas(gcf, file_save_path);
    hold off
    
    PCA_GMM_struct.data_projected_onto_PCs = data_projected_onto_PCA;
    PCA_GMM_struct.GMModel = GMModel;
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

function [matrix_raw, matrix_piezo, cell_of_IDs, cell_of_element_IDs] = get_bioSound_data_as_matrix(loaded_200_file)
    % This function grabs only the 1 x 1 numerical PAF data like in Julie's
    % 2016 Zebra Finch paper. 
    
    matrix_raw = [];
    matrix_piezo = [];
    
    cell_of_IDs = {};
    cell_of_element_IDs = {};
    
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
            cell_of_element_IDs{i,1} = name1;
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
                    [raw_matrix_tmp, piezo_matrix_tmp, batIDs, elementIDs] = get_bioSound_data_as_matrix(loaded_200_file);

                    % Append table to compiled table
                    if file_count == 1
                        compiled_data.raw_matrix = raw_matrix_tmp;
                        compiled_data.piezo_matrix = piezo_matrix_tmp;
                        compiled_data.batIDs = batIDs;
                        compiled_data.elementIDs = elementIDs;
                        
                    else
                        compiled_data.raw_matrix = [compiled_data.raw_matrix; raw_matrix_tmp];
                        compiled_data.piezo_matrix = [compiled_data.piezo_matrix; piezo_matrix_tmp];    
                        compiled_data.batIDs = [compiled_data.batIDs; batIDs];
                        compiled_data.elementIDs = [compiled_data.elementIDs; elementIDs];
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