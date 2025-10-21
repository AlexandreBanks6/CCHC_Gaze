%{
Author: Alexandre Banks Gadbois
Paper: https://pubmed.ncbi.nlm.nih.gov/38888820/
Description: Uses Corner-Contingent Head Compensation (CCHC) calibration
parameters to estimate Point of Gaze (POG) on new data. Requires that
"CCHCCalibration.m" has been run first to generate the calibration model
which is stored in "CCHC_Calibration.mat".
Make sure to run EyeCornerDetector.py on the data that calibration hasn't
been run on (the data to estimate POG on). Do this by setting "IS_ON_CALIB" 
to False in the EyeCornerDetector.py header.
%}

clear
clc
close all

%% Data Directory Definition
%Change the directory if needed
data_root='../../GazeData/'; %The overall data directory

MERGED_CSV_HEADER = {'timestamp', 'framecount', 'pupil_right_x', 'pupil_right_y', 'pupil_right_width', ...
                'pupil_right_height', 'pupil_right_angle', 'pupil_found_right', 'glint0_right_x', ...
                'glint0_right_y', 'glint1_right_x', 'glint1_right_y', 'glint2_right_x', 'glint2_right_y', ...
                'pupil_left_x', 'pupil_left_y', 'pupil_left_width', 'pupil_left_height', 'pupil_left_angle', ...
                'pupil_found_left', 'glint0_left_x', 'glint0_left_y', 'glint1_left_x', 'glint1_left_y', ...
                'glint2_left_x', 'glint2_left_y', 'calib_x', 'calib_y', 'calib_valid', 'right_inner_x', ...
                'right_inner_y', 'right_outer_x', 'right_outer_y', 'left_inner_x', 'left_inner_y', ...
                'left_outer_x', 'left_outer_y'};

%% Estimating POG
folder_list=dir(data_root);
dirnames={folder_list([folder_list.isdir]).name};
dirnames = dirnames(~ismember(dirnames, {'.', '..'})); % Remove '.' (current) and '..' (parent) directories from list
num_dir = length(dirnames);

for m=[1:num_dir] %Looping for all directories of data to estimate
    
    current_dir=fullfile(data_root,dirnames{m});

    %--------------------<Loading Data>-----------------
    % Load the previously saved calibration parameters
    calib_file = fullfile(current_dir,'calibration', 'CCHC_Calibration.mat');
    if ~exist(calib_file, 'file')
        disp(['!!!!Do not have calibration file for: ', dirnames{m}, '; cannot run estimation!!!!']);
        continue;
    end
    
    load(calib_file, 'model_poly', 'dist_cell', 'avg_corners','mdl_right_x','mdl_right_y', 'mdl_left_x', 'mdl_left_y');
    
    gazelog_files=dir(fullfile(current_dir,'gazelog_*.txt')); %Gets all gazelog files in this subdirectory
    %Check that we have gazelog files
    if isempty(gazelog_files)
        disp(['!!!!No gazelog files in: ',current_dir,' cannot run estimation!!!!']);
        continue;
    end


    %Looping through gazelog files, merge with eyecorner data file, and
    %estimate compensated POG
    for j=[1:length(gazelog_files)]
        gazelog_filename = gazelog_files(j).name;
        disp(['Estimating POG For: ',gazelog_filename]);
        % Extract timestamp from gazelog filename to find corresponding eyeCornerData file
        timestamp_str = extractBetween(gazelog_filename, 'gazelog_', '.txt');
        timestamp_str = timestamp_str{1};

        % Construct expected eyeCornerData filename
        eyecorner_filename = ['eyeCornerData_', timestamp_str, '.csv'];
        eyecorner_filepath = fullfile(current_dir, eyecorner_filename);
        
        %Read gazelog data
        gazelog_filepath = fullfile(current_dir, gazelog_filename);
        gazelog_data = readmatrix(gazelog_filepath);

        %Read eyeCornerData
        eyecorner_data = readmatrix(eyecorner_filepath);
        
        %------------<Merging Data>------------
        % Get framecount from gazelog (column 3) and Frame_No from eyeCornerData (column 1)
        gazelog_frames=gazelog_data(:, 3); % framecount column
        eyecorner_frames=eyecorner_data(:, 1); % Frame_No column

        % Find matching frames
        [common_frames, gazelog_idx, eyecorner_idx] = intersect(gazelog_frames, eyecorner_frames,'stable');
        if isempty(common_frames)
            disp(['!!!No matching frames found between: ', gazelog_filename, ' and ', eyecorner_filename,' cannot run POG estimation on these.']);
            merged_files{j} = [];
            continue;
        end

        % Extract relevant columns from gazelog data (matched rows only)
        gazelog_matched=gazelog_data(gazelog_idx, :);
        gazelog_matched=gazelog_matched(:, [1, 3, 13:24,34:48]); % timestamp, framecount, pupil/glint data, calib data
        
        %Extract matched eye corners and get correct order of columns as
        %expected in merged data
        eyecorner_matched=eyecorner_data(eyecorner_idx,:);
         % Original: 'Right_Inner_x', 'Right_Inner_y', 'Right_Outer_x', 'Right_Outer_y', 'Left_Outer_x', 'Left_Outer_y', 'Left_Inner_x', 'Left_Inner_y'
        % Desired: right_inner_x, right_inner_y, right_outer_x, right_outer_y, left_inner_x, left_inner_y, left_outer_x, left_outer_y
        eyecorner_matched=eyecorner_matched(:,[2:5,8:9,6:7]);
        
         %Combine two gazelog and eye corner data
        % Combine gazelog and eyeCornerData
        merged_data = [gazelog_matched, eyecorner_matched];
        %Create the merged filename
        % Create merged filename
        merged_filename = ['mergedData_', timestamp_str, '.csv'];
        merged_filepath = fullfile(current_dir, merged_filename);
        
        % Write merged data to CSV with header
        merged_table = array2table(merged_data, 'VariableNames', MERGED_CSV_HEADER);
        writetable(merged_table, merged_filepath);
        
        %------------<Estimating POG>----------
        % these are calibration parameters: model_poly, dist_cell, avg_corners,mdl_right_x, mdl_right_y, mdl_left_x, mdl_left_y
        pog_estimates = estimateCCHCPOG(merged_data, model_poly, dist_cell, avg_corners,mdl_right_x, mdl_right_y, mdl_left_x, mdl_left_y);


        %------------<Saving Results>----------
        output_file=fullfile(current_dir,['CCHC_POG_',timestamp_str,'.csv']);
        output_table=array2table(pog_estimates,'VariableNames',{'pog_right_x','isCCHC_right_x', 'pog_right_y','isCCHC_right_y',...
            'pog_left_x','isCCHC_left_x', 'pog_left_y','isCCHC_left_y'});
        writetable(output_table, output_file);
        disp(['POG Estimates Saved To: ',output_file])
    end


end



%% Function Definitions
function pog_estimates = estimateCCHCPOG(merged_data, model_poly, dist_cell, avg_corners, mdl_right_x, mdl_right_y, mdl_left_x, mdl_left_y)
    
    % Reformat data to match expected format
    reformatted_data = reformatData(merged_data);
    
    % Estimate POG
    pog_estimates = evalCCHCPOG(model_poly, reformatted_data, dist_cell, avg_corners,mdl_right_x, mdl_right_y, mdl_left_x, mdl_left_y);
end


function pog_results = evalCCHCPOG(model_cell, reformatted_data, dist_cell, avg_corners, ...
                                  comp_model_x_right, comp_model_y_right, comp_model_x_left, comp_model_y_left)
    
    left_headers = {'pg0_left_x','pg0_left_y','pg1_left_x','pg1_left_y','pg2_left_x','pg2_left_y'};
    right_headers = {'pg0_right_x','pg0_right_y','pg1_right_x','pg1_right_y','pg2_right_x','pg2_right_y'};
    check_model_right = any(ismember(model_cell(:,1), right_headers));
    check_model_left = any(ismember(model_cell(:,1), left_headers));

    [row_n, ~] = size(reformatted_data);
    pog_results = [NaN,0,NaN,0,NaN,0,NaN,0]; % [pog_right_x, is_comp_rightx,pog_right_y,is_comp_righty, pog_left_x,is_comp_leftx, pog_left_y,is_comp_lefty]
    
    for i=[1:row_n]
        curr_row = reformatted_data(i, :);
        
        % Calculate corner deltas
        del_corner_inner_x = avg_corners(1) - curr_row(14);
        del_corner_outer_x = avg_corners(3) - curr_row(16);
        del_corner_inner_y = avg_corners(2) - curr_row(15);
        del_corner_outer_y = avg_corners(4) - curr_row(17);
        
        %% Right Eye Processing
        if check_model_right
            pog_right = processSingleEye(curr_row, 2:7, right_headers, model_cell, dist_cell, ...
                                       comp_model_x_right, comp_model_y_right, ...
                                       del_corner_inner_x, del_corner_inner_y, ...
                                       del_corner_outer_x, del_corner_outer_y);
            pog_results(i, 1:4) = pog_right;
        end
        
        %% Left Eye Processing  
        if check_model_left
            pog_left = processSingleEye(curr_row, 8:13, left_headers, model_cell, dist_cell, ...
                                      comp_model_x_left, comp_model_y_left, ...
                                      del_corner_inner_x, del_corner_inner_y, ...
                                      del_corner_outer_x, del_corner_outer_y);
            pog_results(i, 5:8) = pog_left;
        end
    end
end


function pog_xy = processSingleEye(curr_row, pg_indices, headers, model_cell, dist_cell, ...
                                  comp_model_x, comp_model_y, del_corner_inner_x, del_corner_inner_y, ...
                                  del_corner_outer_x, del_corner_outer_y)
    
    pog_xy = [NaN, 0,NaN,0]; % Default return value
    
    pgs = curr_row(pg_indices);
    nan_indexs = isnan(pgs);
    
    if sum(~nan_indexs) >= 4 % At least 2 x,y pairs detected
        stripped_header = headers(~nan_indexs);
        valid_header = findPgWithAssociatedDistance(stripped_header, dist_cell);
        
        if length(valid_header) > 2
            % Find best model
            model_valid_indexes = ismember(model_cell(:,1), valid_header);
            updated_model_cell = model_cell(model_valid_indexes, :);
            [row_new, ~] = size(updated_model_cell);
            
            % Find best model pair
            cur_val = updated_model_cell{1,2};
            cur_ind = 1;
            for j = 1:2:row_new
                if updated_model_cell{j,2} > cur_val
                    cur_val = updated_model_cell{j,2};
                    cur_ind = j;
                end
            end
            
            model_x = updated_model_cell{cur_ind, 3};
            model_y = updated_model_cell{cur_ind+1, 3};
            header_x = valid_header{cur_ind};
            header_y = valid_header{cur_ind+1};
            pg_x_ind = ismember(headers, header_x);
            pg_y_ind = ismember(headers, header_y);
            
            % Get scaling factors
            [d_calib, d_curr] = findScalingFactors(dist_cell, valid_header, headers, pgs);
            
            if ~isnan(d_calib) && ~isnan(d_curr)
                % Scale PGs
                pg_x = (d_calib/d_curr) * pgs(pg_x_ind);
                pg_y = (d_calib/d_curr) * pgs(pg_y_ind);
                
                % Basic polynomial POG estimation
                [predictors_x, predictors_y] = customPolynomial(pg_x, pg_y);
                POG_x_poly = findPOG(model_x, predictors_x);
                POG_y_poly = findPOG(model_y, predictors_y);
                
                % CCHC Compensation
                predictors = [del_corner_inner_x, del_corner_inner_y, del_corner_outer_x, del_corner_outer_y];
                [comp_predictors_x, comp_predictors_y] = compensationPolynomial(predictors);
                
                POG_x_compensated = POG_x_poly;
                POG_y_compensated = POG_y_poly;
                is_comp_x=0;
                is_comp_y=0;
                if ~isempty(comp_model_x) && ~any(isnan(comp_predictors_x)) %Only enter if we have compensation model
                    del_POG_x = findCompensation(comp_model_x, comp_predictors_x);
                    POG_x_compensated = del_POG_x + POG_x_poly;
                    is_comp_x=1;
                end
                
                if ~isempty(comp_model_y) && ~any(isnan(comp_predictors_y)) %Only update POG if we have compensation model
                    del_POG_y = findCompensation(comp_model_y, comp_predictors_y);
                    POG_y_compensated = del_POG_y + POG_y_poly;
                    is_comp_y=1;
                end
                
                pog_xy = [POG_x_compensated,is_comp_x, POG_y_compensated,is_comp_y];

            end
        end
    end
end



function reformatted_data=reformatData(eval_data)
    %Returns the data to evaluate the model in the format of: 
    % frame_no, pg0_rightx, pg0_righty, ..., pg2_rightx, pg2_righty, pg0_leftx, pg0_lefty,..., pg2_leftx, pg2_lefty,
    % right_inner_x,right_inner_y,right_outer_x,right_outer_y,left_inner_x,left_inner_y,left_outer_x,left_outer_y,
    % target_x,target_y, pupil_right_x, pupil_right_y, pupil_left_x,pupil_left_y


    glintspupils_right_ind=[3,4,9,10,11,12,13,14]; %Contains the glints and pupil positions such that pupil_x,pupil_y,glint0_x,glint0_y...
    glintspupils_left_ind=[15,16,21,22,23,24,25,26];

    glintspupils_right=eval_data(:,glintspupils_right_ind);
    glintspupils_left=eval_data(:,glintspupils_left_ind);

    reformatted_data=[eval_data(:,2),glintspupils_right(:,3)-glintspupils_right(:,1),...
        glintspupils_right(:,4)-glintspupils_right(:,2),glintspupils_right(:,5)-glintspupils_right(:,1),...
        glintspupils_right(:,6)-glintspupils_right(:,2),glintspupils_right(:,7)-glintspupils_right(:,1),...
        glintspupils_right(:,8)-glintspupils_right(:,2),...
        glintspupils_left(:,3)-glintspupils_left(:,1),...
        glintspupils_left(:,4)-glintspupils_left(:,2),glintspupils_left(:,5)-glintspupils_left(:,1),...
        glintspupils_left(:,6)-glintspupils_left(:,2),glintspupils_left(:,7)-glintspupils_left(:,1),...
        glintspupils_left(:,8)-glintspupils_left(:,2),...
        eval_data(:,30:37),...
        eval_data(:,27),eval_data(:,28),glintspupils_right(:,1),glintspupils_right(:,2),...
        glintspupils_left(:,1),glintspupils_left(:,2)];

end


function valid_header=findPgWithAssociatedDistance(header,dist_cell)
valid_header=cell(0);
    for i=[1:length(dist_cell(:,1))]
        switch dist_cell{i,1}
            case 'd_01_right'
                check_pgs={'pg0_right_x','pg0_right_y','pg1_right_x','pg1_right_y'};
                check_inds=ismember(header,check_pgs);
                valid_header=[valid_header,header{check_inds}];
            case 'd_02_right'
                check_pgs={'pg0_right_x','pg0_right_y','pg2_right_x','pg2_right_y'};
                check_inds=ismember(header,check_pgs);
                valid_header=[valid_header,header{check_inds}];
            case 'd_12_right'
                check_pgs={'pg1_right_x','pg1_right_y','pg2_right_x','pg2_right_y'};
                check_inds=ismember(header,check_pgs);
                valid_header=[valid_header,header{check_inds}];
            case 'd_01_left'
                check_pgs={'pg0_left_x','pg0_left_y','pg1_left_x','pg1_left_y'};
                check_inds=ismember(header,check_pgs);
                valid_header=[valid_header,header{check_inds}];
            case 'd_02_left'
                check_pgs={'pg0_right_x','pg0_right_y','pg2_left_x','pg2_left_y'};
                check_inds=ismember(header,check_pgs);
                valid_header=[valid_header,header{check_inds}];
            case 'd_12_left'
                check_pgs={'pg1_left_x','pg1_left_y','pg2_left_x','pg2_left_y'};
                check_inds=ismember(header,check_pgs);
                valid_header=[valid_header,header{check_inds}];
        end


    end
    valid_header=unique(valid_header);


end

function [d_calib,d_curr]=findScalingFactors(dist_cell,valid_header,overall_header,pgs_only)
    %Finding distance of calibration

    if any(ismember(valid_header,'pg0_right_x')) && any(ismember(valid_header,'pg1_right_x')) && any(ismember(valid_header,'pg2_right_x'))
        dist_names={'d_01_right','d_02_right','d_12_right'};
    elseif any(ismember(valid_header,'pg0_right_x')) && any(ismember(valid_header,'pg1_right_x'))
        dist_names={'d_01_right'};
    elseif any(ismember(valid_header,'pg0_right_x')) && any(ismember(valid_header,'pg2_right_x'))
        dist_names={'d_02_right'};
    elseif any(ismember(valid_header,'pg1_right_x')) && any(ismember(valid_header,'pg2_right_x'))
        dist_names={'d_12_right'};

    elseif any(ismember(valid_header,'pg0_left_x')) && any(ismember(valid_header,'pg1_left_x')) && any(ismember(valid_header,'pg2_left_x'))
        dist_names={'d_01_left','d_02_left','d_12_left'};
    elseif any(ismember(valid_header,'pg0_left_x')) && any(ismember(valid_header,'pg1_left_x'))
        dist_names={'d_01_left'};
    elseif any(ismember(valid_header,'pg0_left_x')) && any(ismember(valid_header,'pg2_left_x'))
        dist_names={'d_02_left'};
    elseif any(ismember(valid_header,'pg1_left_x')) && any(ismember(valid_header,'pg2_left_x'))
        dist_names={'d_12_left'};

    end
    dist_ind=ismember(dist_cell(:,1),dist_names);
    if all(~dist_ind) %We don't have any corresponding distances in the calibration
        d_calib=nan;
        d_curr=nan;
    else
        dist_vec=cell2mat(dist_cell(dist_ind,2));
        d_calib=mean(dist_vec);
    
        %Finding the current inter-glint distance
    
        valid_inds=ismember(overall_header,valid_header);
        valid_pgs=pgs_only(valid_inds);
        x_vals=[];
        y_vals=[];
        for i=[1:2:length(valid_pgs)]
            x_vals=[x_vals,valid_pgs(i)];
            y_vals=[y_vals,valid_pgs(i+1)];
    
        end
        if length(x_vals)==2
            d_curr=sqrt((y_vals(1)-y_vals(2)).^2+(x_vals(1)-x_vals(2)).^2);
    
        elseif length(x_vals)==3
            diff_1=sqrt((y_vals(1)-y_vals(2)).^2+(x_vals(1)-x_vals(2)).^2);
            diff_2=sqrt((y_vals(1)-y_vals(3)).^2+(x_vals(1)-x_vals(3)).^2);
            diff_3=sqrt((y_vals(2)-y_vals(3)).^2+(x_vals(2)-x_vals(3)).^2);
            d_curr=(diff_1+diff_2+diff_3)/3;
        else
            d_curr=nan;
            d_calib=nan;
    
        end
    end
    

end


function [predictors_x,predictors_y]=customPolynomial(pg_x,pg_y)
    %predictors=[pg_x.^2,pg_x.*pg_y,pg_y.^2,pg_x,pg_y];
    %predictors_x=[pg_x,pg_y,pg_x.^2];
    %predictors_y=[pg_y,pg_x.^2,pg_x.*pg_y,pg_x.^2.*pg_y];
    predictors_x=[pg_x.^2,pg_x.*pg_y,pg_y.^2,pg_x,pg_y];
    predictors_y=[pg_x.^2,pg_x.*pg_y,pg_y.^2,pg_x,pg_y];
    %predictors_x=[pg_x.^2,pg_x.*pg_y,pg_y.^2,pg_x,pg_y];
    %predictors_y=[pg_x.^2,pg_x.*pg_y,pg_y.^2,pg_x,pg_y];

end

function del_POG=findCompensation(model,predictors)
%Generalized function to find the POG at run time
del_POG=model(1)+sum(model(2:end)'.*predictors,2,'omitnan');


end
function POG=findPOG(model,predictors)
%Generalized function to find the POG at run time
POG=model(1)+sum(model(2:end)'.*predictors,'omitnan');


end
function [predictors_x,predictors_y]=compensationPolynomial(predictors)
%Predictors are in the format: del_corner_inner_x,del_corner_inner_y, del_corner_outer_x,del_corner_outer_y,alpha

%predictors_x=[predictors(:,1),predictors(:,3),predictors(:,2),predictors(:,4),predictors(:,1).*predictors(:,2),predictors(:,3).*predictors(:,4)];
%predictors_y=[predictors(:,1),predictors(:,3),predictors(:,2),predictors(:,4),predictors(:,1).*predictors(:,2),predictors(:,3).*predictors(:,4)];
%predictors_x=[predictors(:,1).*predictors(:,2),predictors(:,3).*predictors(:,4)];
%predictors_y=[predictors(:,1).*predictors(:,2),predictors(:,3).*predictors(:,4)];
%predictors_x=[predictors(:,1).^2,predictors(:,3).^2,predictors(:,1),predictors(:,3),predictors(:,2),predictors(:,4)];
%predictors_y=[predictors(:,2).^2,predictors(:,4).^2,predictors(:,2),predictors(:,4),predictors(:,1),predictors(:,3)];
%predictors_x=[predictors(:,1:4),predictors(:,1).*predictors(:,2),predictors(:,3).*predictors(:,4)];
%predictors_y=[predictors(:,1:4),predictors(:,1).*predictors(:,2),predictors(:,3).*predictors(:,4)];
%predictors_x=[predictors(:,1:4)];
%predictors_y=[predictors(:,1:4)];
%predictors_x=[predictors(:,1:4)];
%predictors_y=[predictors(:,1:4)];
predictors_x=predictors(:,[1,3]);
predictors_y=predictors(:,[2,4]);

end


