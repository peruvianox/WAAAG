function [sortedPoints] = sortPoints(Point_xyz, PointLabels, MkrSet_type, alt_str)

% -------------------------------------------------------------------------
% SORT POINTS
% -------------------------------------------------------------------------

% Author: Kate Worster
% Date: June 6, 2011
%
% Updated 07Apr2016
% Kayla Burnim
% Compatible with Score/Sara _M structure
%
% DESCRIPTION:  Pairs the xyz coordinates of a point with point
%               labels for dynamic trials.  Compatible with Gaia, PiG, and AddRadius models.
%
% INPUTS: 
%           Point_xyz               xyz coordinates of points from C3D file
%           PointLabels             labels of all points in C3D file
%           alt_str                 alternative string if subject name was
%                                   written to 3D data fields in c3d file
%           MkrSet_type             identifies which marker set is being
%                                   used
%
% OUTPUTS: 
%           sortedPoints            cells of each point's xyz coordinates
%                                   with corresponding point label
%
% NOTE:     Popup window for alternative marker label name search turned
%           off.

%% DEBUG LOOP
% close all
% clear all
% clc
% 
% % Prompts User for Input
% [FullFileName, pathname] = uigetfile('*.c3d', 'Please select C3D file to be analyzed');
% 
% [C3Dfile] = OpenC3D(FullFileName);
% 
% % Open and Read the user specified C3D file
% % and returns: HeaderInfo, Point_xyz, PointLabels
% [HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);
% 
% % alternative string
% alt_str = '';
%
%%%%% END DEBUG LOOP %%%%%



%% DEFINE SEARCH LABELS
% extract and organize labels based upon marker set type
[search_Pt_label] = defPoints_Labels(MkrSet_type);


%% ALTERNATIVE STRING
% Check to see if an alternative string was provided
if (isempty(alt_str) == 0)    
    for t = 1:length(search_Pt_label)
        % If an alternative string was provided, concatenate with typical labels
        % provide new string labels
        search_Pt_label(t).name = strcat(alt_str, search_Pt_label(t).name);
    end
end


    
%% Find and organize Point_xyz and PointLabels 
for s = 1:length(search_Pt_label)
    % Search for each (default) point label
    default_searchPtLabel = strncmp(PointLabels, search_Pt_label(s).name,8);
    if (sum(default_searchPtLabel) == 1)
        for c1 = 1:length(PointLabels) 
            m1 = 1;
            while m1 ~= 0        
                if default_searchPtLabel(1,c1) == 1
                    % finds the column number where desired point is stored
                    ncolumn_default_pt = c1;
                    m1 = 0;
                else
                    m1 = 0;
                end
            end
        end
        default_pt = Point_xyz(:,ncolumn_default_pt,:);
        
    elseif (sum(default_searchPtLabel) == 0)
        default_pt = [];
    end
        
%     sum_default_pt = sum(default_pt(1,1,:));
        
    
    % ------ TEST DEFAULT NAME  ------ %
    if (isempty(default_pt) == 0)
        % If default search returns data, keep data
        sortedPoints(s).data(:,:,:) = default_pt;
        sortedPoints(s).name = search_Pt_label(s).name;
                
    elseif (isempty(default_pt) == 1)
        % If default label is empty, check for an alternative label
        if (isempty(search_Pt_label(s).alt_name1) == 0)
            % If there is an alternative1 label, search point data
            alt1_searchPtLabel = strncmp(PointLabels, search_Pt_label(s).alt_name1,8);
            if (sum(alt1_searchPtLabel) == 1)
                for c1 = 1:length(PointLabels) 
                    m1 = 1;
                    while m1 ~= 0        
                        if alt1_searchPtLabel(1,c1) == 1
                            % Finds the column number where desired point is stored
                            ncolumn_alt1_pt = c1;
                            m1 = 0;
                        else
                            m1 = 0;
                        end
                    end
                end
                alt1_pt = Point_xyz(:,ncolumn_alt1_pt,:);
            elseif (sum(alt1_searchPtLabel) == 0)
                alt1_pt = [];
            end
                
            % ------ TEST ALTERNATE NAME 1 ------ %
            if (isempty(alt1_pt) == 0)
                % If alternate name 1 search returns data, keep data
                sortedPoints(s).data(:,:,:) = alt1_pt;
                sortedPoints(s).name = search_Pt_label(s).alt_name1;
                
            elseif (isempty(alt1_pt) == 1)
    
                % If 1st alternate name is empty, search for 2nd alternate label
                if (isempty(search_Pt_label(s).alt_name2) == 0)
                    % Check to see if a 2nd alternative label has data
                        alt2_searchPtLabel = strncmp(PointLabels, search_Pt_label(s).alt_name2,8);
                        if (sum(alt2_searchPtLabel) == 1)
                            for c1 = 1:length(PointLabels) 
                                m1 = 1;
                                while m1 ~= 0        
                                    if alt2_searchPtLabel(1,c1) == 1
                                        % finds the column number where desired point is stored
                                        ncolumn_alt2_pt = c1;
                                        m1 = 0;
                                    else
                                        m1 = 0;
                                    end
                                end
                            end
                            alt2_pt = Point_xyz(:,ncolumn_alt2_pt,:);
                            
                        elseif sum(alt2_searchPtLabel == 0)
                            alt2_pt = [];
                        end

   
                        % ------ TEST ALTERNATE NAME 2 ------ %
                        if (isempty(alt2_pt) == 0)
                            % If alternate label 2 has data, keep
                            sortedPoints(s).data(:,:,:) = alt2_current_pt;    
                            sortedPoints(s).name = search_Pt_label(s).alt_name2;

                        elseif (isempty(alt2_pt) == 1)
                            % If can't find the point, offer user chance to type in a label name
%                             first_msg = 'Please enter another label for marker: ';  current_label = search_Pt_label(s).alt_name2;
%                             promt_msg = strcat(first_msg, current_label);
%                             prompt = {promt_msg};
%                             dlg_title = 'Unable to find marker label';
%                             num_lines = 1;
%                             def = {''};
%                             uanswer = inputdlg(prompt,dlg_title,num_lines,def);
                              uanswer = [];
                              
                            % ------ TEST USER NAME ------ %
                            if (isempty(uanswer) == 0)
                                % If user provided an alternative label to search
                                pt_label = char(uanswer);    % convert from cell to string
                                % Search for marker
                                [pt_xyz] = searchPoint(Point_xyz, PointLabels, pt_label);
                                
                                if (isempty(pt_xyz) == 0)
                                    % if there is data in user label, keep
                                    sortedPoints(s).data(:,:,:) = pt_xyz;
                                    sortedPoints(s).name = uanswer;
                                else
                                    % If no data in user's point, zero out
                                    sortedPoints(s).data(:,:,:) = search_Pt_label(s).alt_name2;
                                    sortedPoints(s).name = uanswer;
                                end
                            else
                                % If user did not provide an alternative marker label
                                % No other options to search for, so zero out
                                sortedPoints(s).data(:,:,:) = 0;
                                sortedPoints(s).data(:,:,:) = [];
                                sortedPoints(s).name = search_Pt_label(s).alt_name2;
                            end   
                            % ------ END TEST USER NAME ------ %
                        end
                        % ------ END TEST ALTERNATE NAME 2 ------ %   
                   
                        
                elseif (isempty(search_Pt_label(s).alt_name2) == 1)
                    % If 2nd alternate name field does not exist, offer user chance to type in a label name
%                     first_msg = 'Please enter another label for marker: ';  current_label = search_Pt_label(s).alt_name1;
%                     promt_msg = strcat(first_msg, current_label);
%                     prompt = {promt_msg};
%                     dlg_title = 'Unable to find marker label';
%                     num_lines = 1;
%                     def = {''};
%                     uanswer = inputdlg(prompt,dlg_title,num_lines,def);
                     uanswer = [];
                    
                    % ------ TEST USER NAME ------ %
                    if (isempty(uanswer) == 0)
                        % If user provided an alternative label to search
                        pt_label = char(uanswer);    % convert from cell to string
                        % Search for marker
                        [pt_xyz] = searchPoint(Point_xyz, PointLabels, pt_label);
                        
                        if (isempty(pt_xyz) == 0)
                            % if there is data in user label, keep
                            sortedPoints(s).data(:,:,:) = pt_xyz;
                            sortedPoints(s).name = uanswer; 
                        else
                            % If no data in user's point, zero out
                            sortedPoints(s).data(:,:,:) = [];
                            sortedPoints(s).name = search_Pt_label(s).alt_name2;
                        end
                    else
                        % If user did not provide an alternative marker label
                        % No other options to search, so zero out
                        sortedPoints(s).data(:,:,:) = 0;
                        sortedPoints(s).data(:,:,:) = [];
                        sortedPoints(s).name = search_Pt_label(s).alt_name2;
                    end
                    % ------ END TEST USER NAME ------ %
                end
             
            end   % ------ END TEST ALTERNATE NAME 1 ------ %       
        
                
        elseif (isempty(search_Pt_label(s).alt_name1) == 1)  
            % If there is no alternative label, offer user chance to type in a label name
%             first_msg = 'Please enter another label for marker: ';  current_label = search_Pt_label(s).name;
%             promt_msg = strcat(first_msg, current_label);
%             prompt = {promt_msg};
%             dlg_title = 'Unable to find marker label';
%             num_lines = 1;
%             def = {''};
%             uanswer = inputdlg(prompt,dlg_title,num_lines,def);
            uanswer = [];
            
            if (isempty(uanswer) == 0)
                % If user provided an alternative label to search
                pt_label = char(uanswer);    % convert from cell to string
                % Search for marker
                [pt_xyz] = searchPoint(Point_xyz, PointLabels, pt_label);
                
                if (isempty(pt_xyz) == 0)
                    % if there is data in user label, keep
                    sortedPoints(s).data(:,:,:) = pt_xyz;
                    sortedPoints(s).name = uanswer; 
                else
                    % If no data in user's point, zero out
                    sortedPoints(s).data(:,:,:) = search_Pt_label(s).alt_name1;
                    sortedPoints(s).name = uanswer;
                end
                 
            else
                % If user did not provide an alternative marker label
                % No other options to search for, so zero out
               % sortedPoints(s).data(:,:,:) = [];
                sortedPoints(s).name = search_Pt_label(s).alt_name1;
            end

        end  % ------ END ALTERNATE NAME 1 ------ %       
            
    end  % ------ END TEST DEFAULT NAME ------ %       

end  % ------ END ENTIRE SEARCH LOOP ------ %       



