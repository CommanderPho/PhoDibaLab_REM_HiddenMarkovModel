classdef SpeculatedUnitInfo
    %SPECULATED_UNIT_INFO defines cell type ({'pyramidal', 'contaminated', 'interneurons'}) colors: 
    %   Detailed explanation goes here
    enumeration
          pyramidal
          contaminated
          interneurons
    end

    properties (Constant)
      classColors = [0.8, 0.5, 0.1; 0.5, 0.1, 0.1; 0.0, 0.7, 0.7];
      classNames = {'pyramidal','contaminated','interneurons'};
      classCutoffValues = [0 4 7 9];
    end % end Constant properties block

    methods
        function obj = SpeculatedUnitInfo()
            %SPECULATED_UNIT_INFO Construct an instance of this class
            %   Detailed explanation goes here
        end
    end

    methods (Static)
        function [speculated_unit_type, speculated_unit_contamination_level, SpeculatedUnitInfo] = unitQualityToCellType(unit_quality)
            % unitQualityToCellType: converted from standalone function 'fnUnitQualityToCellType.m' from Diba project
            % unit_quality: active_processing.spikes.quality
	        num_units = length(unit_quality);
        
	        speculated_unit_type = discretize(unit_quality, ...
		        SpeculatedUnitInfo.classCutoffValues, ...
		        'categorical', SpeculatedUnitInfo.classNames);

	        speculated_unit_contamination_level = zeros([num_units 1]);
	        for i = 1:num_units
                temp.curr_type_index = grp2idx(speculated_unit_type(i));
                temp.curr_type_string = SpeculatedUnitInfo.classNames{temp.curr_type_index};
		        if (strcmpi(temp.curr_type_string, 'pyramidal'))
			        speculated_unit_contamination_level(i) = (active_processing.spikes.quality(i) - 1);
			        % value will be either [0, 1, 2, 3] 
		        elseif (strcmpi(temp.curr_type_string, 'contaminated'))
			        %%% speculated_unit_contamination_level(i) = (active_processing.spikes.quality(i) - 5); % A value of 5 corresponds to minimal contamination for a contaminated cell, and 7 is the maximum
			        %%% value will be either [0, 1, 2]
			        % Keep contamination level relative to pyramidal to allow easier filtering
			        speculated_unit_contamination_level(i) = (active_processing.spikes.quality(i) - 1); % A value of 5 corresponds to minimal contamination for a contaminated cell, and 7 is the maximum
			        % value will be either [4, 5, 6]
		        elseif (strcmpi(temp.curr_type_string, 'interneurons'))
			        speculated_unit_contamination_level(i) = (active_processing.spikes.quality(i) - 8); % A value of 8 corresponds to minimal contamination for an interneuron, and 9 is the maximum
			        % value will be either [0, 1]   
		        else
			        error('invalid!');
                end        
	        end
        end
       
%       function out = setgetVar(data)
%          persistent Var;
%          if nargin
%             Var = data;
%          end
%          out = Var;
%       end

    end % end static method block


end

