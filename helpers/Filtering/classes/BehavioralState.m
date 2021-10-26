classdef BehavioralState
    %BEHAVIORALSTATE Summary of this class goes here
    %   Detailed explanation goes here
    
   properties (Constant)
      classColors = [0.5, 0.5, 1.0
               0.7, 0.7, 1.0
               1.0, 0.7, 0.7
               1.0, 0.0, 0.0];
      classNames = {'nrem', 'rem', 'quiet', 'active'};
      classValues = [1:length(BehavioralState.classNames)];
    end % end Constant properties block

    methods
        function obj = BehavioralState()
            %BEHAVIORALSTATE Construct an instance of this class
            %   Detailed explanation goes here
        end
    end

end

