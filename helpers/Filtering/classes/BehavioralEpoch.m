classdef BehavioralEpoch
    %BEHAVIORALEPOCH Summary of this class goes here
    %   Detailed explanation goes here
  
   properties (Constant)
      classColors = [0.0, 0.5, 0.0
               0.2, 1.0, 0.2
               0.0, 0.2, 0.0];
      classNames = {'pre_sleep', 'track', 'post_sleep'};
      classValues = [1:length(BehavioralEpoch.classNames)];
    end % end Constant properties block

    methods
        function obj = BehavioralEpoch()
            %BehavioralEpoch Construct an instance of this class
            %   Detailed explanation goes here
        end
    end

end

