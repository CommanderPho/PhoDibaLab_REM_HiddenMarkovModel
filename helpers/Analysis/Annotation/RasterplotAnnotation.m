classdef RasterplotAnnotation
    %RasterplotAnnotation A user-created annotation for a raster plot
    %   Based off of "ExplicitlyTypedUserAnnotation" in Matlab-Pho-Helper-Tools
    
    properties
        timestamp 
        endTimestamp % Optional end timestamp if we're referring to a range of time, NaN otherwise.
        creationDateTime
        modifiedDateTime
        creatingUser
        comment
        typeName % Type of annotation
        unitIDs = {}; % A list of units being referred to, or an empty list of it's a global note.
    end
    
    methods
        function obj = RasterplotAnnotation(typeName, startTimestamp, endTimestamp, comment, unitIDs, creatingUser, creationDatetime, modifiedDatetime)
            %RasterplotAnnotation Construct an instance of this class
            %   Detailed explanation goes here
            if ~exist('creationDatetime','var')
                creationDatetime = datetime('today');
            end
            
            if ~exist('creatingUser','var')
                creatingUser = 'Unknown';
            end
            
            if ~exist('comment','var')
                comment = '';
            end    
            
            if ~exist('modifiedDatetime','var')
                modifiedDatetime = creationDatetime;
            end   
            
            if ~exist('endTimestamp','var') | isempty(endTimestamp)
               endTimestamp = NaN; 
            end
            
            if ~exist('unitIDs','var')
                unitIDs = {};
            end
            
            obj.typeName = typeName;
            obj.creatingUser = creatingUser;
            obj.creationDateTime = creationDatetime;
            obj.modifiedDateTime = modifiedDatetime;
            obj.timestamp = startTimestamp;
            obj.endTimestamp = endTimestamp;
            obj.comment = comment;
            obj.unitIDs = unitIDs;
            
        end
        
        function obj = modifyComment(obj, comment)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.comment = comment;
            obj.modifiedDateTime = datetime('today');
        end
        
        
        
    end
end

