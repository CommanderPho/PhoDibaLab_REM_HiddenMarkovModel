classdef Pho
    %PHO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        InputArgs
    end
    
    methods
        function obj = Pho(optionalInputArg)
            %Pho Construct an instance of this class
            %   Detailed explanation goes here
            if ~exist('optionalInputArg','var')
                optionalInputArg = '';
            end
                
            obj.InputArgs = optionalInputArg;
        end
        
        
        function [] = insertAtCursorNow(self)
                %insertAtCursorNow: Inserts the text at the current cursor location:
                
                if exist('snippet'), snippet = struct('snippet',snippet); end;      % Save the variable SNIPPET, if it exists
             
                % Type your text here!
                snippet.txt = self.snippetStringRepresentation;
                
                snippet.txt = strrep(snippet.txt, '\n', sprintf('%c',10));          % Replace '\n' by a new line
                
                if verLessThan('matlab', '8.1.0')
                    snippet.activeEditor = editorservices.getActive;                % Get the active document in the editor
                else
                    snippet.activeEditor = matlab.desktop.editor.getActive;         % Get the active document in the editor
                end
                
                snippet.activeEditor.JavaEditor.insertTextAtCaret(snippet.txt);     % Insert text at the current position
                
                if isfield(snippet,'snippet')                                       % Delete SNIPPET (or replace it by its original value)
                    snippet = snippet.snippet;
                else
                    clear snippet
                end;
        end
        
    end
    
    methods (Static)
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
        
        
        function [] = insertFunctionHeaderTemplateNow()
			insertTemplateHeader;
        end
        
        
        function [] = newScript(name, body)
            % Creates a new document from a template (no template functionality yet though).
            if ~exist('name','var')
               name = 'DEFAULT';
               error('no name provided to newScript(...)');
            end
            if ~exist('body','var')
                body = '';
            end
            
           matlab.desktop.editor.newDocument(body);
           % TODO: finish
            
        end
        
        
    end
end

