function [foundFunctionDefinitions, matlabCodeLines] = fnParseMatlabScriptForFunctionDefinitions(matlabCodeText)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% matlabCodeText: a string with multiline text containing matlab code

    % regex.functionDefinitionLine = 'function (?:\[(.*)\]\s*=\s*)([A-Za-z0-9_]+)\((.*)\).*';
    regex.functionDefinitionLine = 'function (?:\[(?<output_args>.*)\]\s*=\s*)(?<fn_name>[A-Za-z0-9_]+)\((?<input_args>.*)\)';
    match_variable_names = {'output_args','fn_name','input_args'};
    
    if ~iscell(matlabCodeText)
       % Split into lines if provided as a single multiline string
       matlabCodeLines = matlab.desktop.editor.textToLines(matlabCodeText);
    else
       matlabCodeLines = matlabCodeText;
    end

    numLines = length(matlabCodeLines);

   
    % tokenNames = regexp(matlabCodeText, regex.functionDefinitionLine, 'match');
%     foundFunctionDefinitions = struct();
    numFoundFunctionDefinitions = 0;

    % Iterate through each line and try to match for a function defininiton
    for i = 1:numLines
        tokenNames = regexp(matlabCodeLines{i}, regex.functionDefinitionLine, 'names');
        if isempty(tokenNames)
            % Nothing found for this line: Keep moving on.
            %% TODO: maybe find the matching 'end' blocks this way?

        else
            numFoundFunctionDefinitions = numFoundFunctionDefinitions + 1;
            foundFunctionDefinitions(numFoundFunctionDefinitions).lineNumber = i;
            for fieldIndex = 1:length(match_variable_names)
                if isfield(tokenNames, match_variable_names{fieldIndex})        
                    foundFunctionDefinitions(numFoundFunctionDefinitions).(match_variable_names{fieldIndex}) = tokenNames.(match_variable_names{fieldIndex});
                else
                    foundFunctionDefinitions(numFoundFunctionDefinitions).(match_variable_names{fieldIndex}) = ''; 
                end % end if field

            end % end for fields

        end


    %     parsed_fn.output_args = tokenNames.output_args;
    %     parsed_fn.fn_name = tokenNames.fn_name;
    %     parsed_fn.input_args = tokenNames.input_args;

    end

end

