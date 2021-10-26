% Pho Hale
% 2021-10-25
% Registers the classses in this folder (currently hardcoded) as appDesigner components so they can be used in the designer.

componentParentFolder = '/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/helpers/ToMoveToPhoLib/CustomAppDesignerComponents';
generatedResourcesFolder = fullfile(componentParentFolder, 'resources');
% componentPath = '/Users/pho/repo/NGP Rotations Repos/PhoDibaLab_DataAnalysis/helpers/ToMoveToPhoLib/CustomAppDesignerComponents/cmpColorSelector.m';
componentPath = fullfile(componentParentFolder, 'cmpColorSelector.m');

appdesigner.customcomponent.configureMetadata(componentPath); % Opens the interactive edit metadata dialog

% For the UI component to appear in the App Designer Component Library, you must add the folder containing the class file and generated resources folder to the MATLAB path.
addpath(genpath(generatedResourcesFolder));


