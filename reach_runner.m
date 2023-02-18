% reach_runner
% This scripts goes through all the necessary steps required to run
% successive simulations with the Arm11x models, all to the same target.

clear all;
close all;

global first_sim

model = 'Arm11d';
load_system(model); % only necessary when model isn't already loaded

logger = find_system('Arm11d','BlockType','ToFile');
% logger is a cell array whose only contents should be a handle to the 
% block that saves the error data

tic;
for reach = 1:8
    if reach == 1
        first_sim = true;  % so params11x starts from scratch
    else
        first_sim = false;
    end
    
    params11d;  % initializing parameters
    FVname = strcat('err11d2T8r',num2str(reach));    
    set_param(logger{1},'FileName',FVname,'MatrixName',FVname);
    sim(model);
end
toc;

GenPlots;