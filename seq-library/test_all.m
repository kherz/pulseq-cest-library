% Set the main directory where the subfolders with .m files are stored
main_dir = pwd;

% Get a list of all subfolders in the main directory
subfolders = dir(main_dir);
subfolders = subfolders([subfolders.isdir]); % Only keep folders

% Loop through each subfolder and run the .m file inside
for i = 1:numel(subfolders)
    folder_name = subfolders(i).name;
    
    % Skip over the "." and ".." directories
    if strcmp(folder_name, '.') || strcmp(folder_name, '..')
        continue;
    end
    
    % Change to the subfolder directory
    folder_path = fullfile(main_dir, folder_name);
    cd(folder_path);
    
    % Run the .m file in the subfolder
    run(folder_name);
    
    % call standard sim
    if 0 % not neccesary in new m files, as it is done there already
    M_z = simulate_pulseqcest(seq_filename,'../../sim-library/WM_3T_default_7pool_bmsim.yaml');
    writematrix(M_z', ['M_z_' seq_filename '.txt']);
    end
    % Change back to the main directory
    cd(main_dir);
end

disp('all seq files genererated');
