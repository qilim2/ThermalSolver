function [epoch,vS,vBody,eclipse] = input_file_parser_R01(input_file)
% Reads the input file given path and name of file.
% Gets:
% Epoch of each time step
% Vector from s/c to Sun in s/c body frame
% Vector from s/c to nearby celestial body in s/c body frame

% Open the file for reading
fid = fopen(input_file, 'r');

% Read the file line by line
ii = 1;
while ~feof(fid)
    tline = fgetl(fid); % Get next line
    % Read and store information
    if contains(tline,'output_All')
        fgetl(fid); tline = fgetl(fid); % Go down 2 lines
    end
    values = str2num(tline); % Convert string to vector
    if ~isempty(values) % Skip empty lines
        data(ii,:) = values;
    end
    ii = ii + 1;
end

% Close the file
fclose(fid);

epoch = data(2:end,1);
vS = data(2:end,2:4).*1000;
vBody = data(2:end,5:7).*1000;
eclipse = data(2:end,8);
end