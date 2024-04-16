function [overview, MESHGRIDS_11, SURFACES_11, SURFACES_21, THERMOPHYSICAL_11, THERMOOPTICAL_11, ELEMENTS_11, CONDUCTANCES_11, HEATLOADS_11] = bulk_data_parser_R02(filename)
    % Parses data inside the bulk data file and converts to MATLAB vars.
    % Version 2.0 completed 10/26/2023

    % Open the file for reading
    fid = fopen(filename, 'r');
    
    % Read the file line by line
    ii = 1;
    while ~feof(fid)
        tline = fgetl(fid); % Get next line
        
        %fprintf('Round %i\n',ii)
        ii = ii + 1;

        % Find matrix sizes and preallocate matrices
        if contains(tline,'OVERVIEW')
            fgetl(fid); tline = fgetl(fid); % Go down 2 lines
            overview = str2num(tline); % Create overview matrix
            %fprintf('Overview: %s\n',mat2str(overview)) % Check output

            MESHGRIDS_11 = zeros(overview(1),6);
            ELEMENTS_11 = zeros(overview(2),9);
            THERMOPHYSICAL_11 = zeros(overview(4),4);
            THERMOOPTICAL_11 = zeros(overview(5),3);
            CONDUCTANCES_11 = zeros(overview(1),overview(1));
            SURFACES_11 = zeros(overview(3),15);
            SURFACES_21 = {};
        end

        % Read and store nodes
        if contains(tline,'NODES')
            fgetl(fid); tline = fgetl(fid); % Go down 2 lines
            for jj = 1:overview(1)
                MESHGRIDS_11(jj,:) = str2num(tline);
                tline = fgetl(fid);
            end
            %fprintf('%s\n',mat2str(MESHGRIDS_11));
        end
        
        % Read and store elements
        if contains(tline,'ELEMENTS')
            fgetl(fid); tline = fgetl(fid); % Go down 2 lines
            for jj = 1:overview(2)
                ELEMENTS_11(jj,:) = str2num(tline);
                tline = fgetl(fid);
            end
            %fprintf('%s\n',mat2str(ELEMENTS_1));
        end

        % Read and store thermophysical properties
        if contains(tline,'THERMOPHYSICAL PROPERTIES')
            fgetl(fid); tline = fgetl(fid); % Go down 2 lines
            for jj = 1:overview(4)
                THERMOPHYSICAL_11(jj,:) = str2num(tline);
                tline = fgetl(fid);
            end
            %fprintf('%s\n',mat2str(THERMOPHYSICAL_11));
        end

        % Read and store thermo-optical properties
        if contains(tline,'THERMO-OPTICAL PROPERTIES')
            fgetl(fid); tline = fgetl(fid); % Go down 2 lines
            for jj = 1:overview(5)
                THERMOOPTICAL_11(jj,:) = str2num(tline);
                tline = fgetl(fid);
            end
            %fprintf('%s\n',mat2str(THERMOOPTICAL_11));
        end
        
        % Read and store surfaces properties
        if contains(tline,'SURFACES PROPERTIES')
            fgetl(fid); tline = fgetl(fid); % Go down 2 lines
            for jj = 1:overview(3)
                SURFACES_11(jj,:) = str2num(tline);
                tline = fgetl(fid);
            end
            %fprintf('%s\n',mat2str(SURFACES_11));
        end
        
        % Read and store surfaces definitions
        if contains(tline,'SURFACES DEFINITIONS')
            fgetl(fid); tline = fgetl(fid); % Go down 2 lines
            for jj = 1:overview(3)
                ID = str2num(tline);
                SURFACES_21{jj}.ID = ID(1);
                SURFACES_21{jj}.nodes = ID(2:end);
                elems = find(ELEMENTS_11(:,1)==ID);
                SURFACES_21{jj}.elems = elems;
                tline = fgetl(fid);
            end
            %fprintf('%s\n',mat2str(SURFACES_21));
        end

        % Read and store conductances
        if contains(tline,'CONDUCTANCES')
            fgetl(fid); tline = fgetl(fid); % Go down 2 lines
            for jj = 1:overview(1)
                CONDUCTANCES_11(jj,:) = str2num(tline);
                tline = fgetl(fid);
            end
            %fprintf('%s\n',mat2str(CONDUCTANCES_11));
        end

        % Read and store heatloads
        if contains(tline,'USER DEFINED HEAT LOADS')
            fgetl(fid); tline = fgetl(fid); % Go down 2 lines
            HEATLOADS_11 = str2num(tline);
            %fprintf('%s\n',mat2str(HEATLOADS_11));
        end
    end
    
    % Close the file
    fclose(fid);
end