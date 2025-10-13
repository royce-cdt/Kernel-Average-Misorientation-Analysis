%% MTEX Demo for Grain + GND analysis
% Demonstrates loading of Bruker based h5 files
% Highlights how to make this similar with the Bruker based eSprit analysis
% Two phase analysis and texture segmentation
%
% Written by Dr Ben Britton
% Edited by Orán Magan (o.magan24@imperial.ac.uk)
% Original data & MTEX code derived from work with Simon Wyatt & Jim Hickey
%
% E: b.britton@imperial.ac.uk
% T: @bmatb
% W: http://www.expmicromech.com
%
% Working with:
% MATLAB Version: 9.3.0.713579 (R2017b)
% MTEX-5.0.beta3
% List of file paths
% List of file paths
% List of file paths (cell array of string scalars)
% 1) File list
%%
clear fnames
fnames = {
  
  "F:\EBSD\Hole12_EBSD\LocAB\locab.ctf"
  "F:\EBSD\Hole12_EBSD\LocA\EBSDHOLE12_100625_LOC1\locA_1.ctf",
  "F:\EBSD\Hole12_EBSD\LocA\EBSDHOLE12_240725_A_REDO\locA_2.ctf",
  "F:\EBSD\Hole12_EBSD\Zeropoint\zero.ctf",
  "F:\EBSD\Hole12_EBSD\LocB\locB.ctf",
  "F:\EBSD\Hole12_EBSD\LocC\locC.ctf", 
  "F:\EBSD\Hole12_EBSD\LocD\locD.ctf"
  };

% Check that all files actually exist
for i = 1:numel(fnames)
    if ~isfile(fnames{i})
        warning('Missing file: %s', fnames{i});
    end
end


% 2) MTEX setup
addpath('C:\Users\oranm\Downloads\mtex-6.0.0\mtex-6.0.0');
startup_mtex; opengl software;
setMTEXpref('xAxisDirection','west'); setMTEXpref('zAxisDirection','outOfPlane');

if ~exist('fnames', 'var') || isempty(fnames)
    fnames = {
    
      "F:\EBSD\Hole12_EBSD\LocAB\locab.ctf"
      "F:\EBSD\Hole12_EBSD\LocA\EBSDHOLE12_100625_LOC1\locA_1.ctf",
      "F:\EBSD\Hole12_EBSD\LocA\EBSDHOLE12_240725_A_REDO\locA_2.ctf",
      "F:\EBSD\Hole12_EBSD\Zeropoint\zero.ctf",
      "F:\EBSD\Hole12_EBSD\LocB\locB.ctf",
      "F:\EBSD\Hole12_EBSD\LocC\locC.ctf", 
      "F:\EBSD\Hole12_EBSD\LocD\locD.ctf"
      };
end


% 3) Loop
for x = 1:numel(fnames)
    fprintf('Processing file %d/%d: %s\n', x, numel(fnames), fnames{x});
    
    % Confirm loop variable is alive
    whos
    disp(['Loop counter x = ' num2str(x)]);
    
    % Current filename
    fname = fnames{x};            % string scalar
    disp(fname)
    assert(isfile(fname), 'File does not exist alert!');
    fname_char = char(fname);     % convert to char for EBSD.load
    
    % Create results folder
    [~, base, ~] = fileparts(fname_char);
    resultsdir = fullfile(pwd, 'results', base);
    if ~exist(resultsdir,'dir')
        mkdir(resultsdir);
    end

    
    % Load EBSD
    [ebsd, header] = EBSD.load(fname_char, 'convertSpatial2EulerReferenceFrame');
    disp(ebsd(1:5))
    

    %filter out any values with MAD greater than a threshold
    % madThreshold = 0.9; % degrees
    % ebsd = ebsd(ebsd.mad < madThreshold);

    % Clean non-indexed points
    %ebsd = fill(ebsd);


    
    
    y_min = min(ebsd.y);  % most negative, bottom of scan
    y_max = max(ebsd.y);  % least negative (or zero), top of scan
    
    % Compute the y cutoff point to remove the bottom 40%
    %y_cutoff = y_min + 0.4 * (y_max - y_min);  % boundary between bottom 40% and upper 60%
    
    % Keep only the upper 60%
    %ebsd = ebsd(ebsd.y >= y_cutoff);  % keep points above this threshold
    
    % Open and read header
    fid = fopen(fname, 'r');
    rawLines = {};
    line = fgetl(fid);
    while ischar(line)
        rawLines{end+1} = line;
        if startsWith(line, 'Phase')  % header ends just before data table
            break
        end
        line = fgetl(fid);
    end
    fclose(fid);
    
    % Now parse XStep and YStep
    xStep = [];
    yStep = [];
    
    for i = 1:numel(rawLines)
        line = strtrim(rawLines{i});
        if contains(line, 'XStep')
            tokens = regexp(line, 'XStep\s+([\d\.]+)', 'tokens');
            if ~isempty(tokens)
                xStep = str2double(tokens{1}{1});
            end
        elseif contains(line, 'YStep')
            tokens = regexp(line, 'YStep\s+([\d\.]+)', 'tokens');
            if ~isempty(tokens)
                yStep = str2double(tokens{1}{1});
            end
        end
    end

    disp(['XStep = ', num2str(xStep), ' µm']);
    disp(['YStep = ', num2str(yStep), ' µm']);
    
    % fid = fopen(fname, 'r');
    % rawLines = {};
    % line = fgetl(fid);
    % while ischar(line)
    %     rawLines{end+1} = line;
    %     line = fgetl(fid);
    % end
    % fclose(fid);
    % 
    % xStep = [];
    % yStep = [];
    % 
    % for i = 1:numel(rawLines)
    %     line = strtrim(rawLines{i});
    %     if contains(line, 'XSTEP')
    %         tokens = regexp(line, 'XSTEP:\s+([\d\.]+)', 'tokens');
    %         if ~isempty(tokens)
    %             xStep = str2double(tokens{1}{1});
    %         end
    %     elseif contains(line, 'YSTEP')
    %         tokens = regexp(line, 'YSTEP:\s+([\d\.]+)', 'tokens');
    %         if ~isempty(tokens)
    %             yStep = str2double(tokens{1}{1});
    %         end
    %     end
    % end
    % 
    % disp(['XStep = ', num2str(xStep), ' µm']);
    % disp(['YStep = ', num2str(yStep), ' µm']);
    
    
    % Plot the same plots as Bruker eSprit
    
    
    % phase map
    % 
    % Fix the colours
    % set(ebsd('titanium-alpha'), 'color', [1 0 0]);
    % set(ebsd('titanium-beta'), 'color', [0 1 0]);
    % 
    % Phase map
    % FigH.phase = figure;
    % plot(ebsd, 'phase');
    % title('Phase Map');
    % saveas(FigH.phase, fullfile(resultsdir, 'PhaseMap.png'));
    % close(FigH.phase)
    
    % % IPF maps
    % phase1 = 'titanium-alpha';
    % phase2 = 'titanium-beta';
    % 
    % oM1 = ipfHSVKey(ebsd(phase1));
    % oM2 = ipfHSVKey(ebsd(phase2));
    % 
    % oM1.inversePoleFigureDirection = xvector;
    % FigH.IPFX = figure;
    % plot(ebsd(phase1), oM1.orientation2color(ebsd(phase1).orientations))
    % title('IPFX');
    % hold on;
    % plot(ebsd(phase2), oM2.orientation2color(ebsd(phase2).orientations))
    % saveas(FigH.IPFX, fullfile(resultsdir, 'IPF_X.png'));
    % close(FigH.IPFX)
    % 
    % oM1.inversePoleFigureDirection = yvector;
    % FigH.IPFY = figure;
    % plot(ebsd(phase1), oM1.orientation2color(ebsd(phase1).orientations))
    % hold on;
    % plot(ebsd(phase2), oM2.orientation2color(ebsd(phase2).orientations))
    % title('IPFY');
    % saveas(FigH.IPFY, fullfile(resultsdir, 'IPF_Y.png'));
    % close(FigH.IPFY)
    % 
    % oM1.inversePoleFigureDirection = zvector;
    % FigH.IPFZ = figure;
    % plot(ebsd(phase1), oM1.orientation2color(ebsd(phase1).orientations))
    % hold on;
    % plot(ebsd(phase2), oM2.orientation2color(ebsd(phase2).orientations))
    % title('IPFZ');
    % saveas(FigH.IPFZ, fullfile(resultsdir, 'IPF_Z.png'));
    % close(FigH.IPFZ)
    % 
    % % Plot the IPF colorkey
    % FigH.Key = figure('Color', [1 1 1]);
    % 
    % direction = vector3d.Z;
    % 
    % s1 = subplot(1, 2, 1);
    % plot(oM1, 'Parent', s1);
    % title('Ti-alpha');
    % mtexColorbar('title', 'IPF Ti-α', 'location', 'southOutside');
    % 
    % s2 = subplot(1, 2, 2);
    % plot(oM2, 'Parent', s2);
    % title('Ti-beta');
    % mtexColorbar('title', 'IPF Ti-β', 'location', 'southOutside');
    % 
    % saveas(FigH.Key, fullfile(resultsdir, 'IPF_Colorkey.png'));
    % close(FigH.Key)

    
    %%
    % KAM
    % Step 1: Gridify and subset
    ebsd_alpha = ebsd('titanium-alpha').gridify;
    
    % % Step 2: Initial grain calculation
    % thresholdInDegrees = 2.5;
    % % Compute grains from the cleaned, gridified ebsd_alpha
    % [grains, ebsd_alpha.grainId] = calcGrains(ebsd_alpha, 'angle', thresholdInDegrees*degree);
    % 
    % grains = grains(grains.area > 1);
    % grains = smooth(grains, 2);
    % 
    % Step 2: Mask invalid points (NaN orientations)
    valid_mask = ~isnan(ebsd_alpha.orientations);
    
    % Step 3: Convolve the mask with a 3x3 kernel to count valid neighbors
    % Padding with 0s ensures edges aren't falsely counted
    neighbor_kernel = [1 1 1; 1 0 1; 1 1 1];  % 8-neighbor connectivity, radius = 1
    
    %neighbor_kernel = ones(5,5); % radius = 2 connectivity
    %neighbor_kernel(3,3) = 0;   % don’t count the center pixel

    valid_neighbor_count = conv2(double(valid_mask), neighbor_kernel, 'same');
    
    % Step 4: Compute KAM using MTEX
    kam = ebsd_alpha.KAM('threshold', 2.5*degree, 'radius', 1);
    
    % Step 5: Mask KAM where fewer than x neighbors are valid
    kam(valid_neighbor_count(:) < 3) = NaN;
    
    % Convert to degrees
    kam_thresholded_deg = kam / degree;
    
    % Plot
    FigH.KAM = figure;
    plot(ebsd_alpha, kam_thresholded_deg, 'micronbar', 'on');
    setColorRange([0, 2]);
    mtexColorbar('title', 'KAM (°)');
    mtexColorMap('parula');
    title('KAM (≥ 4 valid neighbors only)');
    saveas(FigH.KAM, fullfile(resultsdir, 'KAM_Thresholded.png'));
    close(FigH.KAM);


        % ==== BETA PHASE KAM CALCULATION AND EXPORT ====

    % Step 1: Extract and gridify beta phase
    ebsd_beta = ebsd('titanium-beta').gridify;
    
    % % Step 2: Optional MAD filtering (if needed)
    % madThreshold_beta = 1.4;  % adjust if needed
    % ebsd_beta = ebsd_beta(ebsd_beta.mad < madThreshold_beta);
    
    % Step 3: Compute valid mask and neighbors
    valid_mask_beta = ~isnan(ebsd_beta.orientations);
    neighbor_kernel = [1 1 1; 1 0 1; 1 1 1];  % 8-neighbor connectivity, radius = 1
    
    % neighbor_kernel = ones(5,5); % radius = 2 connectivity
    % neighbor_kernel(3,3) = 0;   % don’t count the center pixel

    valid_neighbor_count_beta = conv2(double(valid_mask_beta), neighbor_kernel, 'same');
    
    % Step 4: Compute KAM
    kam_beta = ebsd_beta.KAM('threshold', 2.5*degree, 'radius', 1);
    kam_beta(valid_neighbor_count_beta(:) < 3) = NaN;
    kam_beta_deg = kam_beta / degree;
    
    % Step 5: Plot KAM map
    FigH.KAM_beta = figure;
    plot(ebsd_beta, kam_beta_deg, 'micronbar', 'on');
    setColorRange([0, 2]);
    mtexColorbar('title', 'KAM β (°)');
    mtexColorMap('parula');
    title('KAM - β Phase');
    saveas(FigH.KAM_beta, fullfile(resultsdir, 'KAM_Beta.png'));
    close(FigH.KAM_beta);  

        
    
    % % Grain Orientation Spread (GOS)
    % 
    % % Step 1: Calculate GOS for each grain
    % meanOrientations = grains.meanOrientation;
    % numGrains = length(grains);
    % GOS_values = zeros(numGrains, 1);
    % 
    % for i = 1:numGrains
    %     grainOrientations = ebsd_alpha(grains(i));
    %     misorientations = angle(grainOrientations.orientations * inv(meanOrientations(i)));
    %     GOS_values(i) = mean(misorientations) / degree;
    % end
    % 
    % % Step 2: Assign GOS to grain object
    % grains.GOS = GOS_values;
    % grains.GOS(grains.GOS > 5) = 0;
    % 
    % % Step 3: Map GOS from grains to each EBSD point using grainId
    % grainId_map = ebsd_alpha.grainId;
    % gos = nan(size(grainId_map)); % preallocate same size as ebsd_alpha
    % 
    % for i = 1:numGrains
    %     % Only assign GOS to points belonging to grain i
    %     gos(grainId_map == grains.id(i)) = GOS_values(i);
    % end
    % 
    % % filter any gos point above 5 degree
    % gos(gos > 5) = 0;
    
    % Step 4: Prepare data for export
    
    % Get coordinates as vectors
    x = ebsd_alpha.x(:);  % flatten to column vector, microns
    y = ebsd_alpha.y(:);  % flatten to column vector, microns
    
    % Check if coordinates are in pixels (assuming >200 means pixels)
    if any(x > 200) || any(y > 200)
        % Convert from pixels to microns by multiplying with step sizes
        x = x * xStep;
        y = y * yStep;
    end
    
    kam = kam_thresholded_deg(:); % Ensure column vector
    kam_beta = kam_beta_deg(:);  % ensure column vector
    %gos = gos(:); % Already filtered to ebsd_alpha points
    gos = zeros(size(kam));

    assert(numel(x) == numel(kam) && numel(kam) == numel(gos));
    
    % Now create the table
    T = table(x, y, kam, kam_beta, gos);
    
    % Extract folder and file name
    [folder, filename, ext] = fileparts(fname);
    
    % Create the file name
    outfile_name = filename + '_kam_gos_export.csv';
    
    % Combine with results directory into ONE path
    output_filename = fullfile(resultsdir, outfile_name);  % NO [ ] brackets
    
    % Ensure it's a character vector
    output_filename = char(output_filename);
    
    % Write to CSV
    writetable(T, output_filename);
    
    






    % %% GOS vs Y (pointwise, binning, and smoothing like KAM)
    % 
    % % Assign grain GOS value to each EBSD point
    % ebsd_alpha.prop.GOS = nan(size(ebsd_alpha));
    % 
    % for i = 1:numGrains
    %     idx = (ebsd_alpha.grainId == i);
    %     ebsd_alpha(idx).prop.GOS = grains.GOS(i);
    % end
    % 
    % % Extract Y positions and GOS values
    % y_vals = ebsd_alpha.y * yStep;  % convert to microns
    % GOS_vals = ebsd_alpha.prop.GOS;
    % 
    % % Parameters
    % bin_width_y = 1;              % Bin width in microns (adjust as needed)
    % min_points_per_bin = 10;      % Minimum EBSD points per bin to compute average
    % smooth_window = 5;            % Smoothing window size (points)
    % 
    % % Define bin edges and centers in microns
    % y_edges = min(y_vals):bin_width_y:max(y_vals);
    % y_centers = y_edges(1:end-1) + bin_width_y/2;
    % 
    % % Initialize binned averages with NaNs
    % binned_avg_gos_y = nan(size(y_centers));
    % 
    % % Calculate mean GOS within each bin
    % for i = 1:length(y_centers)
    %     idx = (y_vals >= y_edges(i)) & (y_vals < y_edges(i+1));
    %     if sum(idx) >= min_points_per_bin
    %         binned_avg_gos_y(i) = mean(GOS_vals(idx), 'omitnan');
    %     end
    % end
    % 
    % % Smooth binned averages ignoring NaNs
    % binned_avg_gos_y_smooth = nan(size(binned_avg_gos_y));
    % half_win = floor(smooth_window / 2);
    % 
    % for i = 1:length(binned_avg_gos_y)
    %     win_start = max(1, i - half_win);
    %     win_end = min(length(binned_avg_gos_y), i + half_win);
    %     window_vals = binned_avg_gos_y(win_start:win_end);
    %     if all(isnan(window_vals))
    %         binned_avg_gos_y_smooth(i) = NaN;
    %     else
    %         binned_avg_gos_y_smooth(i) = mean(window_vals, 'omitnan');
    %     end
    % end
    % 
    % % Plot results
    % figure;
    % plot(y_centers, binned_avg_gos_y, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Raw binned avg');
    % hold on;
    % plot(y_centers, binned_avg_gos_y_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed avg');
    % xlabel('Y position (\mum)', 'FontSize', 12);
    % ylabel('Average GOS (°)', 'FontSize', 12);
    % title('Average Grain Orientation Spread (GOS) vs Y Position', 'FontSize', 14);
    % grid on;
    % legend('Location', 'best');
    % set(gca, 'FontSize', 12);
    % hold off;
    % === EVERYTHING ELSE GOES HERE ===
    % Don't delete this block — just place all the rest of your code after it.
end
%%

% %% Select an individual grain for texture based segmentation
% % use grains(x,y) to select a grain from a map based on x and y
% % coordinates
% % grain_selected = grains(200,500);
% %3 check you have the right grain by plotting, then calculate the mean
% %orientation
% 
% %create a container for EBSD phase 1 (alpha - HCP)
% ebsd_2=ebsd(phase1);
% 
% figure('Position',PlotData.ssize);
% subplot(1,3,1);
% oM1.inversePoleFigureDirection = xvector;
% plot(grains(phase1),oM1.orientation2color(grains(phase1).meanOrientation),'parent',gca);
% title('Crystal Orientation - X')
% 
% 
% % take a mouse input
%  [x,y]=ginput(1);
% 
% % % hard code the mouse input to run this - swap the commenting around
% % % if you want a ginput
% % % hard coded position, in um
% % x=340; y=40;
% 
% %add the point to the figure as a black spot
% hold on;
% scatter(x,y,100,'k','filled');
% 
% % find the corresponding grain
% grain_sel = grains(x,y);
% 
% plot(grain_sel.boundary,'linecolor','g','LineWidth',3,'parent',gca)
% plot(grain_sel.boundary,'linecolor','r','LineWidth',1,'parent',gca)
% hold off
% 
% orientation_threshold=20; %in degrees
% % Segment EBSD Cubes based on grain orientation 
% %1 Define a fibre
% % fibre_test = fibre(Miller(1,1,0,ebsd(phase1).CS),yvector);
% %2 Extract Orientations using angle function to within 20 degrees of the fibre axis (generates a logical array) 
% test_orientations_2 = angle(ebsd_2.orientations,grain_sel.meanOrientation,'antipodal')<orientation_threshold*degree; % 20 degrees is arbitary
% % test_orientations_2 = angle(ebsd_2.orientations,fibre_test,'antipodal')<20*degree; 
% 
% %3 Segment original cube based on Phase (can only do this with one phase)
% test_orientations_in = ebsd_2(test_orientations_2 == 1);
% test_orientations_out = ebsd_2(test_orientations_2 == 0);
% 
% %5 Plot to check
% subplot(1,3,2);
% plot(test_orientations_in, oM1.orientation2color(test_orientations_in.orientations),'parent',gca);
% title(['Data close to selected grain (<' int2str(orientation_threshold) '^o)'])
% subplot(1,3,3);
% plot(test_orientations_out, oM1.orientation2color(test_orientations_out.orientations),'parent',gca);
% title(['Data close from selected grain (>' int2str(orientation_threshold) '^o)'])
% %% ODF for segmentation
% % calculate an ODF for the Alpha phase
% odf_in = calcODF(test_orientations_in.orientations,'halfwidth',5*degree);
% odf_out = calcODF(test_orientations_out.orientations,'halfwidth',5*degree);
% 
% % decide on the representations we want in the PDF
% h_plot = Miller({0,0,1},{1,-1,0},{1,-2,0},odf_in.CS);
% 
% pf_range=0.1:0.5:12;
% 
% FigH.ODFin=figure;
% plotPDF(odf_in,h_plot,'antipodal','silent','upper','projection','eangle','contourf',pf_range,'minmax');
% colorbar;
% 
% FigH.ODFout=figure;
% plotPDF(odf_out,h_plot,'antipodal','silent','upper','projection','eangle','contourf',pf_range,'minmax');
% colorbar;

% %% copy this m file over, for archival purposes
% 
% mf_long=mfilename('fullpath');
% [f1,f2,f3]=fileparts(mf_long);
% mf_start=[mf_long '.m'];
% mf_end=fullfile(resultsdir,[f2 '.m']);
% try
% copyfile(mf_start,mf_end)
% catch
%     warning('m file not saved, likely due to spaces in the file name');
% end
% 
% %eos