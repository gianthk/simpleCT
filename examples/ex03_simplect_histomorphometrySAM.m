%ex03_simplect_histomorphometrySAM.m
%   histomorphometric analysis from SAM images
%   
%   This script contains the SAM image processing described in the paper:
%       Iori G, Schneider J, Reisinger A, Heyer F, Peralta L, Wyers C, et al. Large cortical bone pores in the tibia are associated with proximal femur strength. PLOS ONE. doi:10.1371/journal.pone.0215405
%
%   This script is also available at:
%       https://doi.org/10.5281/zenodo.2605365
%
%   preliminary operations:
%   - SAM RF signals are already converted to maps of the surface's acoustic impedance
%   - All images are rotated in a common APML ref. system: ANT looks towards the bottom of the image
%   - Endosteum masks already drawn
%
%   SAM images, masks and master table for the analysis:
%         DOI: 10.5281/zenodo.2593855
%
%   dependencies:
%   - MATLAB R2018a or later
%   - simpleCT                              https://github.com/gianthk/simpleCT
%											DOI: 10.5281/zenodo.2628321
%   - MATLAB and Octave Functions for Computer Vision and Image Processing, Peter Kovesi
%                                           https://www.peterkovesi.com/matlabfns/
%   - export_fig (if graphics TRUE)         https://de.mathworks.com/matlabcentral/fileexchange/23629-export_fig
%   other dependencies that are included in the simpleCT toolbox
%   - dir2                                  https://de.mathworks.com/matlabcentral/fileexchange/21791-search-files-recursively-dir2
%   - fwhm                                  https://de.mathworks.com/matlabcentral/fileexchange/10590-fwhm
%   - max_inscribed_circle                  https://de.mathworks.com/matlabcentral/fileexchange/30805-maximum-inscribed-circle-using-distance-transform
%   - inpoly                                https://github.com/dengwirda/inpoly
%
%   This script does:
%     1.    load Scanning Acoustic Microscopy (SAM) image of cortical bone
%     2.    Convert to acoustic impedance
%     3.    Wellner's adaptive threshold
%     4.    clean the binary mask
%     5.    correct the periosteum (PERI), endosteum (ENDO) and cortical bone (CORT) masks
%     6.    load evaluation (ROI) masks for specific compartments
%     7.    measure macroscopic morph (shaft cross-section geometry)
%     8.    calculate Ct.Th distribution9. calculate cortical thickness (Ct.Th) distribution
%     10.1  calculate pore properties from binary image of the pores (Ct.Po, Po.D, Po.Dm)
%     10.2  calculate pore properties for each eval ROI
%     11.   store single sample results
%     12-13.save results
%   ______________________________________________________
%
%   Author:         Gianluca Iori <gianthk.iori@gmail.com>
%   BSRT - Charite Berlin
%   Created on:   25/04/2018
%   Last update:  25/03/2019
%   ______________________________________________________
%% init
clear all
close all
clc

graphics = true;
currentdir = pwd;

% master table
filename_master = 'tacosound_tibia_left_master.csv';     % the master file is contained in the archive that you downloaded

% path definitions
homedir = currentdir;
resultsdir = [homedir filesep 'results' filesep];
evalmasksdir = [homedir filesep 'masks' filesep '00_eval'];
ROIUSmasksdir = [homedir filesep 'masks' filesep '00_eval' filesep 'ROIUS'];

% create dir struct for analysis
mkdir(resultsdir);
mkdir([homedir filesep 'masks' filesep]);
mkdir([homedir filesep 'masks' filesep 'tissue' filesep]);
mkdir([homedir filesep 'masks' filesep 'endo' filesep]);
mkdir([homedir filesep 'masks' filesep 'peri' filesep]);
mkdir([homedir filesep 'masks' filesep 'cort' filesep]);
mkdir([homedir filesep 'masks' filesep 'pores' filesep]);
mkdir(evalmasksdir);
mkdir(ROIUSmasksdir);

warning('copy the files in the appropriate folders before running the analysis');

resultsdir = [resultsdir filesep datestr(now, 'yymmdd')];
resultsingledir = [resultsdir filesep 'single'];
if exist(resultsdir)==0, mkdir(resultsdir); end
if exist(resultsingledir)==0, mkdir(resultsingledir); end
resultsfile = [resultsdir filesep 'histo.mat'];

% SAM images for the analysis can be downloaded here: DOI
% you should copy the files in the appropriate folders before running the analysis:
%   - SAM images            -> homedir
%   - endosteum masks       -> homedir/masks/endo
%   - ROI masks             -> homedir/masks/00_eval    Each folder contains the specific ROIs for the analysis.
%                                                       ROIUS (anteromedial tibia) are provided.
%                                                       You can draw more ROI masks and put them here.. they will be processed)

% parameters
pixelsize = 0.012;                      % [mm] 100MHz SAM

ctth_BinEdges = 0.00:0.05:10.00;        % [mm]
ctth_bins = 0.025:0.05:9.975;           % [mm]

PoDm_BinEdges = 0.000:0.005:0.600;      % [mm]
PoDm_bins = 0.0025:0.005:0.5975;        % [mm]
PoDmprops = {'mean' 'median' 'std' 'min' 'max' 'Width' 'Mean' 'Peak' 'Variance' 'Skewness' 'Peakedness' '0.05 percentile' '0.95 percentile'};
PoDm_thresh = [0, 0.060 0.090 0.120];   % [mm]
allporesDm = [];
allporesDmUS = [];

PoAprops = {'mean' 'median' 'std' 'min' 'max' 'Width' 'Peak' 'Variance' 'Skewness' 'Kurtosis' '0.10 quantile' '0.90 quantile'};
PoA_thresh = [0 0.005 0.008 0.012 0.016 0.032 0.05];        % [mm^2]

PoDegCirprops = {'mean' 'median' 'std' 'min' 'max'};

sample = 1;
%% load master
if isempty(filename_master) || ~exist(filename_master)
    [FILENAME, PATHNAME] = uigetfile('*.csv;*.CSV', 'Select master table');
    cd(currentdir);
    filename_master = [PATHNAME FILENAME];
end

master = medtooltable.readtable(filename_master);
cd(currentdir);
%% eval all samples
for sample = 1:length(master.specimen)
    close all
    %%  checks
    samplename = master.specimen{sample};
    samplename = regexprep(samplename, '_', '');

    FILENAME        = [homedir filesep samplename '.tif'];
    TISSUEfilename  = [homedir filesep 'masks' filesep 'tissue' filesep samplename '.tif'];
    ENDOfilename    = [homedir filesep 'masks' filesep 'endo' filesep samplename '.tif'];
    PERIfilename    = [homedir filesep 'masks' filesep 'peri' filesep samplename '.tif'];
    CORTfilename    = [homedir filesep 'masks' filesep 'cort' filesep samplename '.tif'];
    
    isright = ~isempty(regexp(samplename, 'R'));    % LEFT (L) or RIGHT (R) flag
    
    % run flag in mastertable
    if ~master.run(sample), continue; end
    
    fprintf('evaluate specimen: \t%s\n', char(master.specimen(sample)));
    %%     1.   load Scanning Acoustic Microscopy (SAM) image of cortical bone
    % already rotated image in the APSI (Anterior-Posterior-Superior-Inferior) system
    SAM = imread(FILENAME);
    %%     2.   convert to acoustic impedance
    IMP = double(SAM)/4369 ;
    %%     3.   Wellner's adaptive threshold
    [BW] = threshbone(IMP, 'adaptive', 3, 22, 'gaussian', 'relative');    % adaptive threshold
    %%     4.   clean the binary mask
    % remove non connected particles
    BW2 = removeparticles(BW, 'connectivity');
    % remove single pixel pores
    BW2 = bwmorph(BW2, 'fill');
    % ..show before and after
    if graphics
        figure;
        imshowpair(BW2, BW); axis image; title('tissue mask cleaned');
    end
    %%     5.   correct the masks PERI, ENDO and CORT
    % Read the endosteum mask which was drawn manually using the set of rules proposed in (Malo et. al. 2013)
    ENDOfilename_in = [homedir filesep 'masks' filesep 'endo' filesep samplename '.tif'];
    endo_BW = imread(ENDOfilename_in);  endo_BW0 = endo_BW;
    % get boundary of the endosteum mask
    endo_bound_BW = bwmorph(~endo_BW, 'remove');
    
    % calculate the periosteum mask
    peri_BW = imfill(BW2 | endo_bound_BW, 'holes');
    % get boundary of the periosteum mask
    peri_bound_BW = bwmorph(peri_BW, 'remove');

    % final version of the endosteum mask from a combination with the tissue mask
    endo_BW = endosteummask(endo_BW0 & BW2 | peri_bound_BW, 'regionprops');

    % cortical bone mask as the combination of the two
    cort_BW = endo_BW & peri_BW;
    % remove non connected objects with an area < 300 pixels from the binary image
    cort_BW = bwareaopen(cort_BW, 300);

    % save the masks
    imwrite(endo_BW , ENDOfilename);
    imwrite(peri_BW , PERIfilename);
    imwrite(cort_BW , CORTfilename);

    if graphics
        figure; endo = gcf;
        imshowpair(IMP, endo_BW); axis image; title('endo');
        figure; peri = gcf;
        imshowpair(IMP, peri_BW); axis image; title('peri');
        figure; cort = gcf;
        imshowpair(IMP, cort_BW); axis image; title('cort');
    end
    %%     6.   load evaluation (ROI) masks
    fprintf('Loading ROI masks...\n');
    evalmasks_content = dir2(evalmasksdir);
    for j=1:length(evalmasks_content)
        evalmasks(j).name = evalmasks_content(j).name;
        evalmasks(j).filename = [evalmasksdir filesep evalmasks_content(j).name filesep samplename '.tif'];
        evalmasks(j).BW = logical(imread(evalmasks(j).filename));
        evalmasks(j).cort_BW = evalmasks(j).BW & cort_BW;
        evalmasks(j).polygon = periosteumcontour(flipud(evalmasks(j).BW));
        evalmasks(j).polygon = periosteumcontour(evalmasks(j).BW);
    end
    %%     7.   macroscopic morph (shaft cross-section geometry): (BA, TA, CA, BV/TV, Ixx, Iyy)
    tic
    
    % 6.2 Tissue Area (only the bone tissue) (T.Ar)
    TAr = sum(BW2(:));                          % [pixels]
    % 6.3 Cortical Area (area of the cortical compartment) (Ct.Ar)
    CtAr = sum(cort_BW(:));                     % [pixels]
    % 6.4 Core Area (the cross-sectional area occupied by the bone) (C.Ar)
    CAr = sum(peri_BW(:));                      % [pixels]
    % 6.5 areal portion of cortex (Ct.Wba)
    % Ratio between the cortical tissue area and the whole bone area (containing both the trabecular and cortical tissue). [Malo et. al. 2013]
    CtWba = sum(cort_BW(:) & BW2(:)) / CAr;     % [%]
    % 6.6 moment of inertia
    MoI = momentofinertia(BW2);                 % [pixels^2]
    % 6.7 second moments of area
    [Ixx, Iyy] = momentofinertia_area(BW2);     % [pixels^4]
    
    % convert to mm
    TAr = TAr * pixelsize^2;                    % [mm^2]
    CtAr = CtAr * pixelsize^2;                  % [mm^2]
    CAr = CAr * pixelsize^2;                    % [mm^2]
    MoI = MoI * pixelsize^2;                    % [mm^2]    
    Ixx = Ixx * pixelsize^4;                    % [mm^4]
    Iyy = Iyy * pixelsize^4;                    % [mm^4]
    
    fprintf('calculation of macroscopic morphological properties done!\t');
    toc
    %%     8.   Ct.Th distribution
    tic
    % get point clouds on the endosteum and periosteum contours
    peri_points = periosteumcontour(imclose(peri_BW, strel('disk', 30)));
    endo_points = endosteumcontour(endo_BW);

    % get only those belonging to a specific ROI 
    for j=1:length(evalmasks_content)
        evalmasks(j).endo_points = endo_points(inpolygon(endo_points(:,1), endo_points(:,2), evalmasks(j).polygon(:,1), evalmasks(j).polygon(:,2)),:);
        evalmasks(j).peri_points = peri_points(inpolygon(peri_points(:,1), peri_points(:,2), evalmasks(j).polygon(:,1), evalmasks(j).polygon(:,2)),:);
    end

    ctth_filename = [resultsdir filesep samplename '_CtTh'];    % output figure
    if graphics
        [ctth, ctth_array] = CtTh(peri_points, endo_points, 'pdist2', 1, 'peak');
        CtTh_fig = gcf;
        savefig(ctth_filename);
        export_fig(ctth_filename, '-dpng');                     % you will need export_fig for this.. you can modify this to 'print' if you don't have export_fig
        close(CtTh_fig);
    else
        [ctth, ctth_array] = CtTh(peri_points, endo_points, 'pdist2', 0, 'peak');
    end
    
    % convert Ct.Th to mm
    ctth = ctth * pixelsize;
    ctth_array = ctth_array * pixelsize;
    
    % histogram of the minimum Ct.Th
    ctth_hist   = histcounts(ctth_array, ctth_BinEdges, 'Normalization', 'pdf');
    
    % repeat for each eval ROI
    for j=1:length(evalmasks_content)
        if graphics
            [evalmasks(j).ctth, evalmasks(j).ctth_array] = CtTh(evalmasks(j).peri_points, evalmasks(j).endo_points, 'pdist2', 1, 'peak');
            CtTh_fig = gcf;
            ctth_filename = [evalmasks(j).filename(1:end-4) '_CtTh'];
            savefig(ctth_filename);
            export_fig(ctth_filename, '-dpng');
            close(CtTh_fig);
        else
            [evalmasks(j).ctth, evalmasks(j).ctth_array] = CtTh(evalmasks(j).peri_points, evalmasks(j).endo_points, 'pdist2', 0, 'peak');
        end
        evalmasks(j).ctth = evalmasks(j).ctth * pixelsize;
        evalmasks(j).ctth_array = evalmasks(j).ctth_array * pixelsize;
        evalmasks(j).ctth_hist  = histcounts(evalmasks(j).ctth_array, ctth_BinEdges, 'Normalization', 'pdf');
    end
    
    fprintf('calculation of cortical thickness done!\t');
    toc
    %%     10.1 pore props from BW of the pores (Ct.Po, Po.D, Po.Dm)
    tic
    
    % pores mask
    pores_BW = ~BW2 & cort_BW;
    % remove pores connected to the periosteum
    pores_BW = pores_BW & ~bwareaopen(pores_BW | ~peri_BW, 5000);
    % remove pores connected to the endosteum
    pores_BW = pores_BW & ~bwareaopen(pores_BW | ~endo_BW, 5000);
    
    % save it
    imwrite(pores_BW , [homedir filesep 'masks' filesep 'pores' filesep samplename '.tif']);
    
    % Ct.Po
    CtPo = 100*(sum(pores_BW(:))/sum(cort_BW(:)));  % [%]
    
    clear PoDmres PoDm_L PoDm_M PoDm_A PoDm_P PoDm_array PoA_array

    % Calculate pore properties
    PoSTATS = regionprops(pores_BW, 'Area');
    
    % Calculate Po.Dm distribution
    for  k=1:length(PoSTATS)
        PoA_array(k) = PoSTATS(k).Area;
    end
    PoDm_array = PoDm(pores_BW);
    PoDm_array = PoDm_array * pixelsize;        % [mm]
    allporesDm = [allporesDm; PoDm_array];
    PoA_array = PoA_array * pixelsize^2;        % [mm^2]
    
    % Degree of circularity
    PoDegCir_array = (pi*(PoDm_array.^2)/4)./PoA_array';
    
    % Po.Dm and Po.A histogram
    [PoDm_hist, PoDm_BinEdges] = histcounts(PoDm_array, 'Normalization', 'pdf');
    PoDm_bins = PoDm_BinEdges(1:end-1) + (PoDm_BinEdges(2)-PoDm_BinEdges(1))/2;
    PoDm_binwidth = PoDm_BinEdges(2:end) - PoDm_BinEdges(1:end-1);
    [PoA_hist, PoA_BinEdges] = histcounts(PoA_array, 'Normalization', 'pdf');
    PoA_bins = PoA_BinEdges(1:end-1) + (PoA_BinEdges(2)-PoA_BinEdges(1))/2;
    PoA_binwidth = PoA_BinEdges(2:end) - PoA_BinEdges(1:end-1);
    
    % relative proportion of pores with DIAMETER over a certain threshold and ammount of surface occupied per square mm
    DmAtot = PoDm_binwidth*PoDm_hist';                          % total area under the curve
    for t = 1:length(PoDm_thresh)
        idx = find(PoDm_array>PoDm_thresh(t));
        idx_hist = find(PoDm_bins>PoDm_thresh(t));
        
        DmAt = PoDm_binwidth(idx_hist)*PoDm_hist(idx_hist)';    % area under the curve after threshold
        
        % relative proportion of pores with diameter over a certain threshold
        PoDm_rel(t) = DmAt / DmAtot;
        
        % ammount of surface occupied
        PoDmA_rel(t) = sum(PoA_array(idx));
    end
    PoDmA_rel = PoDmA_rel / CtAr;                   % 1/BVTV [0-1]
    PoDmCtPo_rel = PoDmA_rel*100;                   % [%]
    PoDmA_rel = (PoDmA_rel/PoDmA_rel(1))*100;       % [%]
    
    % relative proportion of pores with AREA over a certain threshold and ammount of surface occupied
    AAtot = PoA_binwidth*PoA_hist';              % total area under the curve
    for t = 1:length(PoA_thresh)
        idx = find(PoA_array>PoA_thresh(t));
        idx_hist = find(PoA_bins>PoA_thresh(t));
        
        AAt = PoA_binwidth(idx_hist)*PoA_hist(idx_hist)';  % area under the curve after threshold
        
        % relative proportion of pores with diameter over a certain threshold
        PoA_rel(t) = AAt / AAtot;
        
        % ammount of surface occupied
        PoAA_rel(t) = sum(PoA_array(idx));
    end
    PoAA_rel = PoAA_rel / CtAr;                     % 1/BVTV [0-1]
    PoACtPo_rel = PoAA_rel*100;                     % [%]
    PoAA_rel = (PoAA_rel/PoAA_rel(1))*100;          % [%]

    % ksdensity of the Po.Dm distribution
    [prob, idx] = ksdensity(PoDm_array);
    % Po.Dm distr properties
    PoDmres     = [mean(PoDm_array) median(PoDm_array) std(PoDm_array) min(PoDm_array) max(PoDm_array) fwhm(PoDm_bins, PoDm_hist) idx(find(prob == max(prob))) var(PoDm_array) skewness(PoDm_array) kurtosis(PoDm_array) quantile(PoDm_array, 0.1) quantile(PoDm_array, 0.9)];
    
    % ksdensity of the Po.A distribution
    [prob, idx] = ksdensity(PoA_array);
    % Po.A distr properties
    PoA         = [mean(PoA_array) median(PoA_array) std(PoA_array) min(PoA_array) max(PoA_array) fwhm(PoA_bins, PoA_hist) idx(find(prob == max(prob))) var(PoA_array) skewness(PoA_array) kurtosis(PoA_array) quantile(PoA_array, 0.1) quantile(PoA_array, 0.9)];
    
    % Pore degree of circularity results
    PoDegCir = [mean(PoDegCir_array) median(PoDegCir_array) std(PoDegCir_array) min(PoDegCir_array) max(PoDegCir_array)];
        
    if graphics
        % Save figure to results
        PoDmhist_filename = [resultsdir filesep samplename '_PoDmHist'];
        figure;
        histogram(PoDm_array*1000, PoDm_BinEdges*1000, 'Normalization', 'pdf');
        hold on;
        % ksdensity
        plot(PoDm_bins*1000, ksdensity(1000*PoDm_array(:), 1000*PoDm_bins), '-.r','LineWidth',1.5)
        xlabel('Po.Dm [microns]')
        ylabel('Probability');
        title('TOT');
        % axis([0 300 0 1])
        axis 'auto y'
        grid MINOR
        legend('hist', 'ksdensity');
        PoDmhist_fig = gcf;
        savefig(PoDmhist_filename);
        export_fig(PoDmhist_filename, '-dpng');
        close(PoDmhist_fig);
        
        PoAhist_filename = [resultsdir filesep samplename '_PoAHist'];
        figure;
        histogram(PoA_array, PoA_BinEdges, 'Normalization', 'pdf');
        hold on;
        % ksdensity
        plot(PoA_bins, ksdensity(PoA_array(:), PoA_bins), '-.r','LineWidth',1.5)
        xlabel('Po.A [mm^2]')
        ylabel('Probability');
        title('TOT');
        % axis([0 0.04 0 1])
        axis 'auto y'
        grid MINOR
        legend('hist', 'ksdensity');
        PoAhist_fig = gcf;
        savefig(PoAhist_filename);
        export_fig(PoAhist_filename, '-dpng');
        close(PoAhist_fig);
    end
    
    % Po.D of pores with diameter > threshold
    for t = 1:length(PoDm_thresh)
        PoD(t) = sum(PoDm_array>PoDm_thresh(t)) / CtAr;     % [1/mm^2]
        Pon(t) = sum(PoDm_array>PoDm_thresh(t));
    end
    
    % Po.D of pores with diameter > threshold
    for t2 = 1:length(PoA_thresh)
        PoD(t+t2) = sum(PoA_array>PoA_thresh(t2)) / CtAr;   % [1/mm^2]
        Pon(t+t2) = sum(PoA_array>PoA_thresh(t2));
    end
    %%     10.2 pore props for each evalmask
    for j=1:length(evalmasks_content)
        CtAr_evalmask = sum(evalmasks(j).cort_BW(:)) * pixelsize^2;       % [mm^2]
        
        % pores mask
%         evalmasks(j).pores_BW = ~BW2 & imerode(evalmasks(j).cort_BW, strel('disk', 4));
        evalmasks(j).pores_BW = pores_BW & evalmasks(j).cort_BW;
        imwrite(evalmasks(j).pores_BW , [evalmasks(j).filename(1:end-4) '_pores.tif']);
        
        % Ct.Po
        evalmasks(j).CtPo = 100*(sum(evalmasks(j).pores_BW(:))/sum(evalmasks(j).cort_BW(:)));
        
        % Calculate Po.Dm distribution
        PoSTATS = regionprops(evalmasks(j).pores_BW, 'Area');
        evalmasks(j).PoA_array = [];
        for  k=1:length(PoSTATS)
            evalmasks(j).PoA_array(k) = PoSTATS(k).Area;
        end
        evalmasks(j).PoDm_array = PoDm(evalmasks(j).pores_BW);
        evalmasks(j).PoDm_array = evalmasks(j).PoDm_array * pixelsize;
        allporesDmUS = [allporesDmUS; evalmasks(j).PoDm_array];
        evalmasks(j).PoA_array = evalmasks(j).PoA_array * pixelsize^2;

        % Degree of circularity
        evalmasks(j).PoDegCir_array = (pi*(evalmasks(j).PoDm_array.^2)/4)./(evalmasks(j).PoA_array');
        
        % Calculate Po.Dm histogram
        evalmasks(j).PoDm_hist = histcounts(evalmasks(j).PoDm_array, PoDm_BinEdges, 'Normalization', 'pdf');
        evalmasks(j).PoA_hist = histcounts(evalmasks(j).PoA_array, PoA_BinEdges, 'Normalization', 'pdf');
        
        % relative proportion of pores with diameter over a certain threshold
        % and ammount of surface occupied per square mm
        DmAtot = PoDm_binwidth*evalmasks(j).PoDm_hist';          	% total area under the curve
        for t = 1:length(PoDm_thresh)
            idx = find(evalmasks(j).PoDm_array>PoDm_thresh(t));
            idx_hist = find(PoDm_bins>PoDm_thresh(t));
            
            DmAt = PoDm_binwidth(idx_hist)*evalmasks(j).PoDm_hist(idx_hist)';	% area under the curve after threshold
            
            % relative proportion of pores with diameter over a certain threshold
            evalmasks(j).PoDm_rel(t) = DmAt / DmAtot;
            
            % ammount of surface occupied per square mm
            evalmasks(j).PoDmA_rel(t) = sum(evalmasks(j).PoA_array(idx));
        end
        evalmasks(j).PoDmA_rel = evalmasks(j).PoDmA_rel / CtAr_evalmask;                        % [0-1]
        evalmasks(j).PoDmCtPo_rel = evalmasks(j).PoDmA_rel*100;                                 % [%]
        evalmasks(j).PoDmA_rel = (evalmasks(j).PoDmA_rel/evalmasks(j).PoDmA_rel(1))*100;        % [%]

        % relative proportion of pores with Area over a certain threshold
        % and ammount of surface occupied per square mm
        AAtot = PoA_binwidth*evalmasks(j).PoA_hist';          	% total area under the curve
        for t = 1:length(PoA_thresh)
            idx = find(evalmasks(j).PoA_array>PoA_thresh(t));
            idx_hist = find(PoA_bins>PoA_thresh(t));
            
            AAt = PoA_binwidth(idx_hist)*evalmasks(j).PoA_hist(idx_hist)';	% area under the curve after threshold
            
            % relative proportion of pores with diameter over a certain threshold
            evalmasks(j).PoA_rel(t) = AAt / AAtot;
            
            % ammount of surface occupied per square mm
            evalmasks(j).PoAA_rel(t) = sum(evalmasks(j).PoA_array(idx));
        end
        evalmasks(j).PoAA_rel = evalmasks(j).PoAA_rel / CtAr_evalmask;                          % [0-1]
        evalmasks(j).PoACtPo_rel = evalmasks(j).PoAA_rel*100;                                   % [%]
        evalmasks(j).PoAA_rel = (evalmasks(j).PoAA_rel/evalmasks(j).PoAA_rel(1))*100;           % [%]
        
        % ksdensity of the Po.Dm distribution
        [prob, idx] = ksdensity(evalmasks(j).PoDm_array);
        % Po.Dm distr properties
        evalmasks(j).PoDm   = [mean(evalmasks(j).PoDm_array) median(evalmasks(j).PoDm_array) std(evalmasks(j).PoDm_array) min(evalmasks(j).PoDm_array) max(evalmasks(j).PoDm_array) fwhm(PoDm_bins, evalmasks(j).PoDm_hist) idx(find(prob == max(prob))) var(evalmasks(j).PoDm_array) skewness(evalmasks(j).PoDm_array) kurtosis(evalmasks(j).PoDm_array) quantile(evalmasks(j).PoDm_array, 0.1) quantile(evalmasks(j).PoDm_array, 0.9)];
        
        % ksdensity of the Po.Dm distribution
        [prob, idx] = ksdensity(evalmasks(j).PoA_array);
        % Po.A distr properties
        evalmasks(j).PoA    = [mean(evalmasks(j).PoA_array) median(evalmasks(j).PoA_array) std(evalmasks(j).PoA_array) min(evalmasks(j).PoA_array) max(evalmasks(j).PoA_array) fwhm(PoA_bins, evalmasks(j).PoA_hist) idx(find(prob == max(prob))) var(evalmasks(j).PoA_array) skewness(evalmasks(j).PoA_array) kurtosis(evalmasks(j).PoA_array) quantile(evalmasks(j).PoA_array, 0.1) quantile(evalmasks(j).PoA_array, 0.9)];
        
        % Pore degree of circularity results
        evalmasks(j).PoDegCir = [mean(evalmasks(j).PoDegCir_array) median(evalmasks(j).PoDegCir_array) std(evalmasks(j).PoDegCir_array) min(evalmasks(j).PoDegCir_array) max(evalmasks(j).PoDegCir_array)];
        
        if graphics
            % Save figure to results
            PoDmhist_filename = [resultsdir filesep samplename '_PoDmHist_' evalmasks(j).name];
            figure;
            histogram(evalmasks(j).PoDm_array*1000, PoDm_BinEdges*1000, 'Normalization', 'pdf');
            hold on;
            % ksdensity
            plot(PoDm_bins*1000, ksdensity(1000*evalmasks(j).PoDm_array(:), 1000*PoDm_bins), '-.r','LineWidth',1.5)
            xlabel('Po.Dm [microns]')
            ylabel('Probability');
            title(evalmasks(j).name);
            % axis([0 300 0 1])
            axis 'auto y'
            grid MINOR
            legend('hist', 'ksdensity');
            PoDmhist_fig = gcf;
            savefig(PoDmhist_filename);
            export_fig(PoDmhist_filename, '-dpng');
            close(PoDmhist_fig);
            
            PoAhist_filename = [resultsdir filesep samplename '_PoAHist_' evalmasks(j).name];
            figure;
            histogram(evalmasks(j).PoA_array, PoA_BinEdges, 'Normalization', 'pdf');
            hold on;
            % ksdensity
            plot(PoA_bins, ksdensity(evalmasks(j).PoA_array(:), PoA_bins), '-.r','LineWidth',1.5)
            xlabel('Po.A [mm^2]')
            ylabel('Probability');
            title(evalmasks(j).name);
            % axis([0 0.04 0 1])
            axis 'auto y'
            grid MINOR
            legend('hist', 'ksdensity');
            PoAhist_fig = gcf;
            savefig(PoAhist_filename);
            export_fig(PoAhist_filename, '-dpng');
            close(PoAhist_fig);
        end
        
        % Po.D of pores with diameter > threshold
        for t = 1:length(PoDm_thresh)
            evalmasks(j).PoD(t) = sum(evalmasks(j).PoDm_array>PoDm_thresh(t)) / CtAr_evalmask;       % [1/mm^2]
            evalmasks(j).Pon(t) = sum(evalmasks(j).PoDm_array>PoDm_thresh(t));
        end
        
        % Po.D of pores with Area > threshold
        for t2 = 1:length(PoA_thresh)
            evalmasks(j).PoD(t+t2) = sum(evalmasks(j).PoA_array>PoA_thresh(t2)) / CtAr_evalmask;       % [1/mm^2]
            evalmasks(j).Pon(t+t2) = sum(evalmasks(j).PoA_array>PoA_thresh(t2));
        end
        
    end
    %%     11.  save single sample results
    resfile                 = [resultsingledir filesep samplename '.mat'];
    res.script              = mfilename;
    res.sample              = samplename;
    res.evalmasks           = evalmasks;
    res.pixelsize           = pixelsize;
    res.CtTh_BinEdges       = ctth_BinEdges;
    res.CtTh_bins           = ctth_bins;
    
    res.PoDm_array          = PoDm_array;
    res.PoDm_BinEdges       = PoDm_BinEdges;
    res.PoDm_bins           = PoDm_bins;
    res.PoDmprops           = PoDmprops;
    res.PoDm_thresh         = PoDm_thresh;
    res.PoDm                = PoDmres;
    res.PoDm_hist           = PoDm_hist;
    res.PoDm_rel            = PoDm_rel;
    res.PoDmCtPo_rel        = PoDmCtPo_rel;
    res.PoDmA_rel           = PoDmA_rel;
    
    res.PoA_array           = PoA_array;
    res.PoA_BinEdges        = PoA_BinEdges;
    res.PoA_bins            = PoA_bins;
    res.PoAprops            = PoAprops;
    res.PoA_thresh          = PoA_thresh;
    res.PoA                 = PoA;
    res.PoA_hist            = PoA_hist;
    res.PoA_rel             = PoA_rel;
    res.PoACtPo_rel         = PoACtPo_rel;
    res.PoAA_rel            = PoAA_rel;
    
    res.TAr                 = TAr;
    res.CtAr                = CtAr;
    res.CAr                 = CAr;
    res.CtWba               = CtWba;
    res.MoI                 = MoI;
    res.Ixx                 = Ixx;
    res.Iyy                 = Iyy;
    res.peri_points         = peri_points;
    res.endo_points         = endo_points;
    res.CtTh                = ctth;
    res.CtTh_array          = ctth_array;
    res.CtPo                = CtPo;
    res.PoDegCirprops       = PoDegCirprops;
    res.PoDegCir            = PoDegCir;
    res.PoD                 = PoD;
    res.Pon                 = Pon;
    
    save(resfile, 'res', '-v7.3');
    %%     12. store results for all samples
    histo.script                        = mfilename;
    histo.samplename{sample}            = samplename;
    histo.pixelsize                     = pixelsize;
    
    histo.PoDmprops                     = PoDmprops;
    histo.PoDm_thresh                   = PoDm_thresh;
    histo.PoDm_CORT(sample,:)           = PoDmres;
    histo.PoDm_rel_CORT(sample,:)       = PoDm_rel;
    histo.allporesDm_CORT               = allporesDm;
    histo.PoDmCtPo_rel_CORT(sample,:)   = PoDmCtPo_rel;
    histo.PoDmA_rel_CORT(sample,:)      = PoDmA_rel;
    histo.PoDmCtPo_rel_CORT(sample,:)   = PoDmCtPo_rel;
    histo.PoDmA_rel_CORT(sample,:)      = PoDmA_rel;
    
    histo.PoAprops                      = PoAprops;
    histo.PoA_thresh                    = PoA_thresh;
    histo.PoA_rel_CORT(sample,:)        = PoA_rel;
    histo.PoACtPo_rel_CORT(sample,:)    = PoACtPo_rel;
    histo.PoAA_rel_CORT(sample,:)       = PoAA_rel;
    histo.PoA_CORT(sample,:)            = PoA;
    
    histo.TAr(sample)                   = TAr;
    histo.CtAr(sample)                  = CtAr;
    histo.CAr(sample)                   = CAr;
    histo.CtWba(sample)                 = CtWba;
    histo.MoI(sample)                   = MoI;
    histo.Ixx(sample)                   = Ixx;
    histo.Iyy(sample)                   = Iyy;
    
    histo.CtTh_BinEdges                 = ctth_BinEdges;
    histo.CtTh_bins                     = ctth_bins;
    histo.CtTh_CORT(sample)             = ctth;
    
    histo.CtPo_CORT(sample)             = CtPo;
    histo.PoD_CORT(sample,:)            = PoD;
    histo.Pon_CORT(sample,:)            = Pon;    
    histo.PoDegCir_CORT(sample,:)       = PoDegCir;
    histo.PoDegCirprops                 = PoDegCirprops;
    
    for j=1:length(evalmasks_content)
        histo.evalmasks{j} = evalmasks(j).name;
        histo.allporesDm{j} = allporesDmUS;
        histo.CtTh(sample,j) = evalmasks(j).ctth;
        histo.CtPo(sample,j) = evalmasks(j).CtPo;
        histo.PoDm(sample,:,j) = evalmasks(j).PoDm;
        histo.PoDm_rel(sample,:,j) = evalmasks(j).PoDm_rel;
        histo.PoA_rel(sample,:,j) = evalmasks(j).PoA_rel;
        histo.PoDmCtPo_rel(sample,:,j) = evalmasks(j).PoDmCtPo_rel;
        histo.PoDmA_rel(sample,:,j) = evalmasks(j).PoDmA_rel;
        histo.PoACtPo_rel(sample,:,j) = evalmasks(j).PoACtPo_rel;
        histo.PoAA_rel(sample,:,j) = evalmasks(j).PoAA_rel;
        histo.PoA(sample,:,j) = evalmasks(j).PoA;
        histo.PoA(sample,:,j) = evalmasks(j).PoA;
        histo.PoDegCir(sample,:,j) = evalmasks(j).PoDegCir;
        histo.PoD(sample,:,j) = evalmasks(j).PoD;
        histo.Pon(sample,:,j) = evalmasks(j).Pon;
    end
end
%%         13. save results
save(resultsfile, 'histo', '-v7.3');