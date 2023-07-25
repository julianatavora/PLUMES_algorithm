%
%-------------------------------------------------------------------------%
% PLUMES algorithm
%
% Howdy PLUMES user 
% 
% The following function is fully described in Tavora et al., 2023 (DOI
% 10.3389/fmars.2023.1215327). 
%
%-------------------------------------------------------------------------%
%
% you will need: 
% - a satellite scene masked* (im variable), 
% - lat/lon of satellite scene (lat and lon variables), 
% - date of scene (date variable), 
% - and chosen control points** (refer to Tavora et al., 2023 on how to choose control points)
%
% *land masked + and anything inward from the origin point of the
% plume masked as well, refer to examples available
%
% ** you will need to edit LAT/LON (rows 39 and 75)
%
% refer to Examples/ for application
%
%-------------------------------------------------------------------------%
% Juliana Tavora, 17/july/2023 (ITC, University of Twente)
%-------------------------------------------------------------------------%

function [PLUME_distal,PLUME_table_distal,PLUME_proximal,PLUME_table_proximal] = PLUMES_algorithm(date,im,lat,lon)

im_log = log10(im);

box_size = 4; 

%-------------------------------------------------------------------------%
%                        get samples from estuary                         %
%-------------------------------------------------------------------------%

LAT = -32.1917; LON = -52.0729;  %official coordinates

[seed_row, seed_col] = findClosestPixel(lon, lat, LON, LAT);

est_samples_coord = [seed_row,seed_col];

for i = 1:size(est_samples_coord,1)
    
    est_samples(1:box_size+1,(i-1)*(box_size+1)+1:i*(box_size+1)) = ...
        im(est_samples_coord(i,1)-box_size/2:est_samples_coord(i,1)+box_size/2,...
        est_samples_coord(i,2)-box_size/2:est_samples_coord(i,2)+box_size/2);
end

est_samples = reshape(est_samples,[],1);
est_samples_nan = sum(isnan(est_samples))/((box_size+1)*(box_size+1).*size(est_samples_coord,1));

if round(est_samples_nan,2) <= 0.35 %percentage of acceptable NaNs
    m_est = nanmedian(est_samples);        m_est_log = nanmedian(log10(est_samples));  %median of all samples for estuary
    v_est = std(est_samples,'omitnan');    v_est_log = nanvar(log10(est_samples));
else
    m_est = NaN;
    v_est = NaN;
    
    PLUME_distal   = im.*NaN;
    PLUME_proximal = im.*NaN;
    plume_contour_distal   = [NaN,NaN];
    plume_contour_proximal = [NaN,NaN];
    
    fprintf('Percentage of NaNs is above limit: %d', round(est_samples_nan,2));  % Does this work as expected?! %fprintf('Done\n')
    return
end

%-------------------------------------------------------------------------%
%                    get samples from marine water                        %
%-------------------------------------------------------------------------%

 LAT = [-32.1003 -32.3771 -32.3840];  LON = [-51.7958 -51.8076 -51.9521];

for ii=1:size(LAT,2)
    [row(ii), column(ii)] = findClosestPixel(lon, lat, LON(ii), LAT(ii));
end

mar_samples_coord = [row', column'];
for i = 1:size(mar_samples_coord,1)
    
    mar_samples(1:box_size+1,(i-1)*(box_size+1)+1:i*(box_size+1)) = ...
        im(mar_samples_coord(i,1)-box_size/2:mar_samples_coord(i,1)+box_size/2,...
        mar_samples_coord(i,2)-box_size/2:mar_samples_coord(i,2)+box_size/2);
end

mar_samples = reshape(mar_samples,[],1);
mar_samples_nan = sum(isnan(mar_samples))/((box_size+1)*(box_size+1).*size(mar_samples_coord,1));

if round(mar_samples_nan,2) <= 0.35 %percentage of acceptable NaNs
    m_mar = nanmedian(mar_samples);      m_mar_log = nanmedian(log10(mar_samples));
    v_mar = std(mar_samples,'omitnan');  v_mar_log = nanvar(log10(mar_samples));        %variance of all samples of marine waters
else
    m_mar = NaN;
    v_mar = NaN;
    
    PLUME_distal   = im.*NaN;
    PLUME_proximal = im.*NaN;
    plume_contour_distal   = [NaN,NaN];
    plume_contour_proximal = [NaN,NaN];
    
    fprintf('Percentage of NaNs is above limit: %d', round(mar_samples_nan,2));  % Does this work as expected?! %fprintf('Done\n')
    return
end

if m_mar >= m_est
    
    PLUME_distal   = im.*NaN;
    PLUME_proximal = im.*NaN;
    plume_contour_distal   = [NaN,NaN];
    plume_contour_proximal = [NaN,NaN];
    
    fprintf('No plume detected');  % Does this work as expected?! %fprintf('Done\n')
    return
    
end
%-------------------------------------------------------------------------%
%                              phase II                                   %
%-------------------------------------------------------------------------%



if exist('PLUME_distal','var') ==0
    
    
    %-------------------------------------------------------------------------%
    %                        start segmentation (distal)                      %
    %-------------------------------------------------------------------------%
    [lin, col] = size(im);
    scene = zeros(lin,col);
    
    for i = 1:lin
        for j = 1:col
            
            if ~isnan(im(i,j)) == 1
                
                S_est  = ((im(i,j) - (m_est))^2)/(v_est); %measure of similarity of pixel (i,j) with estuarine water
                S_mar  = ((im(i,j) - (m_mar))^2)/(v_mar); %measure of similarity of pixel (i,j) with marine water
                
                if m_est > m_mar
                    
                    if (S_est < S_mar)      %minimun difference %originally was <
                        scene(i,j) = 1;     %estuarine samples
                    else
                        scene(i,j) = 0;     %marine samples
                    end
                    
                else
                    scene(i,j) = 0;     %marine samples
                end
                
            else
                scene(i,j) = 0;       %land or cloud
                
            end
        end
    end
    
    %-------------------------------------------------------------------------%
    %                  start segmentation (proximal = core)                   %
    %-------------------------------------------------------------------------%
    scene_log = zeros(lin,col);
    
    
    for i = 1:lin
        for j = 1:col
            
            if ~isnan(im_log(i,j)) == 1
                
                S_est_log  = ((im_log(i,j) - (m_est_log))^2)/(v_est_log); %measure of similarity of pixel (i,j) with estuarine water
                S_mar_log  = ((im_log(i,j) - (m_mar_log))^2)/(v_mar_log); %measure of similarity of pixel (i,j) with marine water
                
                if m_est_log > m_mar_log
                    
                    if (S_est_log < S_mar_log)  %minimun difference %originally was <
                        scene_log(i,j) = 1;     %estuarine samples
                    else
                        scene_log(i,j) = 0;     %marine samples
                    end
                    
                else
                    scene_log(i,j) = 0;     %marine samples
                end
                
            else
                scene_log(i,j) = 0;       %land or cloud
                
            end
        end
    end
    
    
    %-------------------------------------------------------------------------%
    %                        morphological operations                         %
    %-------------------------------------------------------------------------%
    pixel_center  = scene(seed_row,seed_col);
    
    if pixel_center ~= 0
        
        %-------------------------outer limit-----------------------------%
        PLUME_distal   = reggrow(scene,seed_row,seed_col);
        
        %smoothing segmented boundary
        se = strel('disk',3); PLUME_distal = imclose(PLUME_distal,se);
        
        %Fill holes.
        PLUME_distal = imfill(PLUME_distal, 'holes');
        
        % get boundaries
        boundaries = bwboundaries(PLUME_distal); % Get list of (x,y) coordinates of outer perimeter.
        for k = 1:size(boundaries,1)
            h = roipoly(PLUME_distal,boundaries{k,1}(:,2),boundaries{k,1}(:,1));
            index(k) = h(seed_row,seed_col);
            
            if index(k) ==1
                
                %get LON/LAT of boundary
                plume_contour_distal = bound2coord(boundaries{k,1}, lon, lat);
                
                PLUME_distal = roipoly(PLUME_distal,boundaries{k,1}(:,2),boundaries{k,1}(:,1));
                
                %smoothing segmented boundary
                se = strel('disk',3); PLUME_distal = imclose(PLUME_distal,se);
                break
            end
            
        end
 
        %-------------------------------------------------------------------------------------------------------%
        %                   table of stats (distal plume)
        %-------------------------------------------------------------------------------------------------------%
        
        se          = strel('octagon',3);
        BW2         = imdilate(PLUME_distal,se);
        stats_plume = regionprops('table',BW2,'Centroid','Area','MajorAxisLength','MinorAxisLength','Orientation');
        SPM         = im.*PLUME_distal; SPM(SPM == 0) = NaN;
        stats_SPM   = array2table([nanmin(SPM(:)), nanmax(SPM(:)), nanmean(SPM(:)),nanmedian(SPM(:)),std(SPM(:),'omitnan'),m_est,v_est,m_mar,v_mar],...
            'VariableNames',{'Min SPM plume','Max SPM plume','Mean','Median SPM plume','Stdev SPM plume','control point origin (mean SPM)','control point origin (stdev SPM)',...
            'control point marine (mean SPM)','control point marine (stdev SPM)'});
        date        = table(datetime(date,'Format','dd.MM.yyyy'),'VariableNames',{'dd.mm.yyyy'}); 
        PLUME_table_distal = [date, stats_plume, stats_SPM];
        
        clear BW2 stats_plume stats_SPM  SPM
        %-------------------------------------------------------------------------------------------------------%

        if sum(index) <1
            
            PLUME_distal   = scene.*NaN;
            PLUME_proximal = scene.*NaN;
            plume_contour_distal   = [NaN,NaN];
            plume_contour_proximal = [NaN,NaN];
            
            fprintf('No plume detected!');  % Does this work as expected?! %fprintf('Done\n')
            
            return
        end
        
        
        clear boundaries
        
        %-------------------------inner limit-----------------------------%
        % the core of the plume (proximal) must be whithin the distal plume limits.
        % we first check if the plume core is located in the origin of the
        % of the plume. otherwise, we search for the location of highest
        % turbidity values whithin the segmented proximal plume
        
        scene_log =  PLUME_distal.*scene_log;
        
        pixel_center  = scene_log(seed_row,seed_col);
        
        if pixel_center == 0
            
            %find highest turbidity whithin distal plume limits
            B = ones(4,4)/4^2;
            C = conv2(scene_log.*im_log,B,'same');
            [~,c] = (max(C,[],'all','linear'));
            [seed_row,seed_col] = ind2sub(size(C),c); %clear C c B
        end
        
        PLUME_proximal = reggrow(scene_log,seed_row,seed_col);
        
        %smoothing segmented boundary
        se = strel('disk',3); PLUME_proximal = imclose(PLUME_proximal,se);
        
        PLUME_proximal = imfill(PLUME_proximal, 'holes'); % Fill holes.
        
        % get boundaries
        boundaries = bwboundaries(PLUME_proximal); % Get list of (x,y) coordinates of outer perimeter.
        for k = 1:size(boundaries,1)
            h = roipoly(PLUME_distal,boundaries{k,1}(:,2),boundaries{k,1}(:,1));
            index(k) = h(seed_row,seed_col);
            
            if index(k) ==1
                %get LON/LAT of boundary
                plume_contour_proximal = bound2coord(boundaries{k,1}, lon, lat);
                
                PLUME_proximal = roipoly(PLUME_proximal,boundaries{k,1}(:,2),boundaries{k,1}(:,1));
                
                 %smoothing segmented boundary
                 se = strel('disk',3); PLUME_proximal = imclose(PLUME_proximal,se);
                break
            end
            
        end
        
        %-------------------------------------------------------------------------------------------------------%
        %                   table of stats (proximal plume)
        %-------------------------------------------------------------------------------------------------------%
        se          = strel('octagon',3);
        BW2         = imdilate(PLUME_proximal,se);
        stats_plume = regionprops('table',BW2,'Centroid','Area','MajorAxisLength','MinorAxisLength','Orientation');
        SPM         = im_log.*PLUME_proximal; SPM(SPM == 0) = NaN; SPM = 10.^SPM;
        stats_SPM   = array2table([nanmin(SPM(:)), nanmax(SPM(:)), nanmean(SPM(:)),nanmedian(SPM(:)),std(SPM(:),'omitnan'),10.^m_est_log,10.^v_est_log,10.^m_mar_log,10.^v_mar_log],...
            'VariableNames',{'Min SPM plume','Max SPM plume','Mean','Median SPM plume','Stdev SPM plume','control point origin (mean SPM)','control point origin (stdev SPM)',...
            'control point marine (mean SPM)','control point marine (stdev SPM)'});
        
        PLUME_table_proximal =  [date, stats_plume, stats_SPM];
        
        clear BW2 stats_plume stats_SPM  SPM
        %-------------------------------------------------------------------------------------------------------%
        
        if sum(index) <1
            
            PLUME_distal   = scene.*NaN;
            PLUME_proximal = scene.*NaN;
            plume_contour_distal   = [NaN,NaN];
            plume_contour_proximal = [NaN,NaN];
            
            fprintf('No plume detected!');  
            
        end
        
        
        %-----------------------------------------------------------------%
    else
        
        PLUME_distal   = scene.*NaN;
        PLUME_proximal = scene.*NaN;
        plume_contour_distal   = [NaN,NaN];
        plume_contour_proximal = [NaN,NaN];
        
        fprintf('No plume detected!'); 
        
        
    end
    
else
    
    PLUME_distal   = scene.*NaN;
    PLUME_proximal = scene.*NaN;
    plume_contour_distal   = [NaN,NaN];
    plume_contour_proximal = [NaN,NaN];
    
    fprintf('No plume detected!'); 
    
    
end

end

function PLUME = reggrow(I,x,y)

% This function performs "region growing" in asegmented image from a specified
% seed-pixel (x,y)
%
% PLUME = reggrow(I,x,y) 
% 
% I   : input segmented image 
% PLUME : logical output image of region (PLUME estimated region)
% x,y : the position of the seed-pixel

% Based on function by D. Kroon, University of Twente
% (%https://nl.mathworks.com/matlabcentral/fileexchange/19084-region-growing)

% Adapted by Juliana Tavora, University of Twente



Isizes      = size(I); %size of the image
PLUME         = zeros(Isizes); %output

reg_mean    = I(x,y); % The mean of the segmented region (inicializado com valor do pixel semente)
reg_size    = 1; %number of pixels in region

% Free memory to store neighbours of the (segmented) region
neg_free    = 10000;  
neg_pos     = 0;
neg_list    = zeros(neg_free,2);

pixdist     = 1; % Distance of the region newest pixel to the regio mean

neighbor    = [-1 0; 1 0; 0 -1;0 1];  % Neighbor locations (footprint)


while(pixdist && reg_size < numel(I))
    
    for j = 1:4 %pointer for the neighboring pixels
        
        %get neighor pxel of pixel seed
        xn = x + neighbor(j,1);
        yn = y + neighbor(j,2);  
        

        %simple check if pixel position still inside the image
        check= (xn>=1) && (yn>=1 )&& (xn <= Isizes(1)) && (yn <= Isizes(2)); 
        
        if(check && (PLUME(xn,yn) == 0) && I(xn,yn) == 1) %check if it belongs to the thresholding boundary and if not set yet on the image we want to recreate
            neg_pos = neg_pos+1;
            neg_list(neg_pos,:) = [xn yn]; % add the new pixel
            PLUME(xn,yn)=1;
        end
        
    end
    
    %add new block of free memory
    if (neg_pos + 10 > neg_free)
        neg_free = neg_free + 10000;
        neg_list((neg_pos+1):neg_free,:) = 0;
    end
    
    PLUME(x,y) = 2; 
    reg_size=reg_size+1;
    
    % Save the x and y coordinates of the pixel (for the neighbour add proccess)
    x = neg_list(1,1);
    y = neg_list(1,2);
    

    % Remove the pixel from the neighbour (check) list
    if neg_pos > 0
        neg_list(1,:) = neg_list(neg_pos,:);
        neg_pos = neg_pos - 1;
    end
    
    if neg_pos == 0
        pixdist = 0;
    end

end

PLUME=PLUME>1;

end

function [plume_contour] = bound2coord(boundaries, lon, lat)

LON = []; LAT = [];
for jj = 1:size(boundaries,1)
    long(jj) = lon(boundaries(jj,1),boundaries(jj,2));
    lati(jj) = lat(boundaries(jj,1),boundaries(jj,2));
end
LON = [LON; long']; clear long
LAT = [LAT; lati']; clear lati

plume_contour = [LON, LAT];

end


function [row, column] = findClosestPixel(lon, lat, LON, LAT)

minDist = [abs(lon - LON) + abs(lat - LAT)];
[row, column] = find(minDist == min(abs(minDist(:))));

end
