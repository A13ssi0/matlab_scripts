function [lapMask, weighted_lapMask, chanlocs] = get_laplacianMask(channels, distance, isDistanceMeasure, chanlocs)
% [lapMask, weighted_lapMask] = get_laplacianMask(channels, distance, isDistanceMeasure, chanlocs) 
%      returns the classic laplacian mask, weighted on the distance and the chanlocs .
%      Inputs :     - channels = cell array of channels name strings. The lapMask is based on this channels order
%                       - distance = integer of neighbours (default:all) or float cutoff distance (default:1) to take into account. The channels distance is normalized 
%                                       to 1 (between FPz and Oz or the greatest distance between channels if another chanlocs is given)
%                       - isDistanceMeasure (defaul:false) = flag to select if variable distance is the number of neighbours (default) or the cutoff distance
%                       - chanlocs (defaul:antNeuro 64) = custom struct chanlocs with at least column names 'labels','X','Y','Z'
    
    if nargin<3 || isempty(isDistanceMeasure),      isDistanceMeasure=false;    end
    if nargin<3,    chanlocs = get_standard_chanlocs(); end
    if nargin<1 || isempty(channels),    channels = {chanlocs.labels}; end
    if nargin<2 || isempty(distance),    distance = length(channels)-1; end

    [labels, idx_chan, idx_locs] = intersect(lower(channels), lower({chanlocs.labels})); 
    if length(labels)~=length(channels) 
        not_found = true(1,length(channels));
        not_found(idx_chan) = false;
        warning([num2str(length(channels)-length(labels)), ' channel not found in chanlocs.   Not found: ', strjoin(channels(not_found))])
    end

    % get distance between channels + normalization + only wanted channels
    pos_chanlocs = [chanlocs.X; chanlocs.Y; chanlocs.Z]';
    locs_dist = pdist2(pos_chanlocs,pos_chanlocs);
    locs_dist = locs_dist/max(locs_dist,[],'all');
    locs_dist = locs_dist(idx_locs,idx_locs);

    % reordering the locs based on channels
    chan_dist(idx_chan,idx_chan) = locs_dist;

    % removing channels that does not meet requirements
    if isDistanceMeasure
        chan_dist(chan_dist>distance) = 0;
    else
        n_max = length(channels)-(distance+1);      % take into account channels itself
        [~,idx_max] = maxk(chan_dist,n_max,2);
        for n_chan = 1:length(channels)
            chan_dist(idx_max(n_chan,:),n_chan) = 0;
        end
    end

    weighted_lapMask = zeros(size(chan_dist));
    lapMask = zeros(size(chan_dist));
    % normalize matrix where diag=1 and sum column = 0
    for n_chan = 1:length(channels)
        col = chan_dist(:,n_chan);

        col(col>0) = -1./col(col>0);
        weighted_lapMask(:,n_chan) = -col/sum(col);
        weighted_lapMask(n_chan,n_chan) = 1;

        col(col>0) = -1/sum(col>0);
        lapMask(:,n_chan) = col;
        lapMask(n_chan,n_chan) = 1;
    end

    chanlocs = chanlocs(idx_locs);
    chanlocs(idx_chan) = chanlocs;

end






function chanlocs = get_standard_chanlocs()
        chanlocs = {'FP1',-19.3,0.525,83.9,29.4,-6.99,19.3,-4.49,89.2;'FPZ',0.0729,0.506,88.2,-0.112,-1.71,-0.0729,-1.11,88.3;'FP2',19.4,0.525,84.9,-29.9,-7.08,-19.4,-4.50,90.3;'F7',-58.8,0.544,42.5,70.3,-11.4,58.8,-7.92,82.9;'F3',-43.4,0.333,53.1,50.2,42.2,43.4,30,84.4;'FZ',0.306,0.230,58.5,-0.312,66.5,-0.306,48.6,88.5;'F4',43.7,0.341,54.3,-51.8,40.8,-43.7,28.5,85.5;'F8',58.7,0.544,44.4,-73,-12,-58.7,-7.99,86.3;'FC5',-76.4,0.405,18.6,77.2,24.5,76.4,17.1,83.1;'FC1',-52.6,0.157,26,34.1,80,52.6,61.8,90.7;'FC2',52.8,0.161,26.4,-34.8,78.8,-52.8,61,90.1;'FC6',75.9,0.408,19.9,-79.5,24.4,-75.9,16.6,85.6;'M1',-118,0.694,-45,86.1,-68,118,-35,119;'T7',-101,0.535,-16,84.2,-9.35,101,-6.23,86.2;'C3',-100,0.255,-11.6,65.4,64.4,100,44.1,92.5;'CZ',177,0.0291,-9.17,-0.401,100,-177,84.8,101;'C4',99.2,0.261,-10.9,-67.1,63.6,-99.2,43.1,93.1;'T8',100,0.535,-15,-85.1,-9.49,-100,-6.27,86.9;'M2',118,0.695,-45,-85.8,-68,-118,-35.1,118;'CP5',-120,0.397,-46.6,79.6,30.9,120,18.6,97.3;'CP1',-143,0.183,-47.3,35.5,91.3,143,57.1,109;'CP2',141,0.188,-47.1,-38.4,90.7,-141,56.2,109;'CP6',119,0.399,-46.1,-83.3,31.2,-119,18.1,100;'P7',-135,0.508,-73.5,72.4,-2.49,135,-1.38,103;'P3',-146,0.331,-78.8,53,55.9,146,30.5,110;'PZ',180,0.247,-81.1,-0.325,82.6,-180,45.5,116;'P4',145,0.331,-78.6,-55.7,56.6,-145,30.4,112;'P8',135,0.508,-73.1,-73.1,-2.54,-135,-1.41,103;'POZ',180,0.354,-102,-0.216,50.6,-180,26.3,114;'O1',-165,0.476,-112,29.4,8.84,165,4.35,117;'O2',165,0.476,-112,-29.8,8.80,-165,4.34,116;'AF7',-38.7,0.538,68.6,54.8,-10.6,38.7,-6.88,88.4;'AF3',-23.7,0.421,76.8,33.7,21.2,23.7,14.2,86.5;'AF4',24.7,0.420,77.7,-35.7,22,-24.7,14.4,88.3;'AF8',38.7,0.538,69.7,-55.7,-10.8,-38.7,-6.87,89.9;'F5',-53.3,0.434,48,64.5,16.9,53.3,11.9,82.2;'F1',-25.8,0.257,56.9,27.5,60.3,25.8,43.7,87.4;'F2',27.1,0.263,57.6,-29.5,59.5,-27.1,42.6,87.9;'F6',53.7,0.439,49.8,-67.9,16.4,-53.7,11,85.8;'FC3',-69.3,0.273,22.7,60.2,55.5,69.3,40.8,85;'FCZ',0.787,0.0954,27.4,-0.376,88.7,-0.787,72.8,92.8;'FC4',69.2,0.279,23.7,-62.3,55.6,-69.2,39.8,86.8;'C5',-99.7,0.391,-13.8,80.3,29.2,99.7,19.7,86.5;'C1',-105,0.126,-9.98,36.2,89.8,105,67.3,97.3;'C2',104,0.132,-9.62,-37.7,88.4,-104,66.3,96.6;'C6',98.7,0.394,-12.8,-83.5,29.2,-98.7,19.1,89.3;'CP3',-126,0.279,-47,63.6,65.6,126,39.7,103;'CP4',125,0.284,-46.6,-66.6,65.6,-125,38.9,104;'P5',-139,0.413,-76.3,67.3,28.4,139,15.6,106;'P1',-160,0.270,-80.5,28.6,75.4,160,41.4,114;'P2',158,0.269,-80.5,-31.9,76.7,-158,41.5,116;'P6',138,0.414,-75.9,-67.9,28.1,-138,15.4,106;'PO5',-154,0.439,-99.3,48.4,21.6,154,11.1,113;'PO3',-160,0.394,-101,36.5,37.2,160,19.1,114;'PO4',160,0.396,-101,-36.8,36.4,-160,18.7,113;'PO6',153,0.439,-99.4,-49.8,21.7,-153,11.1,113;'FT7',-80.1,0.543,14.1,80.8,-11.1,80.1,-7.73,82.8;'FT8',79.3,0.543,15.4,-81.8,-11.3,-79.3,-7.75,84;'TP7',-118,0.523,-46,84.8,-7.06,118,-4.18,96.8;'TP8',118,0.523,-45.5,-85.5,-7.13,-118,-4.21,97.2;'PO7',-151,0.492,-97.5,54.8,2.79,151,1.43,112;'PO8',150,0.492,-97.6,-55.7,2.73,-150,1.39,112;'OZ',180,0.460,-115,-0.108,14.7,-180,7.27,116};
        names = {'labels','theta','radius','X','Y','Z','sph_theta','sph_phi','sph_radius'};
        chanlocs = cell2struct(chanlocs,names,2);

        % adding cpz
        cp1 = chanlocs(strcmpi({chanlocs.labels},'CP1'));
        cp2 = chanlocs(strcmpi({chanlocs.labels},'CP2'));
        cz = chanlocs(strcmpi({chanlocs.labels},'Cz'));
        pz = chanlocs(strcmpi({chanlocs.labels},'Pz'));
        cpz_x = mean([cp1.X cp2.X cz.X pz.X]);
        cpz_y = mean([cp1.Y cp2.Y cz.Y pz.Y]);
        cpz_z = mean([cp1.Z cp2.Z cz.Z pz.Z]);
        [cpz_sph_theta,cpz_sph_phi,cpz_sph_r] = cart2sph(cpz_x,cpz_y,cpz_z);
        cpz_sph_theta = cpz_sph_theta*180/pi;
        cpz_sph_phi = cpz_sph_phi*180/pi;
        % point of reference is close to Cz but between it and FCz (theta_Cz~180)
        cpz_theta = 180;    % all in the same line are 180  
        cpz_radius = (0.247+0.029)/2;   %mean between Cz and Pz. It makes sense because is close to the difference between FCz and Cz
        cpz = {'CPz', cpz_theta,cpz_radius, cpz_x, cpz_y, cpz_z, cpz_sph_theta, cpz_sph_phi, cpz_sph_r};
        chanlocs(end+1) = cell2struct(cpz,names,2);
end
