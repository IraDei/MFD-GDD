function [out] = MFD_GDD(img, l1, l2, Sr)
%MFD_GDD 
% Realization Code for"Infrared Small Target Detection Based on Flux Density and Direction Diversity in Gradient Vector Field," 
% by D. Liu, L. Cao, Z. Li, T. Liu and P. Che, 
% in IEEE Journal of Selected Topics in Applied Earth Observations and
% Remote Sensing, vol. 11, no. 7, pp. 2528-2554, July 2018.
% Source paper doi: 10.1109/JSTARS.2018.2828317.
%
% By IraDei, Jul.15.2020
% PARAM LIST:
%   'l1, l2' - Min and max Gradient divergence range in MFD computation.
%   'Sr' - structure tensor scale enums neighbor whose Manhattan dist towards kernel pixel. 

% convert input frame into double format
if(ischar(img))
    img = imread(img);
end
[imgR ,imgC, dimension]=size(img);
if(dimension>2) 
    img = rgb2gray(img);
    dimension = 1;
end

% init global variables
img = im2double(img);
% init Sample Offest of k-scale Girdle surround kernel pixel
lss = l2 - l1 + 1;  % girdle scale size
smp_ofs_l = cell(lss, 1);   % sample grid offset and gradient flag
smp_rng_l = l1:1:l2;
smp_pix_sz = smp_rng_l * 8;
smp_desert_k = 4;
smp_dist_decay = smp_pix_sz - smp_desert_k; % decay factor of Gradient divergence alone radius lk 
for i = 1:lss
    smp_ofs_l{i} = getGirdleSampleOffset(smp_rng_l(i)); % Gradient sample offset on Li
end

% compute Infrared Gradient Vector Field(IGVF), derivative in x and y
% direction infact, according to Eq.3~5 given in section II.
hx = [-1 0 1];
hy = hx';

% I(x,y) = [dy(x,y), dx(x,y)]'
dx = imfilter(img, hx, 'corr', 'same', 'replicate');
dy = imfilter(img, hy, 'corr', 'same', 'replicate');

% compute Gradient Divergence according to Eq.8~13
% The total complexity seems to be quite high under serial mode

% Gradient divergence of each scale
gradDiv_k = zeros(imgR, imgC, lss);
gradDiv_min = Inf(imgR, imgC);    % minimum Gradient divergence among all scales
parfor y = 1:imgR
    for x = 1:imgC
        for k = 1:lss
            % sample girdle pixels at range lk
            smp_pqty = smp_pix_sz(k);
            smp_ofs = smp_ofs_l{k};
            smp_pix = zeros(smp_pqty,3);   % I(x,y) = [dy(x,y), dx(x,y)]' and norm of I(x,y)
            smp_gnc = zeros(smp_pqty, 2);
            smp_idx = 0;
            for i = 1:smp_pqty
                smpy = y + smp_ofs(i,1);
                if(smpy>0 && smpy<=imgR)
                    smpx = x + smp_ofs(i,2);
                    if(smpx>0 && smpx<=imgC)
                        smp_idx = smp_idx + 1;
                        smp_pix(smp_idx, 1) = dy(smpy, smpx);
                        smp_pix(smp_idx, 2) = dx(smpy, smpx);
                        % normalization is time cosuming so we use abs sum
                        % here to boost total procession by 6 secs
                        %smp_pix(smp_idx, 3) = norm(smp_pix(smp_idx, 1:2) ,2);
                        smp_pix(smp_idx, 3) = abs(dy(smpy, smpx)) + abs(dx(smpy, smpx));
                        
                        % Gradient factor 'nc' of sample pixel
                        smp_gnc(smp_idx, 1) = smp_ofs(i,3); % ncx in Eq.10
                        smp_gnc(smp_idx, 2) = smp_ofs(i,4); % ncy in Eq.11
                    end
                end
            end
            smp_pix = smp_pix(1:smp_idx, :);
            smp_gnc = smp_gnc(1:smp_idx, :);
            [smp_pix_sort, smp_ord] = sort(smp_pix(:,3), 'descend');
            
            % calculate gradient divergence
            for i = smp_desert_k + 1 : smp_idx
                idx = smp_ord(i);
                % dot product costs much time here, which does need parallel optimization 
                gradDiv_k(y, x, k) = gradDiv_k(y, x, k) + smp_pix(idx,1:2) * smp_gnc(idx,:)';
            end
        end
    end
end

% Multi-scale flux density
for k = 1:lss
    gradDiv_k(:, :, k) = gradDiv_k(:, :, k) / smp_dist_decay(k);
    for y = 1:imgR
        for x = 1:imgC
            gdk = gradDiv_k(y, x, k);
            if(gradDiv_min(y,x)>gdk)
                gradDiv_min(y,x) = gdk;
            end
        end
    end
end
% test passed, valid MFD Saliency
MFD = gradDiv_k(:,:, 1) - gradDiv_min(:,:); 

% We shall continue realization from GDD saliency of Eq.14 later
% Compute Structural Saliency of Gradient vector field
S = zeros(imgR, imgC, 2, 2);    % 2*2 structural saliency matrice for all piexls
Wr = 1/(Sr*2 + 1)^2;
for y = 1:imgR
    for x = 1:imgC
        Sxx = 0;
        Sxy = 0;
        Syy = 0;
        for ofs_y = -Sr:Sr
            smp_y = max(1, min(y + ofs_y, imgR));
            for ofs_x = -Sr:Sr
                smp_x = max(1, min(x + ofs_x, imgC));
                Sxx = Sxx + dx(smp_y, smp_x)^2;
                Sxy = Sxy + dx(smp_y, smp_x) * dy(smp_y, smp_x);
                Syy = Syy + dy(smp_y, smp_x)^2;
            end
        end
        S(y, x, :, :) = [Sxx Sxy; Sxy Syy];
    end
end
S = S / Wr;

% generate GDD Saliency based on Strcutural saliency
% And Weight MFD map with GDD saliency
GDD = zeros(imgR, imgC);
WMFD = zeros(imgR, imgC);
ptb_e = 10^-9;  % perturbation term given in Eq.19
parfor y = 1:imgR
    for x = 1:imgC
        % compute eigen val and eigen vector of Structural matrix by pixel
        [egval egvec] = eig(reshape(S(y, x, :, :),2,[]));
        GDD(y, x) = exp(-( (egval(1) - egval(2) + ptb_e)/(egval(1) + egval(2) + ptb_e) )^2);
        WMFD(y, x) = GDD(y, x) * MFD(y, x);
    end
end

% Threshold WMFD map
bw = out_bw(WMFD, 7);

% visualization
figure(1);
subplot(1,4,1);
imshow(MFD,[]);
subplot(1,4,2);
imshow(GDD,[]);
subplot(1,4,3);
imshow(WMFD,[]);
subplot(1,4,4);
imshow(bw,[]);
end


