function [zlayt,zlayb,vlay,varargout] = layerise(Z,V,dvmin,ifplot,varargin)


% Code to take a continuous series of velocity (V) with depth (Z) and
% convert it to several layers of constant velocity. 
% New approach (brb2022.06.29) is to traverse down the model and make a new
% layer when either the velocity or depth (or norm of those two) has
% changed suffiently. 
% 
% INPUTS:
%   Z     - [N x 1] vector of depths (km)
%   V     - [N x 1] vector of velocities (km/s)
%   dvmin - brb2022.06.29 removed temporarily. 
%           Tunable variable specifying layer thickness. If this value is
%           smaller, the continuous profile will be better fit, but with
%           more layers. If the value is larger, the fit will be poorer,
%           but there will be fewer layers. Recommended value: 0.1 (km/s)
%   ifplot - option to plot (true) or not (false)
%   varargin - any number of [N x 1] vectors for other parameters (e.g. Vp, 
%           rho) to be discretised onto the same mesh. The discontinuity 
%           depths will be well-preserved, but because the algorithm isn't
%           actually adapting to these parameters, the jumps across the
%           discontinuities will be less precisely recovered. 
% 
% OUTPUTS:
%   zlayt - [Nlay x 1] vector of depths of layer tops (km)
%   zlayb - [Nlay x 1] vector of depths of layer bases (km)
%   vlay  - [Nlay x 1] vector of velocity values for each layer (km/s)
%   varargout - any number of [Nlay x 1] vector of other parameter values
%           in each layer

%     Z. Eilon 08/2016
%     Modified brb2022.06.29. 

break_v    = 0.01 ; % Tend toward break into a layer every this much velocity change. 
break_z    = 10   ; % Tend toward making all layers smaller than this. 
% break_v    = 0.05 ; % Tend toward break into a layer every this much velocity change. 
% break_z    = 30   ; % Tend toward making all layers smaller than this. 
break_grad = 0.025; % If gradient in v between two provided depths is greater than this, divide a layer here. Mostly this value is irrelevant, we just identify grad==inf, i.e. discontinuities. 
% Note: The new model will not get finer than the provided Z array. 

if nargin<4 || isempty(ifplot)
    ifplot=false;
end

dz = diff(Z);
dv = diff(V);
dvdz = dv./dz; % Gradient. 
is_disc = find(isinf(dvdz)); % Top index of discontinuities. dz = 0. 

% Loop to find the depths where we will break up layers. 
% Put in a layer where there is strong gradients. 
% Or if we haven't put in a layer for a long distance. 
% Or a combination of the two. 
% Manually assign first and last Z as layer breaks. 
break_inds = [1]; % Indecies where we will break into layers. 
for iz = 2:(length(Z)-1);

    if dvdz(iz) > break_grad; % Make sure discontinuities (e.g. grad == inf) have a layer break. 
        break_inds = [break_inds; iz]; 
        continue; 
    end 
    
    vchange = abs(V(iz) - V(break_inds(end))); 
    zchange = Z(iz)     - Z(break_inds(end) ); 
    
    tendchange = sqrt( (vchange / break_v)^2 ...
                      +(zchange / break_z)^2); % Will be 1 if vchange or zchange meets our predefined values of when to break layers, but also will be one if a combination of v and z change is high enough. 
    
    if tendchange > 1; 
        break_inds = [break_inds; iz]; 
    end
end
break_inds = unique([break_inds; length(Z)]); % Manually assign last Z as a layer break. 

zlayt = Z(break_inds(1:end-1)); % Top depth of each layer. 
zlayb = Z(break_inds(2:end  )); % Bottom depth of each layer. 

dup_vals = (zlayt == zlayb); % This can catch duplicate depths in one array - Remove accidentally added infinitesimally small layers. 
zlayt = zlayt(~dup_vals); 
zlayb = zlayb(~dup_vals); 

% Look for some obvious signs of mistakes.  
layer_gaps = zlayt(2:end) ~= zlayb(1:end-1); 
if any(layer_gaps); 
    error('The bottom of some layer is not the same as the top of the next layer'); 
elseif Z(end) ~= zlayb(end);
    error('Bottom of new layers not equal to bottom of old layers'); 
elseif Z(1) ~= zlayt(1); 
    error('Top of new layers not equal to top of old layers'); 
end

%%% Get ready to take average of velocity and whatever within each layer. 
nlay               = length(zlayt); % Number of layers in new model
vlay               = nan(size(zlayt)); % Output velocity

% Keep track of whether above or below layer boundary. For discontinuities. 
Z_shift            = Z; 
Z_shift(is_disc  ) = Z_shift(is_disc  ) - 0.00001; 
Z_shift(is_disc+1) = Z_shift(is_disc+1) + 0.00001; 

% Prepare varagin/out for processing.
process_varargin = ~isempty(varargin);  
if process_varargin; % Prepare to process varargin, if it was provided. If needed for testing: varargin = [{V+1} {V+2}]; 
    nargs = length(varargin); 
    varargin_arr = nan(size(Z,1), nargs); % Make into an array for vectorized processing    
    for iv = 1:nargs; 
        varargin_arr(:,iv) = varargin{iv}; 
    end
    varargout_arr = nan(nlay, nargs); 
else
    varargout = {}; 
end

% Important loop - get average of v and varargin within each layer. If needed for testing: disp([Z(in_lay)' ; V(in_lay)']); % See what depths and velocities you are taking. 
for iz = 1:length(zlayt); 
    in_lay = and( ( Z_shift >= zlayt(iz) ) , ( Z_shift <= zlayb(iz) ) ); % Within this layer or not? Note: This realies on Z_shift having discontinuity depths modified. 
    vlay(iz) = mean(V(in_lay)); % Mean of velocity from top to bottom of this layer. 
    if process_varargin; 
        varargout_arr(iz, :) = mean(varargin_arr(in_lay,:),1); % Mean of all varargin varaibles, from top to bottom of this layer. Take mean along first dimension.  
    end
end

if any(isnan(vlay)); 
    error('Did not get velocity for some layer. '); 
end

% Convert varargout_arr to an array of cells. 
if process_varargin; 
    if any(isnan(varargout_arr)); 
        error('A varargout value is nan. Might not have been averaged or provided correctly?'); 
    end
    
    varargout = {}; 
    for iarg = 1:nargs; 
        varargout(:,end+1) = {varargout_arr(:,iarg)}; 
    end
end
         
ifplot = false; if ifplot; warning('Setting ifplot = true'); end; 
if ifplot
    figure(11); clf; set(gcf,'pos', [-581 247 514 796], 'color', 'white'); 
    h = tiledlayout(1, 1+length(varargin),'TileSpacing','compact'); 
    nexttile; hold on; box on; 

    plot(V,Z,'-ko')
    zlayp = reshape([zlayt';zlayb'],2*nlay,1);
    vlayp = reshape([vlay';vlay'],2*nlay,1);
    plot(vlayp,zlayp,'-ro')
    set(gca,'ydir','reverse',...
        'ylim',[0, max(Z)],'xlim',[0.9*min(V) 1.1*max(V)])
    
    % other vars
    for iv = 1:length(varargin)
        nexttile; hold on; box on; 
        plot(varargin{iv},Z,'-ko')
        ivlayp = reshape([varargout{iv}';varargout{iv}'],2*nlay,1);
        plot(ivlayp,zlayp,'-ro')
        set(gca,'ydir','reverse','ylim',[0, max(Z)],...
            'xlim',[0.9*min(varargin{iv}) 1.1*max(varargin{iv})])    
    end
    text(min(V), 250, sprintf('nlay=%1.0f\nbreak z=%1.1f\nbreak v=%1.3f',...
        nlay,break_z,break_v), 'fontsize', 12); 
end

end % End function
    