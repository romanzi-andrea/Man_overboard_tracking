function buoyDistance = getBuoyDistance(pb, K)
%getBuoyDistance Computes the distance of the buoy.
%   It computes the distance of the buoy that is in position pb 

cx = K(1,3); cy = K(2,3); % camera center
ch = 2.5; % height of the camera
earthRadius = 6371000;
dHorizon = sqrt((earthRadius+ch)^2-earthRadius^2);   %horizon distance

% distance in y
y0 = pb(1) - cy;       % diff buoy horizon in the y axis 
gammay = atan(y0/K(2,2));     %angle y between z camera axis and buoy
eta = atan(earthRadius/dHorizon);                    % angle between horizon line and vertical line of the camera
dYBuoy = ch*tan(eta-gammay);                         % distance in y of the buoy

% distance in x
x0 = pb(2) - cx;        % diff buoy horizon in the x axis
dXBuoy = x0 * dYBuoy/K(1,1);    

buoyDistance = sqrt(dXBuoy^2 + dYBuoy.^2);
end