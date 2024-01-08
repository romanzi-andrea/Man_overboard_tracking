%% IPCV PROJECT: Man Overboard
% Group members Berasi Davide, Romanzi Andrea, José

close all
clear variables

% save cameraParams cameraParams;
load cameraParams;
tpl = padarray(ones(3,3), [3 3], 'both');
K = cameraParams.K; % intrinsic matrix
cx = K(1,3); cy = K(2,3); % camera center
depthAvgDist = 5; distBuoy = 0; % for averaging

vidReader = VideoReader("MAH01462.MP4");
H = vidReader.Height; W = vidReader.Width;

% helper function that returns y value of the line going through the points p(1,:) and p(2,:)
Y = @(c, p) p(1,2) + (p(2,2)-p(1,2))/(p(2,1)-p(1,1))*(c-p(1,1));

% ROI and tracking parameters
pb = [585, 668]; % initial position (index) of the buoy (time=1, rotated video frame)
l = 10; L = 30; % size(ROI) = (2*l+1, 2*L+1)
tl = pb-[l L]; br = pb+[l L]; % top-left and bottom-right indeces of ROI
gauss = fspecial("gaussian", 2, 2);
opticalFlow = opticalFlowFarneback;
[ROIx, ROIy] = meshgrid(-L:L, -l:l); dist = sqrt(ROIx.^2+ROIy.^2);

testing = false;
vidReader.CurrentTime = 1; oldTime = 1;
numFrames = 0;
while vidReader.hasFrame
    vidFrame = readFrame(vidReader);

    %%%%%%% Step 1: Stabilization %%%%%%%

    % Detect horizon line with Hough line detection
    edgeMap = edge(rgb2gray(vidFrame), 'canny', [0.2 0.3], 3);
    [HH,T,R] = hough(edgeMap, 'Theta', cat(2, linspace(-90, -85, 50), linspace(84, 89, 50)));
    P = houghpeaks(HH, 1);
    lines = houghlines(edgeMap,T,R,P, 'FillGap', W, 'MinLength', 100);
    horizon = lines(1);
    p = [horizon.point1; horizon.point2];

    % correct pitching (x axis fixed)
    dy = cy - Y(cx, p);
    alpha = atan(dy/K(2,2)); % we pitch the camera by -alpha
    Ralpha = [1     0           0     ;
              0 cos(alpha) -sin(alpha);
              0 sin(alpha)  cos(alpha)];
    
    % correct tilting (z axis fixed)
    M = K*Ralpha'/K;
    newp = (M*[p'; [1 1]])'; newp = newp/newp(1,3);
    newDir = newp(2,1:2) - newp(1,1:2);
    beta = atan( newDir(2) / newDir(1)); % we will tilt the camera by -beta
    Rbeta = [cos(beta) -sin(beta)  0;
             sin(beta)  cos(beta)  0;
                 0          0      1];
    
    % virtually rotate the camera
    M = K*Rbeta'*Ralpha'/K;
    tform = projtform2d(M);
    rotVidFrame = imwarp(vidFrame, tform, 'OutputView', imref2d([H W]));

    tic

    %%%%%%% Step 2: Detection of the buoy %%%%%%% 

    % computing local optical flow
    ROIoldInNew = rgb2gray(rotVidFrame(tl(1):br(1), tl(2):br(2), :)); % ROI of old frame, taken in new frame
    flow = estimateFlow(opticalFlow, ROIoldInNew);
    avgFlowx = mean(flow.Vx(:));
    avgFlowy = mean(flow.Vy(:));

    % define ROI    
    pbPred = pb + round([avgFlowy, avgFlowx]);     % predicted position of the buoy
    tl = pbPred-[l L]; br = pbPred+[l L]; % top-left and bottom-right indeces of ROI
    ROI = rgb2gray(rotVidFrame(tl(1):br(1), tl(2):br(2), :)); % new ROI, centered in pbPred
    %ROI = imfilter(ROI, gauss, 'replicate');

    a = 5; % hyper parameter
    score = (double(ROI)-mean(ROI(:)))./(dist+a);
    [maxScore, I] = max(score, [], 'all');
    buoyDetected = maxScore>4.5;

    % update buoy position
    if buoyDetected
        [r, c] = ind2sub(size(ROI),I);
        pb = tl + [r c];
    else
        pb = pbPred;
    end
    estimateFlow(opticalFlow, ROI); % insert new ROI in current frame in memory buffer


    %%%%%%% Step 3: Compute the distance of the buoy %%%%%%% 

    currDist = getBuoyDistance(pb, K); % current distance
    weight = min(numFrames, depthAvgDist);
    distBuoy = (weight*distBuoy + currDist)/(weight+1);

    1 / toc
    %%%%%%% Step 4: Show the frame %%%%%%%

    t=tiledlayout(3, 2);
    %title(t, sprintf("Time = %0.3g  alpha=%0.3g° beta=%0.3g°, Distance=%0.3g", ...
    %    vidReader.CurrentTime, alpha*180/pi, beta*180/pi, distBuoy));

    nexttile([2 2])
    imshow(rotVidFrame); hold on
    line([1 W], [cy cy], 'Color', 'green'); % target line for horizon
    [newX, newY] = transformPointsForward(tform, p(:,1), p(:,2)); newp = [newX';newY']';
    plot(cx, cy, Marker='+', Color='green') % Focal center of the camera
    
    if buoyDetected color='red'; else color='magenta'; end
    plot(pb(2), pb(1), Marker="+", Color=color);
    rectangle('Position', [tl(2) tl(1) (2*L+1) (2*l+1)], 'EdgeColor', color)
    text(tl(2)-50, tl(1)+70, sprintf("%0.3g m", distBuoy), "FontSize", 15, 'Color', color)
    hold off
    title("Stabilized frame");

    nexttile; 
    imshow(ROI); hold on
    if buoyDetected
        plot(c, r, Marker='+', Color='red');
        viscircles([c r], 2, 'Color', 'red', 'LineWidth', 1);
    end

    %plot(flow, 'DecimationFactor', [5 5], 'ScaleFactor', 1);
    quiver(L+1, l+1, avgFlowy, avgFlowx, 'Color', 'green');     hold off
    title("ROI")

    nexttile;
    imshow(score.^2, []); hold on;
    plot(L+1, l+1, 'Marker','+', 'Color','blue')
    if buoyDetected
        plot(c, r, Marker='+', Color='red');
        viscircles([c r], 2, 'Color', 'red', 'LineWidth', 1);
    end
    hold off;
    title("SCORE");

    t.Padding = 'compact'; %t.TileSpacing = 'compact';
    drawnow

    if testing
        pause;
    else
        pause(0)
    end
    numFrames = numFrames + 1;
    
end