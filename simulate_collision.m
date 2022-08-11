function simulate_collision(varargin)

    %% The main function to simulate collision, as shown in Fig. 1 
    % Here we simulate objects with random position and velocity around an
    % observer that is moving forward, and calculate their future collision
    % risks as the time-dicounted inverse intercepts of their trajectories.
    % We then plot the mean collision risk as a function of instantaneous
    % apparent angular size and velocity of the objects.

    %% Default parameters (used in the paper)
    % These can be overwritten by providing input arguments
    
    random_seed = 1;     % for reproducibility
    n_sample = 5000000;  % #objects simulated
    observer_speed = 10; % in Y dimension, mm/s
    object_size = 2;     % in mm (radius)
    r_max = 200;         % in mm. Max distance from the observer to the object
    vmin = 0;            % in mm/s. Object speed is uniformly sampled between vmin/vmax
    vmax = observer_speed * 2; 
    tau = 10;            % in s. This defines the time horizon for the collision risk calculation
    epsilon = object_size;
    
    % overwrite default parameters by  input arguments 
    for ii = 1:2:length(varargin)
        eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
    end
    
    %% Initialze the random number generator
    rng(random_seed);

    %% Run simulation
    % initialize the position and velocity of the objects
    [x,y]   = seed_position_uniform(r_max,n_sample);
    [vx,vy] = seed_velocity_uniform(vmin,vmax,n_sample); 

    % account for the forward locomotion of the observer
    vy = vy - observer_speed; 

    % calculate apparent angular position, size, and their rates of change
    [phi,phi_dot,psi,psi_dot] = calc_observables(x,y,vx,vy,object_size);
    
    % calculate collision risk
    CR = calc_collision_risk(x,y,vx,vy,tau,epsilon);

    %% Visualize
    % visualize collision risk as a function of position and velocity
    [map,counts,xx,yy] = my2dHeatmap((phi),(phi_dot),CR);
    figure; hold on
    imagesc(xx,yy,map); contour(xx,yy,map,'w','LineWidth',1);
    axis xy; xlabel('\phi'); ylabel('d\phi/dt'); title('Immediate Collision Risk')

    % visualize example trajectories with high risk and high btf velocity
    % in the observer centric coordinate 
    folded_phi_dot = phi_dot .* sign(phi); % + = ftb, - = btf
    temp = find(CR>0.2 & folded_phi_dot<-20);
    plot_example_trajectory(temp,x,y,vx,vy,observer_speed)
    xlabel('x (mm)'); ylabel('y (mm)'); title('Objects with high collision risk and high btf velocity')

    license('inuse')
end


%% Subroutines
%% Sampling functions
function [x,y,r,phi] = seed_position_uniform(r_max,N)
    % seed the object position uniformly within the circular area within
    % r_max distance from the observer
    r = rand(N,1)+rand(N,1);
    r = (2-r).*(r>1) + r.*(r<=1);
    r = r*r_max;
    phi = 360*rand(N,1);
    % convert to cartesian
    x = sind(phi).*r;
    y = cosd(phi).*r;
end

function [vx,vy] = seed_velocity_uniform(vmin,vmax,N)
    % seed the object velocity 
    % direction is uniform, speed is bounded uniform
    speed = (rand(N,1)*(vmax-vmin))+vmin;
    heading = 360*rand(N,1);
    % convert to cartesian
    vx = sind(heading).*speed;
    vy = cosd(heading).*speed;
end


%% Function to calculate collision risk and observables

% Collision risk defined as inverse distance at future frontal path crossing
% under the constant velocity assumption.
function h = calc_collision_risk(x,y,vx,vy,tau,epsilon)
    if nargin<6
        epsilon = 0;
    end
    if nargin<5
        tau = +Inf;
    end
    % (signed) distance at path crossing (= intercept)
    intercept = (y.*vx - x.*vy)./vx;
    % timing of path crossing
    T = -x./vx;
    % risk is the inverse of future positive intercept
    % epsilon is there to avoid zero division
    h = ((T>=0) & (intercept>0)).*(1./(intercept+epsilon)).*exp(-(T.*(T>0))/tau);
end

% Calculate observable variables
function [phi,phi_dot,psi,psi_dot] = calc_observables(x,y,vx,vy,R)
    % current distance
    d = sqrt(x.^2 + y.^2);
    % angular position 
    phi = atan2(x,y);
    % angular velocity
    phi_dot = (vx.*y - vy.*x)./(d.^2);
    % angular extent    
    psi = 2*atan(R./d);
    % angular expansion
    psi_dot = -2*R.*(x.*vx + y.*vy)./d./(d.^2 + R.^2);
    
    % convert to degrees
    phi = phi/pi/2*360;
    phi_dot = phi_dot/pi/2*360;
    psi = psi/pi/2*360;
    psi_dot = psi_dot/pi/2*360;
end

%% visualization function
% calculate mean 2d heatmap
function [M,counts,xmeshcenter,ymeshcenter] = my2dHeatmap(X,Y,C,varargin)

prctileThresh = 2.5; % ignore outliers
nMesh = 21;
xrange = [];
yrange = [];

for ii = 1:2:length(varargin)
   eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
end

if isempty(xrange) 
    xrange = [prctile(X(:),prctileThresh),prctile(X(:),100-prctileThresh)];
end
if isempty(yrange)
    yrange = [prctile(Y(:),prctileThresh),prctile(Y(:),100-prctileThresh)];
end

xmesh = linspace(xrange(1),xrange(2),nMesh);
ymesh = linspace(yrange(1),yrange(2),nMesh);

xmeshcenter = (xmesh(1:end-1)+xmesh(2:end))/2;
ymeshcenter = (ymesh(1:end-1)+ymesh(2:end))/2;

M = nan(nMesh-1,nMesh-1);
counts = nan(nMesh-1,nMesh-1);


for xx = 1:nMesh-1
    for yy = 1:nMesh-1
        if any(X>xmesh(xx) & X<xmesh(xx+1) & Y>ymesh(yy) & Y<ymesh(yy+1))
           M(yy,xx) = nanmean(C(X>xmesh(xx) & X<xmesh(xx+1) & Y>ymesh(yy) & Y<ymesh(yy+1)));
           counts(yy,xx) = sum(X>xmesh(xx) & X<xmesh(xx+1) & Y>ymesh(yy) & Y<ymesh(yy+1));
        end
    end
end

end

function plot_example_trajectory(index,x,y,vx,vy,vo)
% Given sample vectors of position and velocity x, y, vx, vy
% Plot vectors corresponding to samples indicated by the index array
% compensate for observer velocity vo
if nargin<6
    vo = 0;
end
vy = vy + vo;

figure; hold on
quiver(0,0,0,vo/10,'ro','AutoScale','off');
quiver(x(index),y(index),vx(index)/10,vy(index)/10,'bo','AutoScale','off');

end


