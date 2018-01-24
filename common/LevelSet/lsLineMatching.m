function [phi_t1, val_t1, disp_points, protrusion] = lsLineMatching(varargin);
% LSLINEMATCHING matches two lines based on the Level Set method
% 
%
% SYNOPSIS       lsLineMatching
%
% INPUT                    
% 
% OUTPUT         phi_t1             : given distance function at second time step
%                val_t1(:,:,end)    : calculated distance function at second time step
%                disp_points        : are the corresponding points of the edge at the
%                                     next time step (x,y,t)
%                protrusion         : is the integrated marker path length
%                                     it positive for prot and negative for
%                                     retraction
%
% DEPENDENCES    lsLineMatching uses { 
%                                    } 
%
% Matthias Machacek 6/10/04


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=length(varargin);
for i=1:2:l
    in_found=0;
    if strcmp(varargin{i},'mask_img_t0')
        mask_img_t0=varargin{i+1};
        in_found=1; 
    elseif strcmp(varargin(i),'mask_img_t1')
        mask_img_t1=varargin{i+1};
        in_found=1; 
    elseif strcmp(varargin(i),'x_spline_t0')
        x_spline_t0=varargin{i+1}; 
        in_found=1; 
    elseif strcmp(varargin(i),'y_spline_t0')
        y_spline_t0=varargin{i+1};  
        in_found=1; 
    elseif strcmp(varargin(i),'x_spline_t1')
        x_spline_t1=varargin{i+1};  
        in_found=1; 
    elseif strcmp(varargin(i),'y_spline_t1')
        y_spline_t1=varargin{i+1};     
        in_found=1; 
    elseif strcmp(varargin(i),'known_zero_level_points_t0')
        known_zero_level_points_t0=varargin{i+1};     
        in_found=1;       
    elseif strcmp(varargin(i),'known_zero_level_points_t1')
        known_zero_level_points_t1=varargin{i+1};     
        in_found=1;         
    elseif  strcmp(varargin(i),'val_t0')
        val_t0=varargin{i+1};     
        in_found=1;         
    elseif  strcmp(varargin(i),'phi_t0')
        phi_t0=varargin{i+1};     
        in_found=1; 
    elseif  strcmp(varargin(i),'i_nn')
        r_t0=varargin{i+1};     
        in_found=1; 
    elseif  strcmp(varargin(i),'result_dir')
        RESULT_DIR=varargin{i+1};     
        in_found=1;
    elseif  strcmp(varargin(i),'control')
        CONTROL=varargin{i+1};     
        in_found=1;    
    end

    if in_found == 0
        error_string = char(varargin(i));
        error(['Unknown input:   ' , error_string]);
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MOVIE = 0;
INT_VELOCITY = 0;  
% you have to change this variable also in lsH.m !!!!
global residual_last;
global residual;
global residual_i;
residual_last = 1e10;
residual_i = 1;
if  ~exist('CONTROL','var')
    CONTROL = 1;    
end
if  ~exist('RESULT_DIR','var')
    RESULT_DIR = '/lccb/projects/alpha/Level_set_test/res/';
    %RESULT_DIR = 'L:\projects\alpha\Level_set_test\';
end
if ~exist('mask_img_t0','var')
     TEST_CASE = 1; % two ellipses
    % TEST_CASE = 2; % two offset non-intersecting circles
    % TEST_CASE = 3; % two lines    
    % TEST_CASE = 4; % line and protrusion
    % TEST_CASE = 5; % line and spike
    % TEST_CASE = 6; % seven star
    % TEST_CASE = 7; % whole cell
    % TEST_CASE = 8; % W512r_smaller
    % TEST_CASE = 9; % 
    % TEST_CASE = 10; % part cell cut_s399
    
    % computational grid size (pixel)
    x_s = 5;
    y_s = 5;
    
    sp_spacing = 2;
    
    [mask_img_t0, mask_img_t1, x_spline_t0, y_spline_t0,...
        x_spline_t1,y_spline_t1,...
        known_zero_level_points_t0, known_zero_level_points_t1,...
        grid_coordinates, domain] = lsLoadData(TEST_CASE, x_s, y_s, CONTROL);


    if 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Uses the Fast Marching Method to get the Distance Fct.  %%%%
        [img_h, img_w] = size(mask_img_t0);

        domain.x_size = img_w;
        domain.y_size = img_h;

        domain.x_spacing = 3;
        domain.y_spacing = 3; 
        % get the trail points
        % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
        % (find the known points)
        [x_X_i_t0, y_X_i_t0, x_Y_i_t0, y_Y_i_t0] = lsGetGridIntersections(x_spline_t0, y_spline_t0, domain, CONTROL);
        known_zero_level_points_t0(:,1) = [x_X_i_t0'; x_Y_i_t0'];
        known_zero_level_points_t0(:,2) = [y_X_i_t0'; y_Y_i_t0'];
   
        trial_grid_points = lsFindTrailPoints(x_X_i_t0, y_X_i_t0, x_Y_i_t0, y_Y_i_t0, domain);

        % get the distance fct values of the trail points
        trial_grid_points(:,3) = lsGetDistanceFctVec(mask_img_t0,...
            trial_grid_points, known_zero_level_points_t0, domain, 1);

        % Inizialize the distance matrix
        num_x_grid_lines = length(domain.x_grid_lines);
        num_y_grid_lines = length(domain.y_grid_lines);
        LARGE_NUMBER = 100000;
        dist_fct_tmp0(1:num_x_grid_lines,1:num_y_grid_lines) = LARGE_NUMBER;

        % get distance function by fast marching
        dist_fct_tmp = lsFastMarching(dist_fct_tmp0, trial_grid_points, domain, LARGE_NUMBER);
        figure
        surface(domain.x_grid_lines,domain.y_grid_lines,dist_fct_tmp);
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Initialize the calculation domain  %%%%%%%%%%%%%%%%%%%%%%%
    [img_h, img_w] = size(mask_img_t0);

    domain.x_size = img_w;
    domain.y_size = img_h;

    domain.x_spacing = 2;
    domain.y_spacing = 2;

    [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

    % fill the grid_line field in structure domain
    sc = 4;  % supersampling coefficient used to extract z=0 from phi
    domain = lsGenerateGridLines(domain, sc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Get edge at former time step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signed_t0 = 1;
if ~exist('phi_t0','var')
    [val_t0, phi_t0] = lsGetDistanceFct(mask_img_t0, grid_coordinates,...
                            known_zero_level_points_t0, domain, signed_t0);
end

phi_zero_t0 = lsGetZeroLevel(phi_t0, domain); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Get distance function for all grid points at present time step %%
signed_t1 = 1;
[val_t1, phi_t1] = lsGetDistanceFct(mask_img_t1, grid_coordinates,...
                             known_zero_level_points_t1, domain, signed_t1);

phi_zero_t1 = lsGetZeroLevel(phi_t1, domain); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

p_t0 = x_spline_t0.knots(1):x_spline_t0.knots(end);
x_spline_points_t0 = fnval(x_spline_t0, p_t0);
y_spline_points_t0 = fnval(y_spline_t0, p_t0);

p_t1 = x_spline_t1.knots(1) : x_spline_t1.knots(end);
x_spline_points_t1 = fnval(x_spline_t1, p_t1);
y_spline_points_t1 = fnval(y_spline_t1, p_t1);

if CONTROL & 0
    % plot results at previous time step
    h_phi_t0 = figure;
    surface(domain.x_grid_lines,domain.y_grid_lines, phi_t0);
    hold on
    plot(x_spline_points_t0, y_spline_points_t0,'r','LineWidth',2);
    contour(domain.x_grid_lines,domain.y_grid_lines, phi_t0, 40);

    % plot results at present time step
    h_phi_t1 = figure;
    surface(domain.x_grid_lines,domain.y_grid_lines,phi_t1);
    hold on
    plot(x_spline_points_t1, y_spline_points_t1,'r','LineWidth',2);
    contour(domain.x_grid_lines,domain.y_grid_lines,phi_t1, 40);
    
    % plot the two contours
    h_edges = figure;
    plot(x_spline_points_t0, y_spline_points_t0, 'b');
    hold on
    plot(x_spline_points_t1, y_spline_points_t1, 'r');
        
    % Save  results
    hgsave(h_phi_t0, [RESULT_DIR 'phi_t0.fig']); 
    print(h_phi_t0,  [RESULT_DIR 'phi_t0.eps'],'-depsc2','-tiff'); 
    print(h_phi_t0,  [RESULT_DIR 'phi_t0.tif'],'-dtiff');
    hgsave(h_phi_t1, [RESULT_DIR 'phi_t0.fig']); 
    print(h_phi_t1,  [RESULT_DIR 'phi_t0.eps'],'-depsc2','-tiff'); 
    print(h_phi_t1,  [RESULT_DIR 'phi_t0.tif'],'-dtiff');  
    hgsave(h_edges,  [RESULT_DIR 'edges_t0_t1.fig']); 
    print(h_edges,   [RESULT_DIR 'edges_t0_t1.eps'],'-depsc2','-tiff'); 
    print(h_edges,   [RESULT_DIR 'edges_t0_t1.tif'],'-dtiff');    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Advance in time             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_x = domain.x_spacing;
delta_y = domain.y_spacing;
i_end   = size(phi_t1,1);
j_end   = size(phi_t1,2);


% Put into vector
phi_t0_vec = reshape(phi_t0, prod(size(phi_t0)),1);
%track_points_0 = phi_zero_t0';

if ~exist('r_t0','var')
    r_t0=1 : sp_spacing : x_spline_t0.knots(end);
    % r_t0 =floor((x_spline_t0.knots(end)-1)/2);
end
% pixel based
%     for i=1:floor(size(known_zero_level_points_t0,1)/2)
%         track_points_0(i,1) = known_zero_level_points_t0(2*i-1,1);
%         track_points_0(i,2) = known_zero_level_points_t0(2*i-1,2);
%     end

% spline based
track_points_0(:,1) = fnval(x_spline_t0, r_t0)';
track_points_0(:,2) = fnval(y_spline_t0, r_t0)';

num_track_points = size(track_points_0, 1);

if INT_VELOCITY    
    phi_t0_vec = cat(1, phi_t0_vec,track_points_0(:,1));
    phi_t0_vec = cat(1, phi_t0_vec,track_points_0(:,2));
end

%options = odeset('RelTol',1e-4,'AbsTol',1e-7,'OutputFcn',@outputfcn,'Events',@events);
% use this output set if you want to generate a figure of the evolution
options = odeset('OutputFcn',@lsOutputfcn, 'Events', @lsEvents);
%options = odeset('Events',@lsEvents);

figure
plot(phi_zero_t0(1,:), phi_zero_t0(2,:),'r');
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Integrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    tcpu1 = clock;
    [t_steps, Y, TE,YE,IE] = ode45(@lsH,[0 200],phi_t0_vec, options,...
        phi_t1, i_end, j_end, delta_x, delta_y, domain);
    cpu_time = etime(clock,tcpu1)
else
    tcpu1 = clock;
    dt = 0.1;
    [phi, residual] = RK_4_nonTVD(phi_t0, phi_t1, dt, i_end, j_end, delta_x, delta_y, domain);
    cpu_time = etime(clock,tcpu1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(TE)
    display('Solution did not converge in given time limit');
else
    display(['Solution converged at time step  ' num2str(TE(end))]);
end
Y=Y';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%              Display results                  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


axis([0 domain.x_size 0 domain.y_size]);
h_level_sets = figure(gcf);

% Number of time steps
num_time_steps = length(t_steps);

% extract the phi and tracked points from the solution vector
for t=1:num_time_steps
    if INT_VELOCITY
        phi_vec = Y(1:i_end * j_end,t);
        phi(:,:,t) = reshape(phi_vec,i_end, j_end);
        track_points_vec = Y(i_end * j_end+1 : end,t);
        track_points(:,:,t) = reshape(track_points_vec, num_track_points, 2);
    else
        phi(:,:,t) = reshape(Y(:,t),i_end, j_end);
    end
end


if ~INT_VELOCITY 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   Get the marker paths by walking the shortest distance between the
    %%   interpolated zero level sets
    %%   The normal conditions does not work because in the case of edges the
    %%   wrong solution is chosen
    track_points(:,:,1) = track_points_0;
    h_waitbar = waitbar(0,'Integrating marker paths');
    for t=2:num_time_steps
        waitbar(t/num_time_steps, h_waitbar, num_time_steps);

        % extract zero level set using super-sampling
        zero_level_set_t = lsGetZeroLevel(phi(:,:,t),domain,1);
        %zero_level_set_t(:,end+1)=zero_level_set_t(:,1);
        % define a spline parameter based on the zero level set run length
        clear s_p;
        clear s_p_fine;
        s_p(1) = 0;
        for i=2:length(zero_level_set_t)
            s_p(i) = s_p(i-1)+sqrt((zero_level_set_t(1,i)-zero_level_set_t(1,i-1))^2+...
                (zero_level_set_t(2,i)-zero_level_set_t(2,i-1))^2);
        end
        
        % create spline of the z=0 curve at t+1 (track points are at time t
        sp_zero_level = csaps(s_p, [zero_level_set_t(1,:); zero_level_set_t(2,:)], 1);
        s_p_fine = 0:0.03:s_p(end);
        s_p_fine(end+1)=s_p(end);
        p1_sp = ppval(sp_zero_level,s_p_fine);

        for i =1:size(track_points_0,1)
            dx_l = p1_sp(1,:) - track_points(i,1,t-1);
            dy_l = p1_sp(2,:) - track_points(i,2,t-1);

            % calculate all distances
            dist_p =dx_l.^2 +dy_l.^2;

            % get the minimum distance
            [min_dist_p min_dist_p_index] = min(dist_p);

            p=ppval(sp_zero_level,s_p_fine(min_dist_p_index));
            track_points(i,1,t) = p(1);
            track_points(i,2,t) = p(2);
        end
    end
    close(h_waitbar);
    %%%%%%%%%%%%%%   End marker path integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




if CONTROL
    if 1
        % if you use this you have to select the right odeset on line ~281
        h_level_sets = figure(gcf);
        axis equal
        axis([0 domain.x_size 0 domain.y_size]);
        hgsave(h_level_sets, [RESULT_DIR 'zero_level_evol.fig']);
        print(h_level_sets,  [RESULT_DIR 'zero_level_evol.eps'],'-depsc2','-tiff');
        print(h_level_sets,  [RESULT_DIR 'zero_level_evol.tif'],'-dtiff');
    
        h_res = figure;
        plot(residual);
        xlabel('Time step');
        ylabel('Residual');
    end
end 


% extract the zero level sets
if MOVIE == 1
    h_zero_levels =figure;
    %axis([0 domain.x_size 0 domain.y_size]);
    %hold on
    %plot(x_spline_points_t0, y_spline_points_t0,'r','LineWidth',2);
    %plot(x_spline_points_t1, y_spline_points_t1,'r','LineWidth',2);
    
    i_frame = 1;     
    MakeQTMovie('start', [RESULT_DIR 'prot_movie2.mov']);
    for t=1:num_time_steps
        zero_level_set_t = lsGetZeroLevel(phi(:,:,t),domain);
        plot(zero_level_set_t(1,:), zero_level_set_t(2,:));
        axis([0 domain.x_size 0 domain.y_size]);
        MakeQTMovie('addfigure');
    end
    MakeQTMovie('finish');
end

% Initialize the displacement
%protrusion = zeros((size(Y,1) - (i_end * j_end))/2,  1);
protrusion = zeros(size(track_points,1),1);

for t=1:num_time_steps - 1
    % get the displacements
    dx = track_points(:,1,t) - track_points(:,1,t+1);
    dy = track_points(:,2,t) - track_points(:,2,t+1);
    protrusion = protrusion + sqrt(dx .^ 2 + dy .^ 2); 
end

% Get the sign of the protrusion
for i=1:size(track_points,1)
    max_l = round(track_points(i,1,end)) <= size(mask_img_t1,2) &&...
            round(track_points(i,2,end)) <= size(mask_img_t1,1);
    min_l = round(track_points(i,1,end)) >=1  && round(track_points(i,2,end)) >= 1;
    if  max_l && min_l
        if logical(mask_img_t0(round(track_points(i,2,end)), round(track_points(i,1,end)))) == 0
            prot_sign = 1;
        else
            prot_sign = -1;
        end
    else
        prot_sign = 0;
    end
       
    protrusion(i) = prot_sign * protrusion(i);
end

phi_t1_calc = phi(:,:,end);
phi_zero = lsGetZeroLevel(phi_t1_calc, domain);

disp_points = [track_points(:,1,end), track_points(:,2,end)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CONTROL
    % show level sets. This is computationally very expensive
    % so if you do not need it, switch it off.
    if 0
        for t=1:10:num_time_steps
            h_s = figure;
            surface(domain.x_grid_lines, domain.y_grid_lines, phi(:,:,t));
            hold on
            zero_level_set_t0 = lsGetZeroLevel(phi(:,:,1),domain);
            zero_level_set_t  = lsGetZeroLevel(phi(:,:,t),domain);
            zero_level_set_t1 = lsGetZeroLevel(phi(:,:,end),domain);
            plot(zero_level_set_t0(1,:), zero_level_set_t0(2,:),'g','LineWidth',1);
            plot(zero_level_set_t0(1,:), zero_level_set_t0(2,:),'g','LineWidth',1);            
            plot(zero_level_set_t(1,:), zero_level_set_t(2,:),'g','LineWidth',2);
            plot(zero_level_set_t(1,:), zero_level_set_t(2,:),'g','LineWidth',2);
            plot(zero_level_set_t1(1,:), zero_level_set_t1(2,:),'r','LineWidth',2);
            plot(zero_level_set_t1(1,:), zero_level_set_t1(2,:),'r','LineWidth',2);  
            plot([0,0],[0 200]);
            plot([200,200],[0 200]);
            plot([0,200],[0 0]);
            plot([200,0],[200 200]);
            campos(1.0e+03 *[ 1.6544   -1.2042   -0.4867]);
            set(gca,'ZDir','reverse')
            axis vis3d off
            %alpha(0.9)
            text(20,10,40,'\phi_t');
            text(170,120,0,'\Gamma_0');
            text(110,175,0,'\Gamma_1');
            text(180,160,0,'z=0');
            print(h_s, [RESULT_DIR 'level_set',num2str(t,'%04.4d'),'.tif'],'-dtiff');
            print(h_s, [RESULT_DIR 'level_set',num2str(t,'%04.4d'),'.eps'],'-depsc2','-tiff');
        end
    end
    
    % Show the time steps
    if 0
        delta_t_opt = diff(t_steps);
        av_t_steps = mean(delta_t_opt);
        h_time_steps = figure;
        plot(delta_t_opt);
        title(['Time steps,  average time step= ' num2str(av_t_steps)]);
        hgsave(h_time_steps, [RESULT_DIR 'time_steps.fig']);
        print(h_time_steps,  [RESULT_DIR 'time_steps.eps'],'-depsc2','-tiff');
        print(h_time_steps,  [RESULT_DIR 'time_steps.tif'],'-dtiff');
    end
    % Show the curvature of the last distance function
    if 0
        kappa = lsCurvature(phi(:,:,end), delta_x, delta_y, i_end, j_end);
    	figure
    	surface(kappa);
    end
    
    % Superimposed contour plots (not so intressting)
    if 0
        h_contour = figure;
        contour(domain.x_grid_lines, domain.y_grid_lines, phi_t0, 40);
        hold on
        contour(domain.x_grid_lines, domain.y_grid_lines, phi(:,:,end), 40);
        plot(phi_zero_t0(1,:), phi_zero_t0(2,:), 'g','LineWidth',2);
        plot(phi_zero_t1(1,:), phi_zero_t1(2,:), 'r','LineWidth',2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%    Track points plots     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_tracked_points = figure;
    %plot(phi_zero_t0(1,:), phi_zero_t0(2,:),'g');
    %hold on
    %plot(phi_zero(1,:), phi_zero(2,:),'r');
    %plot(phi_zero_t1(1,:), phi_zero_t1(2,:),'--m');
    plot(x_spline_points_t0, y_spline_points_t0,'g');
    hold on
    plot(x_spline_points_t1, y_spline_points_t1,'r');
    for p = 1:size(track_points,1)
        plot(squeeze(track_points(p,1,:)), squeeze(track_points(p,2,:)), '-');
    end
    axis equal
    axis([0 domain.x_size 0 domain.y_size]);
    title('Tracked points');
    %legend('Zero level t0', 'Zero level t1 solution', 'Zero level t1 given',...
    %    'org. curve t0','org. curve t1')
    hgsave(h_tracked_points, [RESULT_DIR 'tracked_points.fig']);
    print(h_tracked_points,  [RESULT_DIR 'tracked_points.eps'],'-depsc2','-tiff');
    print(h_tracked_points,  [RESULT_DIR 'tracked_points.tif'],'-dtiff');     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if MOVIE
        h_tracked_points2 = figure;
        axis([0 domain.x_size 0 domain.y_size]);
        axis equal
        plot(x_spline_points_t0, y_spline_points_t0,'g','LineWidth',1);
        hold on
        plot(x_spline_points_t1, y_spline_points_t1,'r','LineWidth',1);

        i_frame = 1;
        MakeQTMovie('start', [RESULT_DIR 'prot_movie1.mov']);
        p_max = size(track_points,3);

        for p = 1:p_max
            plot(track_points(:,1,p), track_points(:,2,p),'.', 'MarkerSize',3);
            i_frame = i_frame +1;
            MakeQTMovie('addfigure');
        end
        MakeQTMovie('finish');
        axis equal
        title('Tracked points');
        legend('Cell edge at t0','Cell edge at t1');
        hgsave(h_tracked_points2, [RESULT_DIR 'tracked_points2.fig']);
        print(h_tracked_points2,  [RESULT_DIR 'tracked_points2.eps'],'-depsc2','-tiff');
        print(h_tracked_points2,  [RESULT_DIR 'tracked_points2.tif'],'-dtiff');   
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%        Time lines       %%%%%%%%%%%%%%%%%%%%%%%%%
    h_time_lines = figure;
    plot(phi_zero_t0(1,:), phi_zero_t0(2,:),'g');
    hold on
    plot(phi_zero(1,:), phi_zero(2,:),'r');
    for t = 1:5:size(track_points,3)
        plot(track_points(:,1,t), track_points(:,2,t), '-');
    end
    axis equal
    axis([0 domain.x_size 0 domain.y_size]);
    title('Time lines');
    hgsave(h_time_lines, [RESULT_DIR 'time_lines.fig']);
    print(h_time_lines,  [RESULT_DIR 'time_lines.eps'],'-depsc2','-tiff');
    print(h_time_lines,  [RESULT_DIR 'time_lines.tif'],'-dtiff');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_disp = figure;
    plot(protrusion);
    title('Protrusion [pixel/frame rate]');
    hgsave(h_disp, [RESULT_DIR 'protrusion.fig']);
    print(h_disp,  [RESULT_DIR 'protrusion.eps'],'-depsc2','-tiff');
    print(h_disp,  [RESULT_DIR 'protrusion.tif'],'-dtiff');  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
