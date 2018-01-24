function [mask_img_t0, mask_img_t1,x_spline_t0, y_spline_t0, x_spline_t1,y_spline_t1,...  
    known_zero_level_points_t0, known_zero_level_points_t1, grid_coordinates, domain] = lsLoadData(TEST_CASE, x_s, y_s, CONTROL)


if TEST_CASE == 1 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  two ellipses
   domain.x_size = 200;
   domain.y_size = 200;
   
   domain.x_spacing = x_s; 
   domain.y_spacing = y_s;
   
   % supersampling coefficient
   sc = 5;
   
   %create circle
   circle       = rsmak('circle',50,[0, 0]);
   ellipse_t0   = fncmb(circle,[1.5 0;0 0.6]);
   ellipse_t1   = fncmb(circle,[0.6 0;0 1.5]);
   
   ellipse_t0   = fncmb(ellipse_t0,'+', 100);
   ellipse_t1   = fncmb(ellipse_t1,'+', 100);
   fnplt(ellipse_t0);
   hold on
   fnplt(ellipse_t1);
   axis equal
   
   % Create mask
   p = 0:0.05:ellipse_t0.pieces;
   ellipse_points_t0 = fnval(ellipse_t0,p);
   p = 0:0.05:ellipse_t1.pieces;
   ellipse_points_t1 = fnval(ellipse_t1,p);
   
   mask_img_t0 = roipoly(domain.y_size, domain.x_size, ellipse_points_t0(1,:)', ellipse_points_t0(2,:)');
   mask_img_t1 = roipoly(domain.y_size, domain.x_size, ellipse_points_t1(1,:)', ellipse_points_t1(2,:)');

   % create x,y splines
   s_p = 1:length(p);
   x_spline_t0 = fn2fm(spline(s_p, ellipse_points_t0(1,:)),'B-');
   y_spline_t0 = fn2fm(spline(s_p, ellipse_points_t0(2,:)),'B-');
   
   x_spline_t1 = fn2fm(spline(s_p, ellipse_points_t1(1,:)),'B-');
   y_spline_t1 = fn2fm(spline(s_p, ellipse_points_t1(2,:)),'B-');
   
   % End test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

   % fill the grid_line field in structure domain
   domain = lsGenerateGridLines(domain, sc);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   % (find the known points)
   [x_X_i_t0, y_X_i_t0, x_Y_i_t0, y_Y_i_t0] = lsGetGridIntersections(x_spline_t0, y_spline_t0, domain, CONTROL);
   known_zero_level_points_t0(:,1) = [x_X_i_t0'; x_Y_i_t0'];
   known_zero_level_points_t0(:,2) = [y_X_i_t0'; y_Y_i_t0'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   [x_X_i_t1, y_X_i_t1, x_Y_i_t1, y_Y_i_t1] = lsGetGridIntersections(x_spline_t1, y_spline_t1, domain, CONTROL);
   known_zero_level_points_t1(:,1) = [x_X_i_t1'; x_Y_i_t1'];
   known_zero_level_points_t1(:,2) = [y_X_i_t1'; y_Y_i_t1'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif TEST_CASE == 2 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % two offset circles
   
   domain.x_size = 200;
   domain.y_size = 200;
   
   domain.x_spacing = x_s; 
   domain.y_spacing = y_s;
   
   %create circle
   circle1   = rsmak('circle',50,[0, 0]);
   circle2   = rsmak('circle',50,[0, 0]);   

   circle1   = fncmb(circle1,'+', 75);
   circle2   = fncmb(circle2,'+', 175);
   
   fnplt(circle1);
   hold on
   fnplt(circle2);
   axis equal
   
   % Create mask
   p = 0:0.05:circle1.pieces;
   circle1_points = fnval(circle1,p);
   p = 0:0.05:circle2.pieces;
   circle2_points = fnval(circle2,p);
   
   mask_img_t0 = roipoly(domain.y_size, domain.x_size, circle1_points(1,:)', circle1_points(2,:)');
   mask_img_t1 = roipoly(domain.y_size, domain.x_size, circle2_points(1,:)', circle2_points(2,:)');

   % create x,y splines
   s_p = 1:length(p);
   x_spline_t0 = fn2fm(spline(s_p, circle1_points(1,:)),'B-');
   y_spline_t0 = fn2fm(spline(s_p, circle1_points(2,:)),'B-');
   
   x_spline_t1 = fn2fm(spline(s_p, circle2_points(1,:)),'B-');
   y_spline_t1 = fn2fm(spline(s_p, circle2_points(2,:)),'B-');
   
   % End test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

    % fill the grid_line field in structure domain
    domain = lsGenerateGridLines(domain);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   % (find the known points)
   [x_X_i_t0, y_X_i_t0, x_Y_i_t0, y_Y_i_t0] = lsGetGridIntersections(x_spline_t0, y_spline_t0, domain, CONTROL);
   known_zero_level_points_t0(:,1) = [x_X_i_t0'; x_Y_i_t0'];
   known_zero_level_points_t0(:,2) = [y_X_i_t0'; y_Y_i_t0'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   [x_X_i_t1, y_X_i_t1, x_Y_i_t1, y_Y_i_t1] = lsGetGridIntersections(x_spline_t1, y_spline_t1, domain, CONTROL);
   known_zero_level_points_t1(:,1) = [x_X_i_t1'; x_Y_i_t1'];
   known_zero_level_points_t1(:,2) = [y_X_i_t1'; y_Y_i_t1'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
elseif TEST_CASE == 3 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  two lines
   domain.x_size = 200;
   domain.y_size = 200;
   
   domain.x_spacing = x_s; 
   domain.y_spacing = y_s;
   
   %create first line
   line1(:,1) = 1:200;
   line1(:,2) = 40 .* ones(1,200);  
   
   %create second line
   line2(:,1) = 1:200;
   line2(:,2) = 70 .* ones(1,200);   
   
   mask_img_t0 = roipoly(domain.y_size, domain.x_size, [0; line1(:,1); domain.x_size;0], [line1(1,2); line1(:,2); 0; 0]);
   mask_img_t1 = roipoly(domain.y_size, domain.x_size, [0; line2(:,1); domain.x_size;0], [line2(1,2); line2(:,2); 0; 0]);
   
   % create x,y interpolating splines
   s_p = 1:200;

   x_spline_t0 = fn2fm(spline(s_p, line1(:,1)),'B-');
   y_spline_t0 = fn2fm(spline(s_p, line1(:,2)),'B-');
   
   x_spline_t1 = fn2fm(spline(s_p, line2(:,1)),'B-');
   y_spline_t1 = fn2fm(spline(s_p, line2(:,2)),'B-');
   
   % End test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

   % fill the grid_line field in structure domain
   domain = lsGenerateGridLines(domain);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   % (find the known points)
   [x_X_i_t0, y_X_i_t0, x_Y_i_t0, y_Y_i_t0] = lsGetGridIntersections(x_spline_t0, y_spline_t0, domain, CONTROL);
   known_zero_level_points_t0(:,1) = [x_X_i_t0'; x_Y_i_t0'];
   known_zero_level_points_t0(:,2) = [y_X_i_t0'; y_Y_i_t0'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   [x_X_i_t1, y_X_i_t1, x_Y_i_t1, y_Y_i_t1] = lsGetGridIntersections(x_spline_t1, y_spline_t1, domain, CONTROL);
   known_zero_level_points_t1(:,1) = [x_X_i_t1'; x_Y_i_t1'];
   known_zero_level_points_t1(:,2) = [y_X_i_t1'; y_Y_i_t1'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
elseif TEST_CASE == 4
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  line protrusion
   domain.x_size = 200;
   domain.y_size = 200;
   
   domain.x_spacing = x_s; 
   domain.y_spacing = y_s;
   
   line_level1 = 105;
   line_level2 = 100;
   
   %create first line + halfcircle
   line1(:,1) = 1:200;   
   line1(1:44,2) = line_level1 .* ones(1,44);  
   line1(156:200,2) = line_level1 .* ones(1,45);
   xc=-55:55;
   %y_half_circle = line_level1 - sqrt(25^2-xc.^2);
   % gauss curve
   y_half_circle = line_level1 + 50 .* exp(-xc.^2./(2*160));
   line1(45:155,2)  = y_half_circle';
   
   %create second line
   line2(:,1) = 1:200;
   line2(:,2) = line_level2 .* ones(1,200);   
   
   mask_img_t1 = roipoly(domain.y_size, domain.x_size, [0; line1(:,1); domain.x_size;0], [line1(1,2); line1(:,2); 0; 0]);
   mask_img_t0 = roipoly(domain.y_size, domain.x_size, [0; line2(:,1); domain.x_size;0], [line2(1,2); line2(:,2); 0; 0]);
   
   % create x,y interpolating splines
   s_p = 1:200;

   x_spline_t1 = fn2fm(spline(s_p, line1(:,1)),'B-');
   y_spline_t1 = fn2fm(spline(s_p, line1(:,2)),'B-');
   
   x_spline_t0 = fn2fm(spline(s_p, line2(:,1)),'B-');
   y_spline_t0 = fn2fm(spline(s_p, line2(:,2)),'B-');
   
   % End test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

   % fill the grid_line field in structure domain
   domain = lsGenerateGridLines(domain);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   % (find the known points)
   [x_X_i_t0, y_X_i_t0, x_Y_i_t0, y_Y_i_t0] = lsGetGridIntersections(x_spline_t0, y_spline_t0, domain, CONTROL);
   known_zero_level_points_t0(:,1) = [x_X_i_t0'; x_Y_i_t0'];
   known_zero_level_points_t0(:,2) = [y_X_i_t0'; y_Y_i_t0'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   [x_X_i_t1, y_X_i_t1, x_Y_i_t1, y_Y_i_t1] = lsGetGridIntersections(x_spline_t1, y_spline_t1, domain, CONTROL);
   known_zero_level_points_t1(:,1) = [x_X_i_t1'; x_Y_i_t1'];
   known_zero_level_points_t1(:,2) = [y_X_i_t1'; y_Y_i_t1'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
elseif TEST_CASE == 5
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  lines and spike
   domain.x_size = 200;
   domain.y_size = 200;
   
   domain.x_spacing = x_s; 
   domain.y_spacing = y_s;
   
   line_level1 = 105.5;
   line_level2 = 110.5;
   
   %create first line + halfcircle
   line1(:,1) = 1:200;   
   line1(:,2) = line_level1 .* ones(1,200);  
   xc=0:80;
   %y_half_circle = line_level1 - sqrt(25^2-xc.^2);
   % gauss curve
   y_spike = line_level1 - exp(-xc./15) ./0.014;
   line1(20:100,2)  = fliplr(y_spike)';
   line1(101:180,2)  = y_spike(2:end)';
   figure,plot(line1(:,1),line1(:,2));
   
   %create second line
   line2(:,1) = 1:200;
   line2(:,2) = line_level2 .* ones(1,200);   
   
   mask_img_t0 = roipoly(domain.y_size, domain.x_size, [0; line1(:,1); domain.x_size;0], [line1(1,2); line1(:,2); 0; 0]);
   mask_img_t1 = roipoly(domain.y_size, domain.x_size, [0; line2(:,1); domain.x_size;0], [line2(1,2); line2(:,2); 0; 0]);
   
   % create x,y interpolating splines
   s_p = 1:200;

   x_spline_t0 = fn2fm(spline(s_p, line1(:,1)),'B-');
   y_spline_t0 = fn2fm(spline(s_p, line1(:,2)),'B-');
   
   x_spline_t1 = fn2fm(spline(s_p, line2(:,1)),'B-');
   y_spline_t1 = fn2fm(spline(s_p, line2(:,2)),'B-');
   
   % End test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

   % fill the grid_line field in structure domain
   domain = lsGenerateGridLines(domain);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   % (find the known points)
   [x_X_i_t0, y_X_i_t0, x_Y_i_t0, y_Y_i_t0] = lsGetGridIntersections(x_spline_t0, y_spline_t0, domain, CONTROL);
   known_zero_level_points_t0(:,1) = [x_X_i_t0'; x_Y_i_t0'];
   known_zero_level_points_t0(:,2) = [y_X_i_t0'; y_Y_i_t0'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   [x_X_i_t1, y_X_i_t1, x_Y_i_t1, y_Y_i_t1] = lsGetGridIntersections(x_spline_t1, y_spline_t1, domain, CONTROL);
   known_zero_level_points_t1(:,1) = [x_X_i_t1'; x_Y_i_t1'];
   known_zero_level_points_t1(:,2) = [y_X_i_t1'; y_Y_i_t1'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
   
   
elseif TEST_CASE == 6
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Seven-point star
   domain.x_size = 200;
   domain.y_size = 200;
   
   domain.x_spacing = x_s; 
   domain.y_spacing = y_s;
   
   %create star
   s_scale = 25;
   s=0:0.01:1;
   x=(100+s_scale.*((2+sin(7.*2.*pi.*s)).*cos(2.*pi.*s)));
   y=(100+s_scale.*((2+sin(7.*2.*pi.*s)).*sin(2.*pi.*s)));
   
   circle   = rsmak('circle',22,[0, 0]);
   circle   = fncmb(circle,'+', 100);
   
   figure
   fnplt(circle);
   hold on
   plot(x,y);
   axis equal
   
   % Create mask
   p = 0:0.05:circle.pieces;
   circle_points = fnval(circle,p);
   
   mask_img_t0 = roipoly(domain.y_size, domain.x_size, circle_points(1,:)', circle_points(2,:)');
   mask_img_t1 = roipoly(domain.y_size, domain.x_size, x, y);

   % create x,y splines
   s_p = 1:length(p);
   x_spline_t0 = fn2fm(spline(s_p, circle_points(1,:)),'B-');
   y_spline_t0 = fn2fm(spline(s_p, circle_points(2,:)),'B-');
     
   x_spline_t1 = fn2fm(spline(1:length(s), x),'B-');
   y_spline_t1 = fn2fm(spline(1:length(s), y),'B-');
   
   % End test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

    % fill the grid_line field in structure domain
    domain = lsGenerateGridLines(domain);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   % (find the known points)
   [x_X_i_t0, y_X_i_t0, x_Y_i_t0, y_Y_i_t0] = lsGetGridIntersections(x_spline_t0, y_spline_t0, domain, CONTROL);

   known_zero_level_points_t0(:,1) = [x_X_i_t0'; x_Y_i_t0'];
   known_zero_level_points_t0(:,2) = [y_X_i_t0'; y_Y_i_t0'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   [x_X_i_t1, y_X_i_t1, x_Y_i_t1, y_Y_i_t1] = lsGetGridIntersections(x_spline_t1, y_spline_t1, domain, CONTROL);
   known_zero_level_points_t1(:,1) = [x_X_i_t1'; x_Y_i_t1'];
   known_zero_level_points_t1(:,2) = [y_X_i_t1'; y_Y_i_t1'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif TEST_CASE == 7
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   domain.x_size = 439;
   domain.y_size = 345;
   
   domain.x_spacing = x_s; 
   domain.y_spacing = y_s;
   
   %cd L:\projects\rho_protrusion\cell1\protrusion
   cd /lccb/projects/rho_protrusion/cell1/protrusion_1-30_ps2_s40
   load edge_spline
   
   PROJECT_DIR = '/lccb/projects/rho_protrusion/cell1/';
   PROT_DIR = 'protrusion_1-30_ps2_s40/';
   IMG_NAME = 'plane01';
   
   % file name containing the edge pixels
   file_pixel_edge=[PROJECT_DIR  PROT_DIR 'pixel_edge.dat'];
   
   fid_pixel_edge= fopen(file_pixel_edge,'r');
   if fid_pixel_edge == -1
      error('Could not open file pixel edge');
   end
   
   % mask file names
   firstfilename_mask=([PROJECT_DIR PROT_DIR 'mask_' IMG_NAME '.tif']);
   [filelist_mask]=getFileStackNames(firstfilename_mask);
   
   time = 4;
   time_increment = 10;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the pixel edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i = 1 : time
       n_pix_t0     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       y_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       known_zero_level_points_t0 = [x_p_edge_t0', y_p_edge_t0'];
   end
   for i = 1 : time_increment
       n_pix_t1     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       y_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       known_zero_level_points_t1 = [x_p_edge_t1', y_p_edge_t1'];
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   fileName_mask=char(filelist_mask(time));
   mask_img_t0=imread(fileName_mask);
   
   fileName_mask=char(filelist_mask(time+time_increment));
   mask_img_t1=imread(fileName_mask);   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % End data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   x_spline_t0 = edge_sp_array_x(time);
   y_spline_t0 = edge_sp_array_y(time);
  
   x_spline_tb = edge_sp_array_x(time+round(time_increment/2));
   y_spline_tb = edge_sp_array_y(time+round(time_increment/2));
   
   x_spline_t1 = edge_sp_array_x(time+time_increment);
   y_spline_t1 = edge_sp_array_y(time+time_increment);
   
   
    [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

    % fill the grid_line field in structure domain
    domain = lsGenerateGridLines(domain);
    
    
    
elseif TEST_CASE == 8
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   domain.x_size = 439;
   domain.y_size = 345;
   
   domain.x_spacing = x_s; 
   domain.y_spacing = y_s;
   
   cd /lccb/projects/alpha/W512r_smaller/protrusion_01-90_s30_p20
   load edge_spline
   
   PROJECT_DIR = '/lccb/projects/alpha/W512r_smaller/';
   PROT_DIR = 'protrusion_01-90_s30_p20/';
   IMG_NAME = 'cut_cut_cut_W512r01';
   
   % file name containing the edge pixels
   file_pixel_edge=[PROJECT_DIR  PROT_DIR 'pixel_edge.dat'];
   
   fid_pixel_edge= fopen(file_pixel_edge,'r');
   if fid_pixel_edge == -1
      error('Could not open file pixel edge');
   end
   
   % mask file names
   firstfilename_mask=([PROJECT_DIR PROT_DIR 'mask_' IMG_NAME '.tif']);
   [filelist_mask]=getFileStackNames(firstfilename_mask);
   
   time = 1;
   time_increment = 1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the pixel edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i = 1 : time
       n_pix_t0     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       y_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       known_zero_level_points_t0 = [x_p_edge_t0', y_p_edge_t0'];
   end
   for i = 1 : time_increment
       n_pix_t1     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       y_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       known_zero_level_points_t1 = [x_p_edge_t1', y_p_edge_t1'];
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   fileName_mask=char(filelist_mask(time));
   mask_img_t0=imread(fileName_mask);
   
   fileName_mask=char(filelist_mask(time+time_increment));
   mask_img_t1=imread(fileName_mask);   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % End data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   x_spline_t0 = edge_sp_array_x(time);
   y_spline_t0 = edge_sp_array_y(time);
  
   x_spline_tb = edge_sp_array_x(time+round(time_increment/2));
   y_spline_tb = edge_sp_array_y(time+round(time_increment/2));
   
   x_spline_t1 = edge_sp_array_x(time+time_increment);
   y_spline_t1 = edge_sp_array_y(time+time_increment);
   
   
   [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

   % fill the grid_line field in structure domain
   domain = lsGenerateGridLines(domain);
    
elseif TEST_CASE == 9
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   cd /lccb/projects/alpha/PtK1_control/cut_s399/protrusion_1-243_s40_p20
   load edge_spline
   
   PROJECT_DIR = '/lccb/projects/alpha/PtK1_control/cut_s399/';
   PROT_DIR = 'protrusion_1-243_s40_p20/';
   IMG_NAME = 'cut_S399ACTIN001';
   
   % Read image to determine the size
   c_image = imread('img_edge_cut_S399ACTIN001.tif');
   [img_h, img_w, c] = size(c_image);
   
   domain.x_size = img_w;
   domain.y_size = img_h;
   
   domain.x_spacing = x_s; 
   domain.y_spacing = y_s;
   
   
   % file name containing the edge pixels
   file_pixel_edge=[PROJECT_DIR  PROT_DIR 'pixel_edge.dat'];
   
   fid_pixel_edge= fopen(file_pixel_edge,'r');
   if fid_pixel_edge == -1
      error('Could not open file pixel edge');
   end
   
   % mask file names
   firstfilename_mask=([PROJECT_DIR PROT_DIR 'mask_' IMG_NAME '.tif']);
   [filelist_mask]=getFileStackNames(firstfilename_mask);
   
   time = 1;
   time_increment = 1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the pixel edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i = 1 : time
       n_pix_t0     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       y_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       known_zero_level_points_t0 = [x_p_edge_t0', y_p_edge_t0'];
   end
   for i = 1 : time_increment
       n_pix_t1     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       y_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       known_zero_level_points_t1 = [x_p_edge_t1', y_p_edge_t1'];
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   fileName_mask=char(filelist_mask(time));
   mask_img_t0=imread(fileName_mask);
   
   fileName_mask=char(filelist_mask(time+time_increment));
   mask_img_t1=imread(fileName_mask);   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % End data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   x_spline_t0 = edge_sp_array_x(time);
   y_spline_t0 = edge_sp_array_y(time);
  
   x_spline_tb = edge_sp_array_x(time+round(time_increment/2));
   y_spline_tb = edge_sp_array_y(time+round(time_increment/2));
   
   x_spline_t1 = edge_sp_array_x(time+time_increment);
   y_spline_t1 = edge_sp_array_y(time+time_increment);
   
   
    [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

    % fill the grid_line field in structure domain
    domain = lsGenerateGridLines(domain);  
end

