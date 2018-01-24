function [x_X_i, y_X_i, x_Y_i, y_Y_i] = lsGetGridIntersections(x_spline, y_spline, domain, CONTROL)

if 0
   % Get the intersection of the x-gridlines 
   [x_x_intersection, x_y_intersection]= lsIntersectMaskLine(domain.x_grid_lines, mask_img, 1);

   % Improve the localization of the intersections
   [x_x_intersection, x_y_intersection]= lsIntersectSplineLine(domain.x_grid_lines, mask_img, 1);
   
   % Get the intersection of the y-gridlines 
   [y_x_intersection, y_y_intersection] = lsIntersectMaskLine(domain.y_grid_lines, mask_img, 2);


   % all intersection points
   x_intersection = cat(2, x_x_intersection, y_x_intersection);
   y_intersection = cat(2, x_y_intersection, y_y_intersection);

   figure,plot(x_intersection, y_intersection,'.');
   hold on
   p=1:x_spline.knots(end);
   x_spline_points = fnval(x_spline, p);
   y_spline_points = fnval(y_spline, p);
   plot(x_spline_points, y_spline_points,'r');
end


if 1
   % get the intersection of the x-gridlines with the edge
   x1(1:length(domain.x_grid_lines)) = x_spline.knots(1);
   x2(1:length(domain.x_grid_lines)) = x_spline.knots(end);
   [x_approx_intersection_parameter, x_intersecting_lines, ds] = lsIntersectApproxSplineLine(domain.x_grid_lines, x_spline, x1, x2);
   
   x1 = x_approx_intersection_parameter - ds;
   x2 = x_approx_intersection_parameter + ds;
   x_intersection_parameter = lsIntersectSplineLine(x_intersecting_lines, x_spline, x1, x2);

   % get the intersection of the y-gridlines with the edge
   y1(1:length(domain.y_grid_lines)) = y_spline.knots(1);
   y2(1:length(domain.y_grid_lines)) = y_spline.knots(end);
   [y_approx_intersection_parameter, y_intersecting_lines, ds] = lsIntersectApproxSplineLine(domain.y_grid_lines, y_spline, y1, y2);
   
   y1 = y_approx_intersection_parameter - ds;
   y2 = y_approx_intersection_parameter + ds;
   y_intersection_parameter = lsIntersectSplineLine(y_intersecting_lines, y_spline, y1, y2);

   x_X_i = fnval(x_spline,x_intersection_parameter);
   y_X_i = fnval(y_spline,x_intersection_parameter);

   x_Y_i = fnval(x_spline,y_intersection_parameter);
   y_Y_i = fnval(y_spline,y_intersection_parameter);

   % all intersection points
   x_intersection = cat(2, x_X_i, x_Y_i);
   y_intersection = cat(2, y_X_i, y_Y_i);
   
   if CONTROL
       figure,plot(x_X_i, y_X_i,'.');
       hold on
       plot(x_Y_i,y_Y_i,'.g');
       p=x_spline.knots(1) :ds :x_spline.knots(end);
       x_spline_points = fnval(x_spline, p);
       y_spline_points = fnval(y_spline, p);
       plot(x_spline_points, y_spline_points,'r');
   end
end

if CONTROL
    % plot the grid lines
    for i=1:length(domain.x_grid_lines)
        line([domain.x_grid_lines(i) domain.x_grid_lines(i)],[0 domain.y_size]);
    end
    for i=1:length(domain.y_grid_lines)
        line([0 domain.x_size], [domain.y_grid_lines(i) domain.y_grid_lines(i)],'Color','c');
    end
end
