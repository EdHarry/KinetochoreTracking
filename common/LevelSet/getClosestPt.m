function c_point = getClosestPt(point, sp_x, sp_y)

SOLUTION_SPACING = 2;

if 1    
    r1 = sp_x.breaks(1);
    r2 = sp_x.breaks(end);
    x = fminbnd(@d_fun, r1, r2, [], sp_x, sp_y, point(1), point(2));
     c_point = [ppval(sp_x,x); ppval(sp_y,x)]';
else
    
    % generate points
    r_n_lower = sp_x.breaks(1);
    r_n_upper = sp_x.breaks(end);
    r_n  = r_n_lower: SOLUTION_SPACING : r_n_upper;


    %curve_points = [fnval(sp_x,r_n); fnval(sp_y,r_n)]';
    curve_points = [ppval(sp_x,r_n); ppval(sp_y,r_n)]';

    % get distances
    dist_cand = sqrt((curve_points(:,1) - point(1)).^2+ (curve_points(:,2) - point(2)).^2);

    %get the minimal distance
    [dist min_index] = min(dist_cand);

    %get the point
    c_point = curve_points(min_index,:);

    % get the spline parameter of this point
    r_n_min_estimate = r_n(min_index);

    % improve solution
    x = fminbnd(@d_fun, r_n_min_estimate - 5, r_n_min_estimate + 5, [], sp_x, sp_y, point(1), point(2));

    c_point = [ppval(sp_x,x); ppval(sp_y,x)]';
end

function d = d_fun(x, sp_x, sp_y, p1, p2)
c_p = [ppval(sp_x,x); ppval(sp_y,x)]';
d = sqrt((c_p(1) - p1).^2+ (c_p(2) - p2).^2);