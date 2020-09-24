function circles = find_circles_ransac(x,y)
% parameters
sample_size = 4;
max_distance = 50;
min_num_points_in_circle = 50;

circles = {};
counter = 0;
while numel(x) > 100 * sample_size && counter < 10
    try
        [circle, inlierIdx] = ransac([x, y], @fit_circle, @distFcn, ...
            sample_size, max_distance);
    catch
        break
    end
    % refit using only inliers
    circle = fit_circle([x(inlierIdx) y(inlierIdx)]);
    distf = distFcn(circle, [x y]);
    founfit = distf < max_distance;
    if sum(founfit) > min_num_points_in_circle
        % this model fits enough points
        circles{end+1} = circle;
        x(founfit) = [];
        y(founfit) = [];
    else
        counter = counter + 1;
    end
end
circles = vertcat(circles{:});
end

function crc = fit_circle(xy)  % xy is n-by-2 matrix of point coordinates
% fit in least squares sens
x = xy(:, 1);
y = xy(:, 2);
X = [-2*x, -2*y, ones(size(x))];
Y = -x.^2 - y.^2;
crc = (X\Y).';  % least squares solution
r2 = -crc(3) +crc(1).^2 + crc(2).^2;
if r2 <= 0
    crc(3) = 0;
else
    crc(3) = sqrt(r2);
end
% output crc is a 3 vector (cx, cy, r)
end



function distf = distFcn(crc, xy)
% how good a fit circle for points
x = xy(:, 1) - crc(1);
y = xy(:, 2) - crc(2);
distf = abs(sqrt(x.^2 + y.^2) - crc(3));
end

