function [lat, lon] = fibonacci_sphere(samples)
    % Generate points on a sphere using the Fibonacci method
    % Output in degrees: latitude (-90 to 90) and longitude (-180 to 180)

    if nargin < 1
        samples = 1000; % default
    end

    phi = pi * (sqrt(5) - 1); % golden angle in radians

    lat = zeros(samples, 1);
    lon = zeros(samples, 1);

    for i = 0:samples-1
        y = 1 - (i / (samples - 1)) * 2; % y from 1 to -1
        radius = sqrt(1 - y^2);          % radius at y

        theta = phi * i;

        x = cos(theta) * radius;
        z = sin(theta) * radius;

        % Convert to latitude/longitude
        lat(i+1) = asind(y);                   % latitude in degrees
        lon(i+1) = atan2d(z, x);               % longitude in degrees
    end
end