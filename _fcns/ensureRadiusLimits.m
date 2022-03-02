function x_corrected = ensureRadiusLimits(x,x_descriptor)
% x_corrected = ensureRadiusLimits(x,x_descriptor)
% 

x_descriptor = x_descriptor(x_descriptor~=11);

max_radius = 30e3;

arc_indices = find( ... 
                    (x_descriptor ==   3) | ... 
                    (x_descriptor ==  31) | ...
                    (x_descriptor ==  13) | ... 
                    (x_descriptor == 131)   ... 
                  );
r_too_big_indices = arc_indices(abs(x(2,arc_indices)) > max_radius);

x_corrected = x;
x_corrected(2,r_too_big_indices) = sign(x_corrected(2,r_too_big_indices))*max_radius;

end

