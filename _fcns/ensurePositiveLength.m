function x_corrected = ensurePositiveLength(x,x_descriptor)
% x_corrected = ensurePositiveLength(x,x_descriptor)
% 

x_descriptor = x_descriptor(x_descriptor~=11);

arc_selector = ( ... 
                 (x_descriptor ==   3) | ... 
                 (x_descriptor ==  31) | ...
                 (x_descriptor ==  13) | ... 
                 (x_descriptor == 131)   ... 
               );
negative_x_selector = (x < 0);

negative_length_selector = negative_x_selector;
negative_length_selector(:,~arc_selector) = false;
negative_length_selector(2,:) = false;

x_corrected = x;
if any(negative_length_selector(:))
    x_corrected(negative_length_selector) = abs(x(negative_length_selector));
end % if

end

