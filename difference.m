function [true_value] = difference(new_y_1, old_y_1, error_threshold)
    diff = abs(new_y_1 - old_y_1);
    
    % Determine base for relative difference
    base = max(abs(new_y_1), abs(old_y_1));
    
    % Handle division by zero
    if base == 0
        relative = 0;
    else
        relative = (diff/base)*100;
    end
    
    % Return the value if it meets the error criterion, otherwise NaN
    if relative <= error_threshold
        true_value = new_y_1;
    else
        true_value = NaN;
    end
end