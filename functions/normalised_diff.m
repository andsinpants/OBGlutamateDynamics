function norm_value = normalised_diff( data )
% Normalise values of an array to be between -1 and 1
% original sign of the array values is maintained.
if abs(min(data)) > max(data)
      max_range_value = abs(min(data));
      min_range_value = min(data);
  else
      max_range_value = max(data);
      min_range_value = -max(data);
  end
norm_value = 2 .* data ./ (max_range_value - min_range_value);
% norm_value2 = 2 .* data2 ./ (max_range_value - min_range_value);
end