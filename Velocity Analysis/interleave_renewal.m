function out = interleave_renewal(data1, data2)
    out = zeros(size(data1,1), size(data1,2), size(data1,3) + size(data2,3));
    out(:, :, 2:2:end) = data1;
    out(:, :, 1:2:end) = data2;
end
