function norm_rows = normrois(data)
%normalises rois in mapanalysis output file to each row (roi)
norm_rows = zeros(size(data));
for i=1:size(data)
    norm_rows(i,:) = normalised_diff(data(i,:));
end
