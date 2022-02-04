function [data, data2, outlier_present] = REMOVE_OUTLIER_SATELLITES...
    (data, data2, innovation_vectors, measurement_matrices, error_SD, threshold)
    
    if nargin < 5
        error_SD = 5;
        threshold = 6;
    elseif nargin < 6
        threshold = 6;
    end
    
    t = data(2:size(data,1),1);
    sat = transpose(data(1,2:size(data,2)));
    
    residuals_vectors = zeros(size(t,1),size(sat,1));
    residuals_covvecs = zeros(size(t,1),size(sat,1));

    for i = 1:size(t,1)
        dz = innovation_vectors(:,i);
        H = measurement_matrices(:,:,i);

        A = H*inv(H'*H)*H';
        v = (A - eye(size(A))) * dz;
        Cv = diag((eye(size(A)) - A) * error_SD^2);

        residuals_vectors(i,:) = v';
        residuals_covvecs(i,:) = Cv';
    end

    outliers_matrix = abs(residuals_vectors) > residuals_covvecs.^0.5 * threshold;
    [R,C] = find(outliers_matrix == 1);
    lin_idx = sub2ind(size(residuals_vectors),R,C);
    outliers = abs(residuals_vectors(lin_idx));
    hi = find(outliers == max(outliers));
    [R,C] = ind2sub(size(data), lin_idx(hi));

    data(:,C+1) = [];
    data2(:,C+1) = [];
    
    outlier_present = sum(outliers_matrix(:));
end