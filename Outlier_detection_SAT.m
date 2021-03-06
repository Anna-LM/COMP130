function [residual_data, residual_data2, outlier] = Outlier_detection_SAT(residual_data, residual_data2, innovation_vectors, measurement_matrices, errorsd, det_threshold)
    
    
    if nargin < 6
        det_threshold = 6;
        errorsd = 5;
    end
    
    times = residual_data(2:size(residual_data,1),1);%first column
    
    sateillites = transpose(residual_data(1,2:size(residual_data,2)));%first row
    %initialise covariance matrix and residal vectors
    % Workshop 1 task 3
    residual_vectors = zeros(size(times,1),size(sateillites,1));
    residual_cov_mat = zeros(size(times,1),size(sateillites,1));

    for i = 1:size(times,1)
        dz = innovation_vectors(:,i);
        H_g = measurement_matrices(:,:,i);
        

        S = H_g*inv(H_g'*H_g)*H_g';

        Iden_m = eye(size(S));
        v = (S - Iden_m) * dz;
        %must use diag to avoid error in matrix multiplication
        C_v = diag((Iden_m - S)*errorsd^2);

        %store vectors and covariance matrix

        residual_vectors(i,:) = v';
        residual_cov_mat(i,:) = C_v';
    end
    %Workshop 1 task 3 part c - Compute normalised residuals and compare
    outliers_matrix = abs(residual_vectors) > residual_cov_mat.^0.5 * det_threshold;
    
    [R,C] = find(outliers_matrix == 1);
    lin_idx = sub2ind(size(residual_vectors),R,C);
    outliers = abs(residual_vectors(lin_idx));
    hi = find(outliers == max(outliers));
    [R,C] = ind2sub(size(residual_data), lin_idx(hi));

    residual_data(:,C+1) = [];
    residual_data2(:,C+1) = [];
    
    outlier = sum(outliers_matrix(:));
end