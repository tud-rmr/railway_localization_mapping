function nu = simFcn_Correction_Straight_ResiduumFcn(z,z_hat)
    
nu = zeros(size(z));
nu(1) = z(1)-z_hat(1);
nu(2) = z(2)-z_hat(2);
% nu(3) = -angdiff(z(3),z_hat(3));

end % end function



