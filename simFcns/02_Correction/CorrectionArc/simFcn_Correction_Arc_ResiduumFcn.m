function nu = simFcn_Correction_Arc_ResiduumFcn(z,z_hat)
    
nu = zeros(size(z));
nu(1) = z(1)-z_hat(1);
nu(2) = z(2)-z_hat(2);
nu(3) = z(3)-z_hat(3);
% nu(4) = -angdiff(z(4),z_hat(4));

end % end function



