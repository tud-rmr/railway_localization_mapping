function nu = simFcn_Imm_ResiduumFcn(z,z_hat)

if size(z,1) > 2
    
    nu = zeros(size(z));
    nu(1) = z(1)-z_hat(1);
    nu(2) = z(2)-z_hat(2);
    nu(3) = z(3)-z_hat(3);
    nu(4) = -angdiff(z(4),z_hat(4));
    nu(5) = z(5)-z_hat(5);
    nu(6) = z(6)-z_hat(6);
    
else
    
    nu = z-z_hat;
    
end % if
    
% nu = zeros(size(z));
% nu(1) = z(1)-z_hat(1);
% nu(2) = z(2)-z_hat(2);
% nu(3) = z(3)-z_hat(3);
% nu(4) = -angdiff(z(4),z_hat(4));

end % end function



