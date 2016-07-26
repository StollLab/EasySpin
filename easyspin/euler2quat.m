%  euler2quat  Convert set of Euler angles to a quaternion.
%
%  q = euler2quat(alpha, beta, gamma);
%
%  Input:
%     alpha          double or 1xN array, first Euler angle
%     beta           double or 1xN array, second Euler angle
%     gamma          double or 1xN array, third Euler angle
%
%  Output:
%     q              4xN array, normalized quaternion

function q = euler2quat(alpha, beta, gamma)

    q0 = cos(beta/2).*cos((gamma+alpha)/2);
    q1 = sin(beta/2).*sin((gamma-alpha)/2);
    q2 = sin(beta/2).*cos((gamma-alpha)/2);
    q3 = cos(beta/2).*sin((gamma+alpha)/2);
    
    q = [q0; q1; q2; q3];

end