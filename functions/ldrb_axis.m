function Q = ldrb_axis(gradPsi, gradPhi, tol)
% Function 2 in supplement of Bayer 2012:
% https://doi.org/10.1007/s10439-012-0593-5

% gradPsi: gradient in apicobasal direction
% gradPhi: gradient in transmural direction

if(norm(gradPsi) < tol)
    gradPsi = [1;0;0];
end
e1 = gradPsi/norm(gradPsi);

if(norm(gradPhi) < tol)
    gradPhi = [0;0;1];
end
e2 = gradPhi-e1'*gradPhi*e1;
if(norm(e2) < tol)
    e2 = [0;1;0];
end
e2 = e2/norm(e2);
e0 = cross(e1,e2);

Q = [e0 e1 e2];

end
