function coeffs = compute4Dcoeffs(data)
% compute4Dcoeffs - compute 16 multilinear coefficients for each 4D cell
%
% data: N x M x P x Q array of function values at grid points
% Coordinates are assumed to be linspace(0,1,...) along each dimension.
%
% coeffs.A(i,j,k,l,16): coefficients for cell (i,j,k,l)

[N,M,P,Q] = size(data);

x = linspace(0,1,N);
y = linspace(0,1,M);
z = linspace(0,1,P);
w = linspace(0,1,Q);

A = zeros(N-1, M-1, P-1, Q-1, 16);

% Precompute the 16×16 basis inverse (same for every cell)
B = buildBasisMatrix();
Binv = inv(B);

for i = 1:N-1
    for j = 1:M-1
        for k = 1:P-1
            for l = 1:Q-1

                % Corner coordinates
                X = [x(i), x(i+1)];
                Y = [y(j), y(j+1)];
                Z = [z(k), z(k+1)];
                W = [w(l), w(l+1)];

                % Corner values vector (in the matching order)
                f = zeros(16,1);
                idx = 1;
                for ix = 1:2
                    for iy = 1:2
                        for iz = 1:2
                            for iw = 1:2
                                f(idx) = data(i+ix-1, j+iy-1, k+iz-1, l+iw-1);
                                idx = idx + 1;
                            end
                        end
                    end
                end

                % Build evaluation matrix for this cell (16 rows)
                C = cellBasisPoints(X, Y, Z, W);

                % Solve for coefficients: C * coeff = f  →  coeff = C^{-1} f
                % but C = B evaluated at the corner coordinates
                coeff = (C \ f);

                A(i,j,k,l,:) = coeff;

            end
        end
    end
end

coeffs = A;

end

% ------------------------------------------------------------
% Build the universal 16×16 polynomial basis matrix
% ------------------------------------------------------------
function B = buildBasisMatrix()
B = zeros(16,16);

% basis: [x,y,z,w, xy,xz,xw,yz,yw,zw, xyz,xyw,xzw,yzw, xyzw, 1]

% Use generic symbolic corner values (0 or 1)
corners = dec2bin(0:15)-'0';

for r = 1:16
    x = corners(r,1);
    y = corners(r,2);
    z = corners(r,3);
    w = corners(r,4);

    B(r,:) = [
        x, y, z, w, ...
        x*y, x*z, x*w, y*z, y*w, z*w, ...
        x*y*z, x*y*w, x*z*w, y*z*w, ...
        x*y*z*w, ...
        1
    ];
end

end

% ------------------------------------------------------------
% Build the 16×16 matrix evaluated at real corner coordinates
% ------------------------------------------------------------
function C = cellBasisPoints(X,Y,Z,W)

C = zeros(16,16);

corners = dec2bin(0:15)-'0';

for r = 1:16
    ix = corners(r,1)+1;
    iy = corners(r,2)+1;
    iz = corners(r,3)+1;
    iw = corners(r,4)+1;

    x = X(ix);
    y = Y(iy);
    z = Z(iz);
    w = W(iw);

    C(r,:) = [
        x, y, z, w, ...
        x*y, x*z, x*w, y*z, y*w, z*w, ...
        x*y*z, x*y*w, x*z*w, y*z*w, ...
        x*y*z*w, ...
        1
    ];
end

end
