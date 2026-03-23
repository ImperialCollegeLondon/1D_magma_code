function coeffs = piecewiseFit(data)
% Returns explicit polynomial coefficients:
%
% 1D:  f(x) = a*x + b
% 2D:  f(x,y) = a*x + b*y + c*x*y + d
% 3D:  f(x,y,z) = a*x + b*y + c*z + d*x*y + e*x*z + f*y*z + g*x*y*z + h

dims = ndims(data);

if isvector(data)
    %% 1D
    data = data(:);
    N = length(data);
    x = linspace(0,1,N);

    a = zeros(N-1,1);
    b = zeros(N-1,1);

    for i = 1:N-1
        x1 = x(i); x2 = x(i+1);
        f1 = data(i); f2 = data(i+1);

        dx = x2 - x1;

        a(i) = (f2 - f1)/dx;
        b(i) = f1 - a(i)*x1;
    end

    coeffs=[a b];

elseif dims == 2
    %% 2D
    [N,M] = size(data);
    x = linspace(0,1,N);
    y = linspace(0,1,M);

    a = zeros(N-1,M-1);
    b = zeros(N-1,M-1);
    c = zeros(N-1,M-1);
    d = zeros(N-1,M-1);

    for i = 1:N-1
        for j = 1:M-1

            x1 = x(i);   x2 = x(i+1);
            y1 = y(j);   y2 = y(j+1);

            z11 = data(i,j);
            z21 = data(i+1,j);
            z12 = data(i,j+1);
            z22 = data(i+1,j+1);

            dx = x2 - x1;
            dy = y2 - y1;

            % === CORRECT COEFFICIENTS ===
            c(i,j) = (z22 - z21 - z12 + z11)/(dx*dy);

            a(i,j) = (z21 - z11)/dx - c(i,j)*y1;
            b(i,j) = (z12 - z11)/dy - c(i,j)*x1;

            d(i,j) = z11 - a(i,j)*x1 - b(i,j)*y1 - c(i,j)*x1*y1;

        end
    end
    

    coeffs=cat(3,a,b,c,d);

elseif dims == 3
    %% 3D
    [N,M,P] = size(data);

    x = linspace(0,1,N);
    y = linspace(0,1,M);
    z = linspace(0,1,P);

    % 8 coefficients
    A = zeros(N-1,M-1,P-1,8);

    for i = 1:N-1
        for j = 1:M-1
            for k = 1:P-1

                x1 = x(i); x2 = x(i+1);
                y1 = y(j); y2 = y(j+1);
                z1 = z(k); z2 = z(k+1);

                dx = x2 - x1;
                dy = y2 - y1;
                dz = z2 - z1;

                f000 = data(i,j,k);
                f100 = data(i+1,j,k);
                f010 = data(i,j+1,k);
                f110 = data(i+1,j+1,k);
                f001 = data(i,j,k+1);
                f101 = data(i+1,j,k+1);
                f011 = data(i,j+1,k+1);
                f111 = data(i+1,j+1,k+1);

                % highest-order term
                g = (f111 - f110 - f101 - f011 + f100 + f010 + f001 - f000)/(dx*dy*dz);

                % mixed terms
                dxy = (f110 - f100 - f010 + f000)/(dx*dy) - g*z1;
                dxz = (f101 - f100 - f001 + f000)/(dx*dz) - g*y1;
                dyz = (f011 - f010 - f001 + f000)/(dy*dz) - g*x1;

                % linear terms
                ax = (f100 - f000)/dx - dxy*y1 - dxz*z1 - g*y1*z1;
                by = (f010 - f000)/dy - dxy*x1 - dyz*z1 - g*x1*z1;
                cz = (f001 - f000)/dz - dxz*x1 - dyz*y1 - g*x1*y1;

                % constant
                h = f000 - ax*x1 - by*y1 - cz*z1 ...
                          - dxy*x1*y1 - dxz*x1*z1 - dyz*y1*z1 ...
                          - g*x1*y1*z1;

                A(i,j,k,:) = [ax by cz dxy dxz dyz g h];

            end
        end
    end
    
    coeffs=A;
elseif dims==4
    [N,M,P,Q] = size(data);

    x = linspace(0,1,N);
    y = linspace(0,1,M);
    z = linspace(0,1,P);
    w = linspace(0,1,Q);

    A = zeros(N-1,M-1,P-1,Q-1,16);

    for i = 1:N-1
        for j = 1:M-1
            for k = 1:P-1
                for l = 1:Q-1

                    x1 = x(i); x2 = x(i+1);
                    y1 = y(j); y2 = y(j+1);
                    z1 = z(k); z2 = z(k+1);
                    w1 = w(l); w2 = w(l+1);

                    dx = x2 - x1;
                    dy = y2 - y1;
                    dz = z2 - z1;
                    dw = w2 - w1;

                    % === 16 corner values ===
                    f0000 = data(i,j,k,l);
                    f1000 = data(i+1,j,k,l);
                    f0100 = data(i,j+1,k,l);
                    f1100 = data(i+1,j+1,k,l);

                    f0010 = data(i,j,k+1,l);
                    f1010 = data(i+1,j,k+1,l);
                    f0110 = data(i,j+1,k+1,l);
                    f1110 = data(i+1,j+1,k+1,l);

                    f0001 = data(i,j,k,l+1);
                    f1001 = data(i+1,j,k,l+1);
                    f0101 = data(i,j+1,k,l+1);
                    f1101 = data(i+1,j+1,k,l+1);

                    f0011 = data(i,j,k+1,l+1);
                    f1011 = data(i+1,j,k+1,l+1);
                    f0111 = data(i,j+1,k+1,l+1);
                    f1111 = data(i+1,j+1,k+1,l+1);

                    % === Highest order ===
                    a15 = (f1111 - f1110 - f1101 - f1011 - f0111 ...
                        + f1100 + f1010 + f0110 + f1001 + f0101 + f0011 ...
                        - f1000 - f0100 - f0010 - f0001 + f0000) / (dx*dy*dz*dw);

                    % === 3-variable terms ===
                    a11 = (f1110 - f1100 - f1010 - f0110 + f1000 + f0100 + f0010 - f0000)/(dx*dy*dz) - a15*w1;
                    a12 = (f1101 - f1100 - f1001 - f0101 + f1000 + f0100 + f0001 - f0000)/(dx*dy*dw) - a15*z1;
                    a13 = (f1011 - f1010 - f1001 - f0011 + f1000 + f0010 + f0001 - f0000)/(dx*dz*dw) - a15*y1;
                    a14 = (f0111 - f0110 - f0101 - f0011 + f0100 + f0010 + f0001 - f0000)/(dy*dz*dw) - a15*x1;

                    % === 2-variable terms ===
                    a5 = (f1100 - f1000 - f0100 + f0000)/(dx*dy) - a11*z1 - a12*w1 - a15*z1*w1;
                    a6 = (f1010 - f1000 - f0010 + f0000)/(dx*dz) - a11*y1 - a13*w1 - a15*y1*w1;
                    a7 = (f1001 - f1000 - f0001 + f0000)/(dx*dw) - a12*y1 - a13*z1 - a15*y1*z1;

                    a8 = (f0110 - f0100 - f0010 + f0000)/(dy*dz) - a11*x1 - a14*w1 - a15*x1*w1;
                    a9 = (f0101 - f0100 - f0001 + f0000)/(dy*dw) - a12*x1 - a14*z1 - a15*x1*z1;
                    a10 = (f0011 - f0010 - f0001 + f0000)/(dz*dw) - a13*x1 - a14*y1 - a15*x1*y1;

                    % === linear terms ===
                    a1 = (f1000 - f0000)/dx - a5*y1 - a6*z1 - a7*w1 ...
                        - a11*y1*z1 - a12*y1*w1 - a13*z1*w1 - a15*y1*z1*w1;

                    a2 = (f0100 - f0000)/dy - a5*x1 - a8*z1 - a9*w1 ...
                        - a11*x1*z1 - a12*x1*w1 - a14*z1*w1 - a15*x1*z1*w1;

                    a3 = (f0010 - f0000)/dz - a6*x1 - a8*y1 - a10*w1 ...
                        - a11*x1*y1 - a13*x1*w1 - a14*y1*w1 - a15*x1*y1*w1;

                    a4 = (f0001 - f0000)/dw - a7*x1 - a9*y1 - a10*z1 ...
                        - a12*x1*y1 - a13*x1*z1 - a14*y1*z1 - a15*x1*y1*z1;

                    % === constant ===
                    a16 = f0000 ...
                        - a1*x1 - a2*y1 - a3*z1 - a4*w1 ...
                        - a5*x1*y1 - a6*x1*z1 - a7*x1*w1 ...
                        - a8*y1*z1 - a9*y1*w1 - a10*z1*w1 ...
                        - a11*x1*y1*z1 - a12*x1*y1*w1 ...
                        - a13*x1*z1*w1 - a14*y1*z1*w1 ...
                        - a15*x1*y1*z1*w1;

                    A(i,j,k,l,:) = [a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 ...
                        a11 a12 a13 a14 a15 a16];

                end
            end
        end
    end

    coeffs=A;
else
    error('Input must be 1-4D array.');
end

end