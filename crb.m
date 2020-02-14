function bound = crb(Nt, K, M, theta, theta_d, phi, phi_d, SNR, u)
sigma_s = 10 ^ (SNR * 0.1);%power of signal
sigma_n = 1;%power of noise
Mx = sqrt(M);
My = Mx;
R = sigma_n * eye(M);
a = zeros(M, 1);
D_store = zeros(M * K, M);
B_store = zeros(M * K, M);
theta = theta / 180 * pi;
theta_d = theta_d / 180 * pi;
phi = phi / 180 * pi;
phi_d = phi_d / 180 * pi;
for k = 1 : K
    for mx = 1 : Mx
        for my = 1 : My
            m = mx + (my - 1) * Mx;
            a(m, 1) = exp(1i * u * sin(phi(k, 1)) * ((mx - 1) * cos(theta(k, 1)) + (my - 1) * sin(theta(k, 1))));
        end
    end
    D = diag(a);
    D_store((k - 1) * M + 1 : k * M, :) = D;
    B = eye(M);
    for m = 2 : M
        for n = 1 : m - 1
            mx = rem(m, Mx);
            if 0 == mx
                mx = Mx;
            else
            end
            my = (m - mx) / Mx + 1;
            nx = rem(n, Mx);
            if 0 == nx
                nx = Mx;
            else
            end
            ny = (n - nx) / Mx + 1;
            dx = mx - nx;
            dy = my - ny;
            Bx = phi_d(k, 1) ^ 2 * (cos(phi(k, 1)))^2 * (dx * cos(theta(k, 1)) + dy * sin(theta(k, 1))) ^ 2;
            By = theta_d(k, 1) ^ 2 * (sin(phi(k, 1)))^2 * (-dx * sin(theta(k, 1)) + dy * cos(theta(k, 1))) ^ 2;
            B(m, n) = exp(-0.5 * u ^ 2 * (Bx + By));
            B(n, m) = B(m, n);
        end
    end
    B_store((k - 1) * M + 1 : k * M, :) = B;
    R = R + sigma_s * D * B * D';
end
parti_theta = zeros(M * K, M);
parti_theta_d = zeros(M * K, M);
parti_phi = zeros(M * K, M);
parti_phi_d = zeros(M * K, M);
parti_s = zeros(M * K, M);
D_theta = zeros(M, M);
D_phi = zeros(M, M);
B_theta = zeros(M, M);
B_phi = zeros(M, M);
B_theta_d = zeros(M, M);
B_phi_d = zeros(M, M);
for k = 1 : K
    for mx = 1 : Mx
        for my = 1 : My
            m = mx + (my - 1) * Mx;
            D_theta(m, m) = 1i * u * sin(phi(k, 1)) * (-(mx - 1) * sin(theta(k, 1)) + (my - 1) * cos(theta(k, 1)));
            D_phi(m, m) = 1i * u * cos(phi(k, 1)) * ((mx - 1) * cos(theta(k, 1)) + (my - 1) * sin(theta(k, 1)));
        end
    end
    for m = 2 : M
        for n = 1 : m - 1
            mx = rem(m, Mx);
            if 0 == mx
                mx = Mx;
            else
            end
            my = (m - mx) / Mx + 1;
            nx = rem(n, Mx);
            if 0 == nx
                nx = Mx;
            else
            end
            ny = (n - nx) / Mx + 1;
            dx = mx - nx;
            dy = my - ny;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Bx = -(phi_d(k, 1)) ^ 2 * (cos(phi(k, 1)))^2 + (theta_d(k, 1)) ^ 2 * (sin(phi(k, 1)))^2;
            By = (dx^2 - dy^2) * sin(2 * theta(k, 1)) - 2 * dx * dy * cos(2 * theta(k, 1));
            B_theta(m, n) = -0.5 * u ^ 2 * (Bx * By);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Bx = -(phi_d(k, 1)) ^ 2 * (dx * cos(theta(k, 1)) + dy * sin(theta(k, 1)))^2;
            By = (theta_d(k, 1)) ^ 2 * (-dx * sin(theta(k, 1)) + dy * cos(theta(k, 1)))^2;
            B_phi(m, n) = -0.5 * u ^ 2 * sin(2 * phi(k, 1)) * (Bx + By);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            B_theta_d(m, n) = -u ^ 2 * theta_d(k, 1) * (sin(phi(k, 1)))^2 * (-dx * sin(theta(k, 1)) + dy * cos(theta(k, 1)))^2;
            B_phi_d(m, n) = -u ^ 2 * phi_d(k, 1) * (cos(phi(k, 1)))^2 * (dx * cos(theta(k, 1)) + dy * sin(theta(k, 1)))^2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            B_theta(n, m) = B_theta(m, n);
            B_phi(n, m) = B_phi(m, n);
            B_theta_d(n, m) = B_theta_d(m, n);
            B_phi_d(n, m) = B_phi_d(m, n);
        end
    end
    parti_theta((k - 1) * M + 1: k * M, :) = sigma_s * (D_theta * D_store((k - 1) * M + 1 : k * M, :) * B_store((k - 1) * M + 1 : k * M, :) * (D_store((k - 1) * M + 1 : k * M, :))' ...
        -D_store((k - 1) * M + 1 : k * M, :) * B_store((k - 1) * M + 1 : k * M, :) * (D_store((k - 1) * M + 1 : k * M, :))' * D_theta ...
        +D_store((k - 1) * M + 1 : k * M, :) * (B_store((k - 1) * M + 1 : k * M, :) .* B_theta) * (D_store((k - 1) * M + 1 : k * M, :))');
    parti_phi((k - 1) * M + 1: k * M, :) = sigma_s * (D_phi * D_store((k - 1) * M + 1 : k * M, :) * B_store((k - 1) * M + 1 : k * M, :) * (D_store((k - 1) * M + 1 : k * M, :))' ...
        -D_store((k - 1) * M + 1 : k * M, :) * B_store((k - 1) * M + 1 : k * M, :) * (D_store((k - 1) * M + 1 : k * M, :))' * D_phi ...
        +D_store((k - 1) * M + 1 : k * M, :) * (B_store((k - 1) * M + 1 : k * M, :) .* B_phi) * (D_store((k - 1) * M + 1 : k * M, :))');
    parti_theta_d((k - 1) * M + 1: k * M, :) = sigma_s * D_store((k - 1) * M + 1 : k * M, :) * (B_store((k - 1) * M + 1 : k * M, :) .* B_theta_d) * (D_store((k - 1) * M + 1 : k * M, :))';
    parti_phi_d((k - 1) * M + 1: k * M, :) = sigma_s * D_store((k - 1) * M + 1 : k * M, :) * (B_store((k - 1) * M + 1 : k * M, :) .* B_phi_d) * (D_store((k - 1) * M + 1 : k * M, :))';
    parti_s((k - 1) * M + 1: k * M, :) = D_store((k - 1) * M + 1 : k * M, :) * B_store((k - 1) * M + 1 : k * M, :) * (D_store((k - 1) * M + 1 : k * M, :))';
end

J_u_u = zeros(4 * K, 4 * K);
J_u_v = zeros(4 * K, K + 1);
J_v_v = zeros(K + 1, K + 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1 : K
    for m = 1 : n - 1
        J_u_u(n, m) =  Nt * trace(R \ parti_theta((n - 1) * M + 1: n * M, :) * (R \ parti_theta((m - 1) * M + 1: m * M, :)));
        J_u_u(m, n) = J_u_u(n, m);
    end
    m = n;
    J_u_u(n, m) =  Nt * trace(R \ parti_theta((n - 1) * M + 1: n * M, :) * (R \ parti_theta((m - 1) * M + 1: m * M, :)));
end
for n = K + 1 : 2 * K
    for m = 1 : K
        J_u_u(n, m) =  Nt * trace(R \ parti_phi((n - K - 1) * M + 1: (n - K) * M, :) * (R \ parti_theta((m - 1) * M + 1: m * M, :)));
        J_u_u(m, n) = J_u_u(n, m);
    end
    for m = K + 1 : n - 1
        J_u_u(n, m) =  Nt * trace(R \ parti_phi((n - K - 1) * M + 1: (n - K) * M, :) * (R \ parti_phi((m - K - 1) * M + 1: (m - K) * M, :)));
        J_u_u(m, n) = J_u_u(n, m);
    end
    m = n;
    J_u_u(n, m) =  Nt * trace(R \ parti_phi((n - K - 1) * M + 1: (n - K) * M, :) * (R \ parti_phi((m - K - 1) * M + 1: (m - K) * M, :)));
end
for n = 2 * K + 1 : 3 * K
    for m = 1 : K
        J_u_u(n, m) =  Nt * trace(R \ parti_theta_d((n - 2 * K - 1) * M + 1: (n - 2 * K) * M, :) * (R \ parti_theta((m - 1) * M + 1: m * M, :)));
        J_u_u(m, n) = J_u_u(n, m);
    end
    for m = K + 1 : 2 * K
        J_u_u(n, m) =  Nt * trace(R \ parti_theta_d((n - 2 * K - 1) * M + 1: (n - 2 * K) * M, :) * (R \ parti_phi((m - K - 1) * M + 1: (m - K) * M, :)));
        J_u_u(m, n) = J_u_u(n, m);
    end
    for m = 2 * K + 1 : n - 1
        J_u_u(n, m) =  Nt * trace(R \ parti_theta_d((n - 2 * K - 1) * M + 1: (n - 2 * K) * M, :) * (R \ parti_theta_d((m - 2 * K - 1) * M + 1: (m - 2 * K) * M, :)));
        J_u_u(m, n) = J_u_u(n, m);
    end
    m = n;
    J_u_u(n, m) =  Nt * trace(R \ parti_theta_d((n - 2 * K - 1) * M + 1: (n - 2 * K) * M, :) * (R \ parti_theta_d((m - 2 * K - 1) * M + 1: (m - 2 * K) * M, :)));
end
for n = 3 * K + 1 : 4 * K
    for m = 1 : K
        J_u_u(n, m) =  Nt * trace(R \ parti_phi_d((n - 3 * K - 1) * M + 1: (n - 3 * K) * M, :) * (R \ parti_theta((m - 1) * M + 1: m * M, :)));
        J_u_u(m, n) = J_u_u(n, m);
    end
    for m = K + 1 : 2 * K
        J_u_u(n, m) =  Nt * trace(R \ parti_phi_d((n - 3 * K - 1) * M + 1: (n - 3 * K) * M, :) * (R \ parti_phi((m - K - 1) * M + 1: (m - K) * M, :)));
        J_u_u(m, n) = J_u_u(n, m);
    end
    for m = 2 * K + 1 : 3 * K
        J_u_u(n, m) =  Nt * trace(R \ parti_phi_d((n - 3 * K - 1) * M + 1: (n - 3 * K) * M, :) * (R \ parti_theta_d((m - 2 * K - 1) * M + 1: (m - 2 * K) * M, :)));
        J_u_u(m, n) = J_u_u(n, m);
    end
    for m = 3 * K + 1 : n - 1
        J_u_u(n, m) =  Nt * trace(R \ parti_phi_d((n - 3 * K - 1) * M + 1: (n - 3 * K) * M, :) * (R \ parti_phi_d((m - 3 * K - 1) * M + 1: (m - 3 * K) * M, :)));
        J_u_u(m, n) = J_u_u(n, m);
    end
    m = n;
    J_u_u(n, m) =  Nt * trace(R \ parti_phi_d((n - 3 * K - 1) * M + 1: (n - 3 * K) * M, :) * (R \ parti_phi_d((m - 3 * K - 1) * M + 1: (m - 3 * K) * M, :)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1 : K
    for m = 1 : K
        J_u_v(n, m) =  Nt * trace(R \ parti_theta((n - 1) * M + 1: n * M, :) * (R \ parti_s((m - 1) * M + 1: m * M, :)));
    end
    J_u_v(n, K + 1) =  Nt * trace(R \ parti_theta((n - 1) * M + 1: n * M, :) / (R));
end
for n = K + 1 : 2 * K
    for m = 1 : K
        J_u_v(n, m) =  Nt * trace(R \ parti_phi((n - K - 1) * M + 1: (n - K) * M, :) * (R \ parti_s((m - 1) * M + 1: m * M, :)));
    end
    J_u_v(n, K + 1) =  Nt * trace(R \ parti_phi((n - K - 1) * M + 1: (n - K) * M, :) / (R));
end
for n = 2* K + 1 : 3 * K
    for m = 1 : K
        J_u_v(n, m) =  Nt * trace(R \ parti_theta_d((n - 2 * K - 1) * M + 1: (n - 2 * K) * M, :) * (R \ parti_s((m - 1) * M + 1: m * M, :)));
    end
    J_u_v(n, K + 1) =  Nt * trace(R \ parti_theta_d((n - 2 * K - 1) * M + 1: (n - 2 * K) * M, :) / (R));
end
for n = 3 * K + 1 : 4 * K
    for m = 1 : K
        J_u_v(n, m) =  Nt * trace(R \ parti_phi_d((n - 3 * K - 1) * M + 1: (n - 3 * K) * M, :) * (R \ parti_s((m - 1) * M + 1: m * M, :)));
    end
    J_u_v(n, K + 1) =  Nt * trace(R \ parti_phi_d((n - 3 * K - 1) * M + 1: (n - 3 * K) * M, :) / (R));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1 : K
    for m = 1 : n - 1
        J_v_v(n, m) =  Nt * trace(R \ parti_s((n - 1) * M + 1: n * M, :) * (R \ parti_s((m - 1) * M + 1: m * M, :)));
        J_v_v(m, n) = J_v_v(n, m);
    end
    m = n;
    J_v_v(n, m) =  Nt * trace(R \ parti_s((n - 1) * M + 1: n * M, :) * (R \ parti_s((m - 1) * M + 1: m * M, :)));
end
n = K + 1;
for m = 1 : K
    J_v_v(n, m) =  Nt * trace(R \ (R \ parti_s((m - 1) * M + 1: m * M, :)));
    J_v_v(m, n) = J_v_v(n, m);
end
m = K + 1;
J_v_v(n, m) =  Nt * trace((inv(R))^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = [J_u_u, J_u_v; J_u_v.', J_v_v];
CRB_eta = diag(inv(F));
% CRB_eta = diag(inv([J_u_u- J_u_v / J_v_v* J_u_v.']));
bound1 = (sqrt(real(CRB_eta(1 : 4 * K))) * 180 / pi).';
bound = [mean(bound1(1 : K)), mean(bound1(K+1 : 2*K)), mean(bound1(2*K+1 : 3*K)), mean(bound1(3*K+1 : 4*K))];





