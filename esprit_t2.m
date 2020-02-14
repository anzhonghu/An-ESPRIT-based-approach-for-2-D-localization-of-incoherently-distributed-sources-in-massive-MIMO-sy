function eta = esprit_t2(RY, K, u, Mx, My)



[U, D] = eig(RY);
D = diag(D);
[D, IX] = sort(D, 'descend');
Us = U(:, IX(1 : 3 * K, 1));
U1 = zeros((Mx - 1) * (My - 1), 3 * K);
U2 = zeros((Mx - 1) * (My - 1), 3 * K);
U3 = zeros((Mx - 1) * (My - 1), 3 * K);
for ii = 1 : Mx - 1
    for jj = 1 : My - 1
        U1((jj - 1) * (Mx - 1) + ii, :) = Us((jj - 1) * Mx + ii, :);
        U2((jj - 1) * (Mx - 1) + ii, :) = Us((jj - 1) * Mx + ii + 1, :);
        U3((jj - 1) * (Mx - 1) + ii, :) = Us(jj * Mx + ii, :);
    end
end

U112 = [U1, U2];
U212 = [U1, U3];
[E1, D1] = eig(U112' * U112);
[E2, D2] = eig(U212' * U212);
D1 = diag(D1);
D2 = diag(D2);
[~, IX] = sort(D1, 'descend');
Psi1 = -E1(1 : 3 * K, IX(3 * K + 1 : 6 * K, 1)) / E1(3 * K + 1 : 6 * K, IX(3 * K + 1 : 6 * K, 1));
[~, IX] = sort(D2, 'descend');
Psi2 = -E2(1 : 3 * K, IX(3 * K + 1 : 6 * K, 1)) / E2(3 * K + 1 : 6 * K, IX(3 * K + 1 : 6 * K, 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [U,~] = eig(Psi1);
% Phi1_t = U \ Psi1 * U;
% Phi2_t = U \ Psi2 * U;
% Phi1 = log(diag(Phi1_t));
% Phi2 = log(diag(Phi2_t));
% Phi1 = imag(Phi1);
% Phi2 = imag(Phi2);
%%%%%%%%%%%%%%%%%%%%%%
[~, Phi1_t] = eig(Psi1);
[~, Phi2_t] = eig(Psi2);
[~, Phi3_t] = eig(Psi1 * Psi2);
Phi2_tn = zeros(3 * K, 3* K);
for k = 1 : 3 * K
    beta_kk1 = zeros(1, 3 * K);
    for k1 = 1 : 3 * K
        beta_kk1(1, k1) = Phi1_t(k, k) * Phi2_t(k1, k1);
    end
    miu_k1k2 = zeros(3 * K * 3 * K, 1);
    for k1 = 1 : 3 * K
        for k2 = 1 : 3 * K
            miu_k1k2((k1 - 1) * 3 * K + k2) = abs(beta_kk1(1, k1) - Phi3_t(k2, k2))^2;
        end
    end
    [~, index] = min(miu_k1k2);
    if 0 == rem(index, 3 * K)
        index = index / 3 / K;
    else
        index = floor(index / 3 / K) + 1;
    end
    Phi2_tn(k, k) = Phi2_t(index, index);
end
Phi1 = log(diag(Phi1_t));
Phi2 = log(diag(Phi2_tn));
Phi1 = imag(Phi1);
Phi2 = imag(Phi2);
%%%%%%%%%%%%%%%%%



[Phi1,Ix]=sort(Phi1,'ascend');
Phi2 = Phi2(Ix);
theta_t = atan(Phi2 ./ Phi1);
phi_t = asin(sqrt(Phi1.^2 + Phi2.^2) /  u);
theta = zeros(K, 1);
phi = zeros(K, 1);
for k = 1 : K
    theta(k, 1) = mean(theta_t((k - 1) * 3 + 1 : k * 3, 1));
    phi(k, 1) = mean(phi_t((k - 1) * 3 + 1 : k * 3, 1));
    if theta(k, 1) < 0
        theta(k, 1) = theta(k, 1) + pi;
    else
    end
end


A = zeros((My) * (Mx), 3 * K);
for k = 1 : K
    for ii = 1 : My
        for jj = 1 : Mx
            A((ii - 1) * (Mx) + jj, k) = exp(1i * u * sin(phi(k, 1)) * ((jj-1) * cos(theta(k, 1)) + (ii-1) * sin(theta(k, 1))));
            A((ii - 1) * (Mx) + jj, k + K) = A((ii - 1) * (Mx) + jj, k) * (1i * u * sin(phi(k, 1)) * (-(jj - 1) * sin(theta(k, 1)) + (ii - 1) * cos(theta(k, 1))));
            A((ii - 1) * (Mx) + jj, k + 2 * K) = A((ii - 1) * (Mx) + jj, k) * (1i * u * cos(phi(k, 1)) * ((jj - 1) * cos(theta(k, 1)) + (ii - 1) * sin(theta(k, 1))));
        end
    end
end
sigma_n2 = mean(D(3 * K + 1 : Mx * My, 1));
Lambdas = (A' * A) \ A' * (RY - sigma_n2 * eye((My) * (Mx))) * A / (A' * A);
Lambdas = real(diag(Lambdas));
theta_d = zeros(K, 1);
phi_d = zeros(K, 1);
for k = 1 : K
    theta_d(k, 1) = sqrt(Lambdas(K + k, 1) ./ Lambdas(k, 1));
    phi_d(k, 1) = sqrt(Lambdas(2 * K + k, 1) ./ Lambdas(k, 1));
end

theta = theta * 180 / pi;
phi = phi * 180 / pi;
theta_d = theta_d * 180 / pi;
phi_d = phi_d * 180 / pi;
eta = real([theta'; theta_d'; phi'; phi_d']);


