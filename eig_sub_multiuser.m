function xi = eig_sub_multiuser(M, K, theta, phi, RY, preci, theta_d_max, theta_range, phi_d_max, phi_range, u)
[phi_sort, IX] = sort(phi(1:K),'ascend');
theta_sort = theta(IX);
Mx = sqrt(M);
Eta = zeros(4, 1);
%%%%%%%%%%%%%%%%%%

[E, D] = eig(RY);
lambda = diag(D);
temp = sum(lambda) * 0.95;
accu = 0;
for jj = M : -1 : 1
    accu = accu + lambda(jj, 1);
    if accu >= temp
        break;
    else
    end
end
E_noise = E(:, 1 : max(1, jj - 1));
xi = zeros(4, 2*K);
for k_user = 1 : K
    alpha_min1 = theta_sort(k_user, 1) - theta_range;
    alpha_max1 = theta_sort(k_user, 1) + theta_range;
    alpha_num1 = (alpha_max1 - alpha_min1) / preci + 1;
    alpha_num1 = ceil(alpha_num1);
    alpha_dnum = theta_d_max / preci;%doa spread range 1`10 degrees
    %%%%%%%%%%%%%%%%%%%%%%
    beta_min1 = phi_sort(k_user, 1) - phi_range;
    beta_max1 = phi_sort(k_user, 1) + phi_range;
    beta_num1 = (beta_max1 - beta_min1) / preci + 1;
    beta_num1 = ceil(beta_num1);
    beta_dnum = phi_d_max / preci;%doa spread range 1`10 degrees
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Lambda1 = zeros(alpha_num1 * alpha_dnum * beta_num1 * beta_dnum, 1);%store of the correlations
    Lambda_sub1 = zeros(alpha_num1 * alpha_dnum * beta_num1 * beta_dnum, 1);%store of the correlations
    for ii = 1 : alpha_num1
        thetai = (alpha_min1 + (ii - 1) * preci) / 180 * pi;
        for jj = 1 : alpha_dnum
            theta_d = jj * preci / 180 * pi;
            for kk = 1 : beta_num1
                phii = (beta_min1 + (kk - 1) * preci) / 180 * pi;
                for ll = 1 : beta_dnum
                    phi_d = ll * preci / 180 * pi;
                    Psi = eye(M);
                    for k = 2 : M
                        for l = 1 : k - 1
                            mx = rem(k, Mx);
                            if 0 == mx
                                mx = Mx;
                            else
                            end
                            my = (k - mx) / Mx + 1;
                            nx = rem(l, Mx);
                            if 0 == nx
                                nx = Mx;
                            else
                            end
                            ny = (l - nx) / Mx + 1;
                            dx = mx - nx;
                            dy = my - ny;
                            Psi(k, l) = exp(1i * u * sin(phii) * (dx * cos(thetai) + dy * sin(thetai))) ...
                                * exp(-0.5 * u ^ 2 * (phi_d ^ 2 * cos(phii) ^ 2 * (dx * cos(thetai) + dy * sin(thetai)) ^ 2 ...
                                + theta_d ^ 2 * sin(phii) ^ 2 * (-dx * sin(thetai) + dy * cos(thetai)) ^ 2));%gaussian
                            Psi(l, k) = conj(Psi(k, l));
                        end
                    end
                    Psi1 = E_noise' * Psi;
                    Lambda1((((ii - 1) * alpha_dnum + jj - 1) * beta_num1 + kk - 1) * beta_dnum + ll, 1) = sum(diag(Psi1' * Psi1));
                    Psi1 = RY \ Psi;
                    Lambda_sub1((((ii - 1) * alpha_dnum + jj - 1) * beta_num1 + kk - 1) * beta_dnum + ll, 1) = sum(diag(Psi1' * Psi1));
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    [~, IX] = sort(Lambda1);
    IX1 = IX(1, 1);
    ll = rem(IX1, beta_dnum);
    if 0 == ll
        ll = beta_dnum;
    else
    end
    IX1 = (IX1 - ll) / beta_dnum;
    kk = rem(IX1, beta_num1) + 1;
    IX1 = (IX1 - kk + 1) / beta_num1;
    jj = rem(IX1, alpha_dnum) + 1;
    IX1 = (IX1 - jj + 1) / alpha_dnum;
    ii = IX1 + 1;
    Eta(1, 1) =  alpha_min1 + (ii - 1) * preci;
    Eta(2, 1) = jj * preci;
    Eta(3, 1) =  beta_min1 + (kk - 1) * preci;
    Eta(4, 1) = ll * preci;
    xi(:, k_user) = Eta;
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, IX] = sort(Lambda_sub1);
    IX1 = IX(1, 1);
    ll = rem(IX1, beta_dnum);
    if 0 == ll
        ll = beta_dnum;
    else
    end
    IX1 = (IX1 - ll) / beta_dnum;
    kk = rem(IX1, beta_num1) + 1;
    IX1 = (IX1 - kk + 1) / beta_num1;
    jj = rem(IX1, alpha_dnum) + 1;
    IX1 = (IX1 - jj + 1) / alpha_dnum;
    ii = IX1 + 1;
    Eta(1, 1) =  alpha_min1 + (ii - 1) * preci;
    Eta(2, 1) = jj * preci;
    Eta(3, 1) =  beta_min1 + (kk - 1) * preci;
    Eta(4, 1) = ll * preci;
    xi(:, k_user+K) = Eta;
end
