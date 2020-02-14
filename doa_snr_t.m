%This file is for DOA estimation in LS-MIMO angular spread scenarios

clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt = 500;%number of snapshots
preci = 0.2;%resolution of 0.1 degree, preci>=180/M
Sk = zeros(1, Nt);%transmitted signal for each terminal
rmse1_bound = zeros(2, 1);
rmse2_bound = zeros(2, 1);
Ite_num = 2e2;%number of iterations
Ite_num_store = [50;1e2;2e2;2e2;];
theta_range = 1;%the range of searched angle, degree
theta_d_max = 2;
phi_range = 1;
phi_d_max = 2;
u = 2 * pi * 0.5;%d/lambda = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 2;%user number
Nk = 50 * ones(K, 1);%number of paths
snr_n = 4;
snr_s = [-5;0;5;10];%dB
M = 100;
Mx = sqrt(M);
My = Mx;
bound_store = zeros(snr_n, 4);
rmse_store_theta = zeros(snr_n, 3);
rmse_store_phi = zeros(snr_n, 3);
rmse_store_theta_d = zeros(snr_n, 3);
rmse_store_phi_d = zeros(snr_n, 3);
theta = [10; 50];%nominal angles, in ascending order, 0~180 degree
theta_d = [1; 1];%angular deviation
phi = [30; 40];%nominal angles, in ascending order, 0~90 degree
phi_d = [1; 1];%angular deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn = 1 : snr_n
    Ite_num = Ite_num_store(nn, 1);
    snr = snr_s(nn, 1);
    %     Y = zeros(M, Nt);%received signal
    %     A = zeros(M, 3 * K);%steering vector
    a = zeros(M, 1);%steering vector
    N = zeros(M, Nt);%received noise at the BS
    rmse = zeros(4, 3);%4 parameters, 3 approaches
    for ii = 1 : Ite_num
        %generate the received signal
        Y1 = zeros(M, Nt);%received signal
        %         Y2 = zeros(M, Nt);%received signal
        N = (randn(M, Nt) + 1i * randn(M, Nt)) / sqrt(2);
        Y1 = Y1 + N;
        %Y2 = Y2 + N;
        S = zeros(3 * K, Nt);
        for k = 1 : K
            amp_k = sqrt(10 ^ (snr * 0.1) / (Nk(k, 1)));
            Sk = sign(randn(1, Nt));%BPSK modulation
            %%%%%%%%%%
            %             for mx = 1 : Mx
            %                 for my = 1 : My
            %                     m = mx + (my - 1) * Mx;
            %                     A(m, k) = exp(1i * u * sin(phi(k, 1) / 180 * pi) * ((mx - 1) * cos(theta(k, 1) / 180 * pi) + (my - 1) * sin(theta(k, 1) / 180 * pi)));
            %                     A(m, k + K) = A(m, k) * 1i * u * sin(phi(k, 1) / 180 * pi) * (-(mx - 1) * sin(theta(k, 1) / 180 * pi) + (my - 1) * cos(theta(k, 1) / 180 * pi));
            %                     A(m, k + 2 * K) = A(m, k) * 1i * u * cos(phi(k, 1) / 180 * pi) * ((mx - 1) * cos(theta(k, 1) / 180 * pi) + (my - 1) * sin(theta(k, 1) / 180 * pi));
            %                 end
            %             end
            %%%%%%%%%%%%%%%
            for jj = 1 : Nt
                if 0 == Sk(1, jj)
                    Sk(1, jj) = 1;
                else
                end
                Sk(1, jj) = Sk(1, jj) * amp_k;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                alpha_k = (randn(Nk(k, 1), 1) + 1i * randn(Nk(k, 1), 1)) / sqrt(2);%small scale fading
                theta_k = randn(Nk(k, 1), 1) * theta_d(k, 1) / 180 * pi;%Gaussian distribution
                phi_k = randn(Nk(k, 1), 1) * phi_d(k, 1) / 180 * pi;%Gaussian distribution
                %                 for ll = 1 : Nk(k, 1)
                %                     S(k, jj) = S(k, jj) + alpha_k(ll, 1);
                %                     S(k + K, jj) = S(k + K, jj) + alpha_k(ll, 1) * theta_k(ll, 1);
                %                     S(k + 2 * K, jj) = S(k + 2 * K, jj) + alpha_k(ll, 1) * phi_k(ll, 1);
                %                 end
                %                 S(k, jj) = S(k, jj) * Sk(1, jj);
                %                 S(k + K, jj) = S(k + K, jj) * Sk(1, jj);
                %                 S(k + 2 * K, jj) = S(k + 2 * K, jj) * Sk(1, jj);
                %%%%%%%%%%%%%%%%
                theta_k = theta_k * 180 / pi + theta(k, 1) * ones(Nk(k, 1), 1);
                phi_k = phi_k * 180 / pi + phi(k, 1) * ones(Nk(k, 1), 1);
                for ll = 1 : Nk(k, 1)
                    for mx = 1 : Mx
                        for my = 1 : My
                            m = mx + (my - 1) * Mx;
                            a(m, 1) = exp(1i * u * sin(phi_k(ll, 1) / 180 * pi) * ((mx - 1) * cos(theta_k(ll, 1) / 180 * pi) + (my - 1) * sin(theta_k(ll, 1) / 180 * pi)));
                        end
                    end
                    Y1(:, jj) = Y1(:, jj) + alpha_k(ll, 1) * a * Sk(1, jj);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        %         Y2 = Y2 +  A * randn(3*K, Nt);
        %         Y2 = Y2 +  A * S;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %estimation
        RY = zeros(M, M);
        for jj = 1 : Nt
            RY = RY + Y1(:, jj) * Y1(:, jj)';
        end
        RY = RY / Nt;
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %DISPARE [37]&Subspace
        xi = eig_sub(M, theta, phi, RY, preci, theta_d_max, theta_range, phi_d_max, phi_range, u);
        for k = 1 : K
            rmse(:, 1) = rmse(:, 1) + (xi(:, k) - [theta(k, 1); theta_d(k, 1); phi(k, 1); phi_d(k, 1)]).^2;
            rmse(:, 2) = rmse(:, 2) + (xi(:, K + k) - [theta(k, 1); theta_d(k, 1); phi(k, 1); phi_d(k, 1)]).^2;
        end
        %Proposed approach
        eta = esprit_t2(RY, K, u, Mx, My);
        rmse(:, 3) = rmse(:, 3) + (abs(eta(:, 1) - [theta(2, 1); theta_d(2, 1); phi(2, 1); phi_d(2, 1)])).^2;
        rmse(:, 3) = rmse(:, 3) + (abs(eta(:, 2) - [theta(1, 1); theta_d(1, 1); phi(1, 1); phi_d(1, 1)])).^2;
        sprintf('%d,%d', nn, ii)
    end
    %%%%%%%%%%%%%%%%%%%
    bound_store(nn, :) = crb(Nt, K, M, theta, theta_d, phi, phi_d, snr, u);%theta,phi,theta_d,phi_d
    rmse_store_theta(nn, :) = sqrt(rmse(1, :) / Ite_num / K);
    rmse_store_theta_d(nn, :) = sqrt(rmse(2, :) / Ite_num / K);
    rmse_store_phi(nn, :) = sqrt(rmse(3, :) / Ite_num / K);
    rmse_store_phi_d(nn, :) = sqrt(rmse(4, :) / Ite_num / K);
    sprintf('%d', nn)
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
set(h,'PaperType','A4');
axes('FontSize',16);
subplot(2,2,1)
semilogy(snr_s,rmse_store_theta(1:snr_n,1),'b-.o','LineWidth',1,'MarkerSize',6)
hold on
semilogy(snr_s,rmse_store_theta(1:snr_n,2),'m-.x','LineWidth',1,'MarkerSize',6)
semilogy(snr_s,rmse_store_theta(1:snr_n,3),'k-.*','LineWidth',1,'MarkerSize',6)
semilogy(snr_s,bound_store(1:snr_n,1),'rs-','LineWidth',1,'MarkerSize',6)
grid on
le = legend('DISPARE [37]','Subspace [40]','Proposed','CRB (72)', 'Location','Northeast');
set(le,'Fontsize',8,'Fontname','Times')
set(gca,'XTick',snr_s)
xlim([min(snr_s), max(snr_s)])
xlabel({'Average received SNR (dB)','(a)'},'Fontsize',10,'Fontname','Times')
ylabel('RMSE of \theta estimate  (degrees)','Fontsize',10,'Fontname','Times')
%%%%%%%%%%%
subplot(2,2,2)
semilogy(snr_s,rmse_store_phi(1:snr_n,1),'b-.o','LineWidth',1,'MarkerSize',6)
hold on
semilogy(snr_s,rmse_store_phi(1:snr_n,2),'m-.x','LineWidth',1,'MarkerSize',6)
semilogy(snr_s,rmse_store_phi(1:snr_n,3),'k-.*','LineWidth',1,'MarkerSize',6)
semilogy(snr_s,bound_store(1:snr_n,2),'rs-','LineWidth',1,'MarkerSize',6)
grid on
le = legend('DISPARE [37]','Subspace [40]','Proposed','CRB (72)', 'Location','Northeast');
set(le,'Fontsize',8,'Fontname','Times')
set(gca,'XTick',snr_s)
xlim([min(snr_s), max(snr_s)])
xlabel({'Average received SNR (dB)','(b)'},'Fontsize',10,'Fontname','Times')
ylabel('RMSE of \phi estimate  (degrees)','Fontsize',10,'Fontname','Times')
%%%%%%%%%%%
subplot(2,2,3)
semilogy(snr_s,rmse_store_theta_d(1:snr_n,1),'b-.o','LineWidth',1,'MarkerSize',6)
hold on
semilogy(snr_s,rmse_store_theta_d(1:snr_n,2),'m-.x','LineWidth',1,'MarkerSize',6)
semilogy(snr_s,rmse_store_theta_d(1:snr_n,3),'k-.*','LineWidth',1,'MarkerSize',6)
semilogy(snr_s,bound_store(1:snr_n,3),'rs-','LineWidth',1,'MarkerSize',6)
grid on
le = legend('DISPARE [37]','Subspace [40]','Proposed','CRB (72)', 'Location','Northeast');
set(le,'Fontsize',8,'Fontname','Times')
set(gca,'XTick',snr_s)
xlim([min(snr_s), max(snr_s)])
xlabel({'Average received SNR (dB)','(c)'},'Fontsize',10,'Fontname','Times')
ylabel('RMSE of \sigma_{\theta} estimate (degrees)','Fontsize',10,'Fontname','Times')
%%%%%%%%%%%
subplot(2,2,4)
semilogy(snr_s,rmse_store_phi_d(1:snr_n,1),'b-.o','LineWidth',1,'MarkerSize',6)
hold on
semilogy(snr_s,rmse_store_phi_d(1:snr_n,2),'m-.x','LineWidth',1,'MarkerSize',6)
semilogy(snr_s,rmse_store_phi_d(1:snr_n,3),'k-.*','LineWidth',1,'MarkerSize',6)
semilogy(snr_s,bound_store(1:snr_n,4),'rs-','LineWidth',1,'MarkerSize',6)
grid on
le = legend('DISPARE [37]','Subspace [40]','Proposed','CRB (72)', 'Location','Northeast');
set(le,'Fontsize',8,'Fontname','Times')
set(gca,'XTick',snr_s)
xlim([min(snr_s), max(snr_s)])
xlabel({'Average received SNR (dB)','(d)'},'Fontsize',10,'Fontname','Times')
ylabel('RMSE of \sigma_{\phi} estimate (degrees)','Fontsize',10,'Fontname','Times')
% print(h,'-dpdf','RMSE_DOA_SNR')

