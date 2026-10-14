clear; clc; close all;
Lx=1; Ly=1;
Nx=200; Ny=200;
lambdas=[5 10 20 30 40];
zetas=[0.5 1 3 5];
tol=1e-8;
alpha=1;               % use ~0.7 if the iteration oscillates

for a=1:length(lambdas)
    for b=1:length(zetas)
        lambda=lambdas(a);
        zeta=zetas(b);
        [psi, x, y, iter]=PB_FVM_solver(Lx, Ly, Nx, Ny, lambda, zeta, tol, alpha);
        cp=exp(-psi);      % c+/c0
        cm=exp(psi);       % c-/c0
        fprintf('lambda=%g, zeta=%g : %d iterations\n', lambda, zeta, iter);

        figure('Name', sprintf('lambda=%g zeta=%g', lambda, zeta));
        subplot(1,3,1); contourf(x, y, psi', 30); colorbar; axis square;
        title('\psi_*'); xlabel('x_*'); ylabel('y_*');
        subplot(1,3,2); contourf(x, y, cp', 30); colorbar; axis square;
        title('c_{+}/c_0'); xlabel('x_*'); ylabel('y_*');
        subplot(1,3,3); contourf(x, y, cm', 30); colorbar; axis square;
        title('c_{-}/c_0'); xlabel('x_*'); ylabel('y_*');
        sgtitle(sprintf('FVM: \\lambda=%g, \\zeta_*=%g', lambda, zeta));

        pos=[0.1 0.25 0.5];
        figure('Name', sprintf('lines lambda=%g zeta=%g', lambda, zeta));
        subplot(1,2,1); hold on;
        for m=1:length(pos)
            [~, jj]=min(abs(y-pos(m)*Ly));
            plot(x, psi(:,jj), 'DisplayName', sprintf('\\psi, y_*=%.2f', y(jj)));
            plot(x, cp(:,jj), '--', 'DisplayName', sprintf('c_+, y_*=%.2f', y(jj)));
            plot(x, cm(:,jj), ':', 'DisplayName', sprintf('c_-, y_*=%.2f', y(jj)));
        end
        xlabel('x_*'); title('profiles vs x_*'); legend('Location','best'); grid on;
        subplot(1,2,2); hold on;
        for m=1:length(pos)
            [~, ii]=min(abs(x-pos(m)*Lx));
            plot(y, psi(ii,:), 'DisplayName', sprintf('\\psi, x_*=%.2f', x(ii)));
            plot(y, cp(ii,:), '--', 'DisplayName', sprintf('c_+, x_*=%.2f', x(ii)));
            plot(y, cm(ii,:), ':', 'DisplayName', sprintf('c_-, x_*=%.2f', x(ii)));
        end
        xlabel('y_*'); title('profiles vs y_*'); legend('Location','best'); grid on;
        sgtitle(sprintf('FVM: \\lambda=%g, \\zeta_*=%g', lambda, zeta));
    end
end
