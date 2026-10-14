function [psi, x, y, iter]=PB_solver(Lx, Ly, Nx, Ny, lambda, zeta, tol, alpha)
% Solves  d2psi/dx2 + d2psi/dy2 = lambda^2*sinh(psi),  psi = zeta on all four walls
% Source term is linearised about the previous iterate psi_old:
%   lambda^2*sinh(psi) ~ lambda^2*cosh(psi_old)*psi + lambda^2*(sinh(psi_old)-psi_old*cosh(psi_old))
% so the equation solved each iteration is   Lap(psi) - f1*psi = f2
    x=linspace(0, Lx, Nx);
    y=linspace(0, Ly, Ny);
    dx=x(2)-x(1);
    dy=y(2)-y(1);

    psi_old=zeta*ones(Nx, Ny);      % initial guess
    iter=0;
    err=1;
    while err>tol && iter<500
        iter=iter+1;
        A=spalloc(Nx*Ny, Nx*Ny, 5*Nx*Ny);
        rhs_val=zeros(Nx*Ny, 1);

        for i=1:Nx % bottom boundary
            k=i;
            A(k,k)=1;
            rhs_val(k,1)=zeta;
        end
        for i=1:Nx % top boundary
            k=(Ny-1)*Nx+i;
            A(k,k)=1;
            rhs_val(k,1)=zeta;
        end
        for j=2:Ny-1 % left boundary
            k=(j-1)*Nx+1;
            A(k,k)=1;
            rhs_val(k,1)=zeta;
        end
        for j=2:Ny-1 % right boundary
            k=j*Nx;
            A(k,k)=1;
            rhs_val(k,1)=zeta;
        end
        for i=2:Nx-1
            for j=2:Ny-1
                k=(j-1)*Nx+i;
                p=psi_old(i,j);
                f1=lambda^2*cosh(p);
                f2=lambda^2*(sinh(p)-p*cosh(p));
                A(k,k)=-2*dy^2/dx^2-2-f1*dy^2;
                A(k,k-1)=dy^2/dx^2;
                A(k,k+1)=dy^2/dx^2;
                A(k,k-Nx)=1;
                A(k,k+Nx)=1;
                rhs_val(k,1)=f2*dy^2;
            end
        end

        psi_vec=A\rhs_val;
        psi_new=reshape(psi_vec, Nx, Ny);
        err=max(max(abs(psi_new-psi_old)));
        psi_old=psi_old+alpha*(psi_new-psi_old);   % under-relaxation (alpha=1: none)
    end
    psi=psi_old;
end
