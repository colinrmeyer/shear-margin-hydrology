function [out,q] = Nfun2d(dx,dy,Nx,Ny,I,EffP,phi,alpha,eta,delta,g,kappa,basalN)
N = Nx*Ny; % Total number of entries

% Diffusion matrix
lambda.matrix = kappa*delta*(phi.^(alpha));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
Tran.xplus = (2*dy/dx)*([lambda.matrix(2:Nx,:);zeros(1,Ny)].^(-1)+lambda.matrix.^(-1)).^(-1);
Tran.xminus = (2*dy/dx)*([zeros(1,Ny);lambda.matrix(1:(Nx-1),:)].^(-1)+lambda.matrix.^(-1)).^(-1);
Tran.yplus = (2*dx/dy)*([lambda.matrix(:,2:Ny),zeros(Nx,1)].^(-1)+lambda.matrix.^(-1)).^(-1);
Tran.yminus = (2*dx/dy)*([zeros(Nx,1),lambda.matrix(:,1:(Ny-1))].^(-1)+lambda.matrix.^(-1)).^(-1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Vectorize % %
TvecXp = reshape(Tran.xplus',N,1);
TvecXm = reshape(Tran.xminus',N,1);
TvecYp = reshape(Tran.yplus',N,1);
TvecYm = reshape(Tran.yminus',N,1);
% % add in the boundary condition % %
Tran.yplusBC = (2*dx/dy)*([lambda.matrix(:,2:Ny),Inf*ones(Nx,1)].^(-1)+lambda.matrix.^(-1)).^(-1);
Tran.yminusBC = (2*dx/dy)*([Inf*ones(Nx,1),lambda.matrix(:,1:(Ny-1))].^(-1)+lambda.matrix.^(-1)).^(-1);
TvecYpBC = reshape(Tran.yplusBC',N,1);
TvecYmBC = reshape(Tran.yminusBC',N,1);
SigmaT = TvecXp+TvecXm+TvecYpBC+TvecYmBC;
% % Compile into sparse matrix % %
S = [TvecXm((Ny+1):N)', TvecYm(2:N)', -SigmaT', TvecYp(1:(N-1))', TvecXp(1:(N-Ny))'];
Ind = [(Ny+1):N, 2:N, 1:N, 1:(N-1), 1:(N-Ny)];
Jnd = [1:(N-Ny), 1:(N-1), 1:N, 2:N, (Ny+1):N];
D = sparse(Ind,Jnd,S,N,N);

% % Sheat heating matrix % %
S = (reshape(phi',N,1).*dx.*dy./(reshape(eta',N,1)))';
Ind = 1:N; Jnd = 1:N;
H = sparse(Ind,Jnd,S,N,N);

% % Drainage % %
alam.matrix = (phi.^(alpha));
phialpha_plus = 2*([alam.matrix(:,2:Ny),alam.matrix(:,Ny)].^(-1)+alam.matrix.^(-1)).^(-1);
phialpha_minus = 2*([alam.matrix(:,1),alam.matrix(:,1:(Ny-1))].^(-1)+alam.matrix.^(-1)).^(-1);
Drain = (g*kappa*dx)*(phialpha_plus-phialpha_minus);
Drainvec = reshape(Drain',N,1);

% % boundary conditions % %
bcsmat = zeros(Nx,Ny);
bcsmat(:,1) = basalN.*Tran.yminusBC(:,1);
bcsmat(:,Ny) = Tran.yplusBC(:,Ny);
bcsvec = reshape(bcsmat',N,1);
% % compute effective pressure % %
M = D-H;
EffP(I) = M(I,I)\(-Drainvec(I)-bcsvec(I)+M(I,~I)*EffP(~I));
EffPmat = reshape(EffP,Ny,Nx);
out = EffPmat; 
q = (g*kappa*phialpha_minus(:,1))+(Tran.yminusBC(:,1)/dx).*(EffPmat(1,:)'-basalN); % not integrated across dy
end

