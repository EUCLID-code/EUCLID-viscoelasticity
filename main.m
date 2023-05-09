%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Automated identification of linear viscoelastic constitutive laws    %%
%%                           with EUCLID                                 %%
%       Mechanics of Materials, Volume 181, June 2023, 104643             %
%            https://doi.org/10.1016/j.mechmat.2023.104643                %
%             Preprint: https://arxiv.org/abs/2212.10969                  %
%     Enzo Marino, Moritz Flaschel, Siddhant Kumar, Laura De Lorenzis     %
%                                                                         %
%                                                                         %
%                              Code by:                                   %
%                                                                         %
%                             Enzo Marino                                 %
%                    University of Florence (Italy)                       %
%                      email: enzo.marino@unifi.it                        %
%                                                                         %
%               -------------------------------------------               %
%                                                                         %
%              ________ ___  ___ _______ ___     ___ ______               %
%             /  _____//  / /  //  ____//  /    /  //  _   \              %
%            /  _____//  /_/  //  /___ /  /___ /  //  /_|  |              %
%           /_______//_______//______//______//__//_______/               %
% Efficient Unsupervised Constitutive Law Identification and Discovery    %
% https://euclid-code.github.io/                                          %
%                                                                         %
%                                                                         %
% Note: this code is not meant for professional use and it contains only  %
% very few comments. Interested users can contact me for additional       %
% information and support.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear all; close all;
addpath(genpath('routines'))
set(0,'DefaultFigureWindowStyle','docked')
filename = mfilename;

%% Set specimen geometry and mesh
Lx = 1; % Length in x direction [dm] (for symmetry I will model only Lx/2)
nx = 6; % Number of nodes in x, referred to only one half of the structure since symmetry is exploited
Ly = 5; % Length in y direction [dm]
ny = 51;% Number of nodex in y

[mesh_data]  = get_mesh_data(Lx,nx,Ly,ny);
coord = mesh_data.coord;
bottomnodes = mesh_data.bottomnodes;
intnodes = mesh_data.intnodes;
dofTy = mesh_data.dofTy;
dofLx = mesh_data.dofLx;
dofB = mesh_data.dofB;
conne = mesh_data.conne;
nx = mesh_data.nx;
ny = mesh_data.ny;
NN = size(coord,1);
NE = size(conne,1);

%% Set material
matselect = 3;% 00; %1 2 3 4 5

[target_mat] = get_true_material(matselect);
Ginf = target_mat.Ginf;
G    = target_mat.G;
tauG = target_mat.tauG;
Kinf = target_mat.Kinf;
K    = target_mat.K;
tauK = target_mat.tauK;
Gvec = [Ginf;G];
Kvec = [Kinf;K];
theta_true = [target_mat.Ginf;target_mat.G;target_mat.tauG;target_mat.Kinf;target_mat.K;target_mat.tauK];


%% Set TD and FD simulation settings
[ntimestep,Tsim,uc_th,time,dt,omevec,freqvec,indvec,fmin,fmax,df,phase_uybc,Amp_uybc] = TD_sim_setting_v04(target_mat);
std_ucth = std(uc_th);
%% True data generation: switch between FD or TD solver
solvertype_fwd_pb = 'FD'; % 'TD' time-domain solver not supported in this verion, please contact the author if you need it
add_noise = 0; % case with noise (add_noise = 1) not supported in this version, please contact the author if you need it
switch solvertype_fwd_pb
    case 'FD'
        % FD solver
        num_exc_ome = length(omevec);
        amp_omevec = Amp_uybc(indvec);
        pha_omevec = phase_uybc(indvec);
        [hatux, hatuy, hatrx, hatry, hatuxc, hatuyc, hatrxc, hatryc] =...
            main_FEM_2D_LinViscoElast_FD_v02(theta_true,omevec,amp_omevec,pha_omevec,...
            length(target_mat.tauG),length(target_mat.tauK),nx,ny,NN,coord,conne,bottomnodes,dofLx,dofTy);
        if add_noise == 1

            % Add noise here

        end

    case 'TD'

        % TD solver here

end

%% Discovery procedure
% Check zero cost
check_zero_cost = 0;
if check_zero_cost == 1
    tau_G = target_mat.tauG;
    tau_K = target_mat.tauK;
    [A_total,b_total,Acmplx,bcmplx] = get_system_enzo_v02(omevec,NN,coord,conne,bottomnodes,hatuxc,hatuyc,intnodes,tau_G,tau_K,hatrxc,hatryc);


    cost = A_total*[target_mat.Ginf;target_mat.G;target_mat.Kinf;target_mat.K] - b_total;
    mean(cost)

    linsolve(Acmplx,bcmplx)
    coeff = [target_mat.Ginf;target_mat.G;target_mat.Kinf;target_mat.K];
    ls = Acmplx*coeff;
    magls = abs(ls)
    magb  = abs(bcmplx)
    rls = real(ls)
    rb  = real(bcmplx)
    plot(real(ls));
    hold on
    plot(imag(ls))
    cost2 = Acmplx*[target_mat.Ginf;target_mat.G;target_mat.Kinf;target_mat.K] - bcmplx;
    mean(cost2)
end

% Optimization problem
% Assumed initinal number of Maxwell elements for G and K
NMeG = 300;
NMeK = 300;
incltruetimes = 0; % 1: includes true relaxation times; 0 does NOT include true relaxation times

tau_G = logspace(floor(log10(1/fmax)),ceil(log10(1/fmin)),NMeG)';
tau_K = logspace(floor(log10(1/fmax)),ceil(log10(1/fmin)),NMeK)';


if incltruetimes == 1
    if NMeG == length(target_mat.G) && NMeK == length(target_mat.K)
        tau_G = target_mat.tauG;
        tau_K = target_mat.tauK;
    else
        tau_G = union(tau_G,target_mat.tauG,'sorted');
        tau_K = union(tau_K,target_mat.tauK,'sorted');
    end
    NMeG = length(tau_G);
    NMeK = length(tau_K);% Just to test the capabilities of the identif method I include the true values in the relax times vecotr
end

tau_G
tau_K

% Get the linear system of equations (Eqs. (21) and (22) of the paper) 
[A_total, b_total, Acmplx, bcmplx] = get_system_enzo_v02(omevec,NN,coord,conne,bottomnodes,hatuxc,hatuyc,intnodes,tau_G,tau_K,hatrxc,hatryc);

Nlam = 8*1e3;
lammin = -12;
lammax = -1;
Lambda = logspace(lammin,lammax,Nlam);

% Solve linear system via Lasso regularization (Eq. (24) of the paper) 
[result,info] = lasso(A_total,b_total,'Lambda',Lambda,'MaxIter',1e8,'RelTol',1e-8,'Standardize',true);

hFig= figure;
for lam_idx = 1:size(result,2)
    num_nonzero_feats(lam_idx) =  nnz(result(:,lam_idx));
end
semilogx(squeeze(info.Lambda),num_nonzero_feats,'LineWidth',1); % 1/10 is needed for units consistency
hold on
xlabel('$\lambda$','Interpreter','LaTeX','FontSize',18)
ylabel('Number of nonzero features','Interpreter','latex','FontSize',18)
ylim([0 605])
xlim(10.^[lammin,lammax])
grid on
grid minor
box on
saveas(hFig,'num_nonzero_feats.fig')
% exportgraphics(hFig,'num_nonzero_feats.pdf','Resolution',600)
% exportgraphics(hFig,'num_nonzero_feats.png','Resolution',600)
close


MSE = squeeze(info.MSE);
lambdavec = squeeze(info.Lambda);
hFig = figure;
semilogx(lambdavec,MSE,'k-','LineWidth',1)
hold on
MSE_thr = 1e-5;
semilogx(squeeze(info.Lambda),MSE_thr*ones(1,Nlam),'k--','LineWidth',0.5)
xlabel('$\lambda$','Interpreter','LaTeX','FontSize',18)
ylabel('Mean Squared Error','Interpreter','LaTeX','FontSize',18)
ylim([0 MSE_thr*10])
xlim(10.^[lammin,lammax])
grid on
grid minor
legend('MSE','MSE threshold $e_\lambda$','Location','northwest', 'interpreter','LaTeX','FontSize',16)
box on
savefig(hFig,'MSE+threshold_vs_lam.fig')
% exportgraphics(hFig,'MSE+threshold_vs_lam.pdf','Resolution',600)
% exportgraphics(hFig,'MSE+threshold_vs_lam.png','Resolution',600)
close




hFig= figure;
for lam_idx = 1:size(result,2)
    num_nonzero_feats(lam_idx) =  nnz(result(:,lam_idx));
end
yyaxis left
semilogx(squeeze(info.Lambda),num_nonzero_feats,'LineWidth',1); %il /10 serve per porterle in N/mm^2
semilogx(squeeze(info.Lambda),num_nonzero_feats,'LineWidth',1); %il /10 serve per porterle in N/mm^2
ylim([0 605])
ylabel('Number of nonzero features','Interpreter','latex','FontSize',18)
yyaxis right
semilogx(squeeze(info.Lambda),info.MSE,'k-','LineWidth',1)
hold on
semilogx(squeeze(info.Lambda),MSE_thr*ones(1,Nlam),'k--','LineWidth',0.5)
ylabel('Mean Squared Error','Interpreter','LaTeX','FontSize',18)
ylim([0 MSE_thr*10])
xlabel('$\lambda$','Interpreter','latex','FontSize',18)
xlim(10.^[lammin,lammax])
legend('Number of nonzero features','MSE','MSE threshold $e_\lambda$','Location','best', 'interpreter','LaTeX','FontSize',16)
grid on
grid minor
box on
ax = gca;
ax.YAxis(2).Color = 'k';
saveas(hFig,'num_nonzero_feats_MSE.fig')
% exportgraphics(hFig,'num_nonzero_feats_MSE.pdf','Resolution',600)
% exportgraphics(hFig,'num_nonzero_feats_MSE.png','Resolution',600)
close

%% Activated features
err = squeeze(info.MSE); % note: already called MSE before. Fix it!
theta_idxopt = find(err > MSE_thr,1,'first');
theta_opt = result(:,theta_idxopt)
theta_opt2 = (theta_opt > 0.05*mean(theta_opt)).*theta_opt % Remove unsignificant features
G_activated = theta_opt2(2:NMeG+1);
tau_G_activated = tau_G.*(theta_opt2(2:NMeG+1)~=0);
K_activated = theta_opt2(NMeG+3:NMeG+NMeK+2);
tau_K_activated = tau_K.*(theta_opt2(NMeG+3:NMeG+NMeK+2)~=0);
Ginf_dsc = theta_opt2(1)
Kinf_dsc = theta_opt2(NMeG+2)


% Activated features plots
hFig = figure;
stem([Ginf_dsc;G_activated]/10,'-*','LineWidth',0.5,'MarkerSize',9)
ylabel('Shear moduli $\mathrm{[N/mm^2]}$','Interpreter','latex','FontSize',18)
xlabel('Features position','Interpreter','latex','FontSize',18)
xlim([0 300]);
grid on
grid minor
box on
ax = gca;
ax.YAxis.Exponent = 3;
savefig(hFig,'Activated G moduli.fig')
% exportgraphics(hFig,'Activated G moduli.pdf','Resolution',600)
% exportgraphics(hFig,'Activated G moduli.png','Resolution',600)
close

hFig = figure;
stem([Kinf_dsc;K_activated]/10,'-*','LineWidth',0.5,'MarkerSize',9)
ylabel('Bulk moduli $\mathrm{[N/mm^2]}$','Interpreter','latex','FontSize',18);
xlabel('Features position','Interpreter','latex','FontSize',18);
xlim([0 300]);
grid on
grid minor
box on
ax = gca;
ax.YAxis.Exponent = 3;
savefig(hFig,'Activated K moduli.fig')
% exportgraphics(hFig,'Activated K moduli.pdf','Resolution',600)
% exportgraphics(hFig,'Activated K moduli.png','Resolution',600)
close


hFig = figure;
stem(nonzeros(tau_G_activated),nonzeros(G_activated)/10,'*-','LineWidth',2,'MarkerSize',9)
hold on
stem(tauG,G/10,'-','LineWidth',2,'MarkerSize',9)
set(gca,'xscal','log')
legend('Discovered', 'True','Interpreter','latex','FontSize',16)
ylabel('Shear moduli $\mathrm{[N/mm^2]}$','Interpreter','latex','FontSize',18)
xlabel('Relaxation times  $\mathrm{[s]}$','Interpreter','latex','FontSize',18)
ax = gca;
ax.YAxis.Exponent = 3;
grid on
grid minor
box on
savefig(hFig,'Activated G vs relaxation times.fig')
% exportgraphics(hFig,'Activated G vs relaxation times.pdf','Resolution',600)
% exportgraphics(hFig,'Activated G vs relaxation times.png','Resolution',600)
close

hFig = figure;
stem(nonzeros(tau_K_activated),nonzeros(K_activated)/10,'*-','LineWidth',2,'MarkerSize',9)
hold on
stem(tauK,K/10,'-','LineWidth',2,'MarkerSize',9)
legend('Discovered', 'True','Interpreter','latex','FontSize',16)
set(gca,'xscal','log')
ylabel('Bulk moduli $\mathrm{[N/mm^2]}$','Interpreter','latex','FontSize',18)
xlabel('Relaxation times  $\mathrm{[s]}$','Interpreter','latex','FontSize',18)
ax = gca;
ax.YAxis.Exponent = 3;
grid on
grid minor
box on
savefig(hFig,'Activated K vs relaxation times.fig')
% exportgraphics(hFig,'Activated K vs relaxation times.pdf','Resolution',600)
% exportgraphics(hFig,'Activated K vs relaxation times.png','Resolution',600)
close


% Clustering
[S,cost_vec] = clustering(tau_G_activated,G_activated,tau_K_activated,K_activated,Ginf_dsc,Kinf_dsc,...
    omevec,NN,coord,conne,bottomnodes,hatuxc,hatuyc,intnodes,hatrxc,hatryc);

TrueParam = [[Ginf;G;Kinf;K],[NaN;tauG;NaN;tauK]];
S(5).param
return