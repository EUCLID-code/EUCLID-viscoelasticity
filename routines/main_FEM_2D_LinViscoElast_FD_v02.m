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
% function [deformationx_amp_ome, deformationy_amp_ome, RBottomx_amp_ome, RBottomy_amp_ome,...
% deformationx_ome,deformationy_ome, RBottomx_ome, RBottomy_ome] = main_FEM_2D_LinViscoElast_FD_v02(theta,omevec,amp_omevec,pha_omevec,NMeG,NMeK,nx,ny,NN,coord,conne,bottomnodes,dofLx,dofTy)

function [deformationx_amp_ome, deformationy_amp_ome, RBottomx_amp_ome, RBottomy_amp_ome,...
   deformationx_ome,deformationy_ome, RBottomx_ome, RBottomy_ome] = main_FEM_2D_LinViscoElast_FD_v02(theta,omevec,amp_omevec,pha_omevec,NMeG,NMeK,nx,ny,NN,coord,conne,bottomnodes,dofLx,dofTy)


%% Material
% Deviatoric
G     = theta(2:NMeG+1); 
tauG  = theta(NMeG+2:2*NMeG+1); 
Ginf = theta(1);
G0 = Ginf + sum(G);
g = G/G0; % Abaqus form
Gvec = [Ginf;G];

% Volumetric
K     = theta(2*NMeG+3:2*NMeG+NMeK+2);
tauK  = theta(2*NMeG+2+NMeK+1:2*NMeG+2*NMeK+2); 
Kinf = theta(2*NMeG+2); 
K0 = Kinf + sum(K);
k = K/K0; % Abaqus form
Kvec = [Kinf;K];




%% Projectors in 2D for plane strain problm
m = [1 1 1 0]';
PrD =  eye(4) - 1/3*(m*m');
PrV =1/3*(m*m');
Dmu = [2 0 0 0;
    0 2 0 0;
    0 0 2 0;
    0 0 0 1];




%% FD solver
for iome = 1:length(omevec)
    omega = omevec(iome);
   
    BGs = [1;omega.^2.*tauG.^2./(1+omega.^2.*tauG.^2)];
    BGl = [0;omega.*tauG./(1+omega^2.*tauG.^2)];
    
    BKs = [1;omega.^2*tauK.^2./(1+omega.^2.*tauK.^2)];
    BKl = [0;omega.*tauK./(1+omega.^2.*tauK.^2)];
    
   
    %% FD solver
    KK = zeros(2*NN);
    Fvol = zeros(2*NN,1);
    for ie = 1:size(conne,1)
        vrtx = conne(ie,1:3); % vertexes array
        X = [coord(vrtx(1),2),coord(vrtx(2),2),coord(vrtx(3),2)]';
        Y = [coord(vrtx(1),3),coord(vrtx(2),3),coord(vrtx(3),3)]';
        N_xi  = [-1 1 0]';
        N_eta = [-1 0 1]';
        Jmat(1,1) = N_xi'*X; Jmat(1,2) = N_xi'*Y;
        Jmat(2,1) = N_eta'*X; Jmat(2,2) = N_eta'*Y;
        
        
        % One Gauss point over the element
        xi_gp =  1/3; % not used actually
        eta_gp = 1/3; % not used actually
        w_gp = 0.5;
        detJ = det(Jmat);
        Ae(ie) = w_gp*detJ;
        Jm1 = 1/detJ*[Jmat(2,2), -Jmat(1,2); -Jmat(2,1), Jmat(1,1)];
        N1_x = Jm1(1,:)*[N_xi(1) N_eta(1)]';
        N1_y = Jm1(2,:)*[N_xi(1) N_eta(1)]';
        N2_x = Jm1(1,:)*[N_xi(2) N_eta(2)]';
        N2_y = Jm1(2,:)*[N_xi(2) N_eta(2)]';
        N3_x = Jm1(1,:)*[N_xi(3) N_eta(3)]';
        N3_y = Jm1(2,:)*[N_xi(3) N_eta(3)]';
        
        B = [N1_x 0    N2_x 0    N3_x 0;
            0     N1_y 0    N2_y 0    N3_y;
            0     0    0    0    0    0;
            N1_y  N1_x N2_y N2_x N3_y N3_x];
        BD = PrD*B;
        b = m'*B;
        % Element matrix (no summation since yhere is only 1 GP)
        
        keG = B'*( Gvec'*BGs + 1i*Gvec'*BGl )*Dmu*BD*w_gp*detJ;
        keK = b'*( Kvec'*BKs + 1i*Kvec'*BKl )*b*w_gp*detJ; 
        ke = keG + keK;
        
        bvol = zeros(6,1);
        
        
        % Assembly
        for i = 1:3
            I = conne(ie,i);
            for j = 1:3
                J = conne(ie,j);
                KK(2*(I-1)+1:2*(I-1)+2,2*(J-1)+1:2*(J-1)+2) = KK(2*(I-1)+1:2*(I-1)+2,2*(J-1)+1:2*(J-1)+2) +...
                    ke(2*(i-1)+1:2*(i-1)+2,2*(j-1)+1:2*(j-1)+2);
            end
            Fvol(2*(I-1)+1:2*(I-1)+2,1) = Fvol(2*(I-1)+1:2*(I-1)+2,1) +...
                bvol(2*(i-1)+1:2*(i-1)+2,1);
        end
    end
    KKold = KK;
    F = Fvol;
    
    dofBx = 2*(bottomnodes-1)+1;
    dofBy = 2*bottomnodes;
    dofB = [dofBx,dofBy];
    % Clamped case
    for k = 1:length(dofB)
        KK(dofB(k),:) = zeros(1,2*NN);
        KK(dofB(k),dofB(k)) = 1;
    end
    F(dofB) = 0;
    clear ndof
%     %% Lateral elements at x = 0 = symmetry axis!
%     kk = 1;
%     for i = 1:ny
%         j = 1;
%         ndof(kk) = (i-1)*nx + j;
%         kk = kk + 1;
%     end
%     dofLx = 2*(ndof-1)+1;
    for k = 1:length(dofLx)
        KK(dofLx(k),:) = zeros(1,2*NN);
        KK(dofLx(k),dofLx(k)) = 1;
    end
    F(dofLx) = 0;
    clear ndof
    %% Lateral elements at x = Lx;
    %{
    kk = 1;
    for i = 1:ny
        j = nx;
        ndof(kk) = (i-1)*nx + j;
        kk = kk + 1;
    end
    dofoffRx = 2*ndof-1;
    for k = 1:length(dofoffRx)
        K(dofoffRx(k),:) = zeros(1,2*NN);
        K(dofoffRx(k),dofoffRx(k)) = 1;
    end
    clear ndof
    %}
%     %% Top elements at y = Ly;
%     % Unit displacement at beam top
%     kk = 1;
%     i = ny;
%     for j = 1:nx
%         ndof(kk) = (i-1)*nx + j;
%         kk = kk+1;
%     end
%     %     dofTx = 2*(ndof-1)+1;
%     dofTy = 2*ndof;
    % Unit displacement enforced
    for k = 1:length(dofTy)
        KK(dofTy(k),:) = zeros(1,2*NN);
        KK(dofTy(k),dofTy(k)) = 1;
    end
    F(dofTy) = amp_omevec(iome)*exp(1i*pha_omevec(iome)); 
    clear ndof
    
    
    %% Solve system
    Z = KK\F;
    
    
    %% Post process
    % Deformed coordinates
    Xdisp = Z(1:2:length(Z)-1);
    Ydisp = Z(2:2:length(Z));
    deformation = [Xdisp,Ydisp];
    deformation_amp = abs(deformation);
    
    
    % Reaction forces
    RBottom = KKold(dofB,:)*Z; 
    RBottom_amp = abs(RBottom);
    
    %         RLeft   = Kold(dofoffLx,:)*Z;
 
    
    deformationx_amp_ome(:,iome) = deformation_amp(:,1);
    deformationy_amp_ome(:,iome) = deformation_amp(:,2);
    RBottomx_amp_ome(:,iome) =  RBottom_amp(1:2*nx/2); 
    RBottomy_amp_ome(:,iome) =  RBottom_amp(2*nx/2+1:2*nx);
    
    deformationx_ome(:,iome) = deformation(:,1);
    deformationy_ome(:,iome) = deformation(:,2);
    RBottomx_ome(:,iome) =  RBottom(1:2*nx/2); 
    RBottomy_ome(:,iome) =  RBottom(2*nx/2+1:2*nx);

    
end
