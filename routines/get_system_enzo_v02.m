%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Automated identification of linear viscoelastic constitutive laws    %%
%%                           with EUCLID                                 %%
%       Mechanics of Materials, Volume 181, June 2023, 104643             %  
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
%                                                                         %
% Note: this code is not meant for professional use and it contains only  %
% very few comments. Interested users can contact me for additional       %
% information and support.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AMAT,bMAT,A_total,b_total] = ...
%     get_system_enzo_v02(omevec, ...
%     NN, ...
%     coord, ...
%     conne, ...
%     bottomnodes, ...
%     hatuxc, ... % INSERT COMPLEX VALUES HERE
%     hatuyc, ... % INSERT COMPLEX VALUES HERE
%     intnodes, ...
%     tauG, ... % INSERT TRUE RELAXATION TIMES (OR LIBRARY OF RELAXATION TIMES)
%     tauK, ... % INSERT TRUE RELAXATION TIMES (OR LIBRARY OF RELAXATION TIMES)
%     hatrxc, ... % INSERT COMPLEX VALUES HERE
%     hatryc) % INSERT COMPLEX VALUES HERE
%
% This function assembles the linear system of equations (Eqs. (21) and (22) of the paper). 
function [AMAT,bMAT,A_total,b_total] = ...
    get_system_enzo_v02(omevec, ...
    NN, ...
    coord, ...
    conne, ...
    bottomnodes, ...
    hatuxc, ... % INSERT COMPLEX VALUES HERE
    hatuyc, ... % INSERT COMPLEX VALUES HERE
    intnodes, ...
    tauG, ... % INSERT TRUE RELAXATION TIMES (OR LIBRARY OF RELAXATION TIMES)
    tauK, ... % INSERT TRUE RELAXATION TIMES (OR LIBRARY OF RELAXATION TIMES)
    hatrxc, ... % INSERT COMPLEX VALUES HERE
    hatryc) % INSERT COMPLEX VALUES HERE



%% Projectors in 2D for plane strain problm
m = [1 1 1 0]';
PrD =  eye(4) - 1/3*(m*m');
PrV =1/3*(m*m');
Dmu = [2 0 0 0;
    0 2 0 0;
    0 0 2 0;
    0 0 0 1];



dof_intx = 2*(intnodes-1)+1;
dof_inty = 2*intnodes;


bottomdof_intx = 2*(bottomnodes-1)+1;
bottomdof_inty = 2*bottomnodes;

hatRx = sum(hatrxc,1).'; % Summation over the nodes since we measure the resultant only
hatRy = sum(hatryc,1).'; % Summation over the nodes since we measure the resultant only


%% FD solver
for iome = 1:length(omevec)
    Z(1:2:2*NN-1,1) = hatuxc(:,iome);
    Z(2:2:2*NN,1)   = hatuyc(:,iome);
    
    omega = omevec(iome);
    
    BGs = [1;omega.^2.*tauG.^2./(1+omega.^2.*tauG.^2)];
    BGl = [0;omega.*tauG./(1+omega^2.*tauG.^2)];
    
    BKs = [1;omega.^2*tauK.^2./(1+omega.^2.*tauK.^2)];
    BKl = [0;omega.*tauK./(1+omega.^2.*tauK.^2)];
    
    AA  = zeros(2*NN,length(tauG) + length(tauK) + 2);
    for ie = 1:size(conne,1)
        vrtx = conne(ie,1:3); % vertexes array
        X = [coord(vrtx(1),2),coord(vrtx(2),2),coord(vrtx(3),2)]';
        Y = [coord(vrtx(1),3),coord(vrtx(2),3),coord(vrtx(3),3)]';
        N_xi  = [-1 1 0]';
        N_eta = [-1 0 1]';
        Jmat(1,1) = N_xi'*X; Jmat(1,2) = N_xi'*Y;
        Jmat(2,1) = N_eta'*X; Jmat(2,2) = N_eta'*Y;
        dofglob = reshape([2*vrtx-1; 2*vrtx],6,1);
        
        % One Gauss point over the element
%         xi_gp =  1/3; % not used actually
%         eta_gp = 1/3; % not used actually
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
        
        vloce = Z(dofglob);
        
        AeG = ((B'*Dmu*BD)*vloce) * (BGs + 1i*BGl).' * w_gp*detJ;
        AeK = ((b'*b)     *vloce) * (BKs + 1i*BKl).' * w_gp*detJ;
        
        
        AAe = [AeG,AeK];
        %         Assembly
        for i = 1:3
            I = conne(ie,i);
            AA(2*(I-1)+1:2*(I-1)+2,:) = AA(2*(I-1)+1:2*(I-1)+2,:) +...
                AAe(2*(i-1)+1:2*(i-1)+2,:);
        end
    end
    
    
    AA_intx = AA(dof_intx,:);
    AA_inty = AA(dof_inty,:);
    
    AA_bottomx = AA(bottomdof_intx,:);
    AA_bottomy = AA(bottomdof_inty,:);
    
    AA_intx_ome{iome} = AA_intx;
    AA_inty_ome{iome} = AA_inty;
    
    
    AA_bottomx_ome{iome} = sum(AA_bottomx,1);
    AA_bottomy_ome{iome} = sum(AA_bottomy,1);
    
    
    
end

A_intx = cat(1,cell2mat(AA_intx_ome'));
A_inty = cat(1,cell2mat(AA_inty_ome'));


A_bottomx = cat(1,cell2mat(AA_bottomx_ome'));
A_bottomy = cat(1,cell2mat(AA_bottomy_ome'));

Aint = [A_intx;
    A_inty];


lambda_r = 1;
Abtm = lambda_r*[A_bottomx;
    A_bottomy];
A_total = [Aint;Abtm];

bint = zeros(2*length(omevec)*length(intnodes),1);
bbtm = lambda_r*[hatRx;hatRy];
b_total = [bint;bbtm];

sci = 1e0;
scr = 1e0;
AMAT = [scr*real(A_total);sci*imag(A_total)];
bMAT = [scr*real(b_total);sci*imag(b_total)];
return