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
% function [mesh_data] = get_mesh_data(Lx,nx,Ly,ny)

function [mesh_data] = get_mesh_data(Lx,nx,Ly,ny)

NN = nx*ny; % Total number of nodes

% Nodes and coordinates
coord = zeros(NN,2);
kk = 1;
for i = 1:ny
    y = Ly/(ny-1)*(i-1);
    for j = 1:nx
        nn = (i-1)*nx + j;
        %         x = Lx/(nx-1)*(j-1);
        x = (Lx/2)/(nx-1)*(j-1); % due to symmetry
        coord(kk,1:3) = [nn,x,y];
        kk = kk + 1;
    end
end

% Select bottom nodes only
kk = 1;
i = 1;
for j = 1:nx
    bottomnodes(kk) = (i-1)*nx + j;
    kk = kk + 1;
end

kk = 1;
for i = 2:ny-1
    for j = 2:nx-1
        nn = (i-1)*nx + j;
        intnodes(kk,1) = nn;
        kk = kk + 1;
    end
end


% Bottom dof
kk = 1;
i = 1;
for j = 1:nx
    ndof(kk) = (i-1)*nx + j;
    kk = kk + 1;
end
dofBx = 2*(ndof-1)+1;
dofBy = 2*ndof;
dofB = [dofBx,dofBy];







% Lateral dofs at x = 0 = symmetry axis!
kk = 1;
for i = 1:ny
    j = 1;
    ndof(kk) = (i-1)*nx + j;
    kk = kk + 1;
end
dofLx = 2*(ndof-1)+1;
clear ndof

% Top elements at y = Ly;
kk = 1;
i = ny;
for j = 1:nx
    ndof(kk) = (i-1)*nx + j;
    kk = kk+1;
end
dofTy = 2*ndof;
clear ndof


% Elements connectivity
conneL = zeros((nx-1)*(ny-1),3); % lower triangle
conneU = zeros((nx-1)*(ny-1),3); % upper triangle
ie = 1;
for i = 1:ny-1
    for j = 1:nx-1
        conneL(ie,1) = (i-1)*nx + j;
        conneL(ie,2) = (i-1)*nx + j + 1;
        conneL(ie,3) = i*nx + j + 1;
        conneU(ie,1) = (i-1)*nx + j;
        conneU(ie,2) = i*nx + j + 1;
        conneU(ie,3) = i*nx + j;
        ie = ie+1;
    end
end
conne = [conneL;conneU];

% Strip of Bottom elements (say at y = 0)
ie = 1;
i = 1;
for j = 1:nx-1
    conneBL(ie,1) = (i-1)*nx + j;
    conneBL(ie,2) = (i-1)*nx + j + 1;
    conneBL(ie,3) = i*nx + j + 1;
    ie = ie+1;
end
conneB = [conneBL];

% Strip of Top elements (say at y = Ly)
ie = 1;
i = ny-1;
for j = 1:nx-1
    conneTU(ie,1) = (i-1)*nx + j;
    conneTU(ie,2) = i*nx + j + 1;
    conneTU(ie,3) = i*nx + j;
    ie = ie+1;
end
conneT = conneTU;

% Strip of Left elements (say at x = 0)
ie = 1;
for i = 1:ny-1
    j = 1;
    conneLU(ie,1) = (i-1)*nx + j;
    conneLU(ie,2) = i*nx + j + 1;
    conneLU(ie,3) = i*nx + j;
    ie = ie+1;
end
conneL = conneLU;

% Strip of Right elements (say at x = Lx)
ie = 1;
for i = 1:ny-1
    j = nx-1;
    conneRL(ie,1) = (i-1)*nx + j;
    conneRL(ie,2) = (i-1)*nx + j + 1;
    conneRL(ie,3) = i*nx + j + 1;
    ie = ie+1;
end
conneR = conneRL;


% Check num of elements
if size(conne,1) ~= 2*(nx-1)*(ny-1)
    return
end

% Plot initial mesh
figure
trimesh(conne,coord(:,2),coord(:,3),zeros(NN,1),'LineStyle','-','edgecolor',[.5 .5 ,.5],'FaceAlpha',0);
view([0 0 1])
grid on
axis equal
hold off


mesh_data.nx = nx;
mesh_data.ny = ny;
mesh_data.coord = coord;
mesh_data.bottomnodes = bottomnodes;
mesh_data.intnodes = intnodes;
mesh_data.dofTy = dofTy;
mesh_data.dofLx = dofLx;
mesh_data.dofB = dofB;

mesh_data.conne = conne;
mesh_data.NN = NN;
% save ('mesh_data','mesh_data', '-v7.3');







