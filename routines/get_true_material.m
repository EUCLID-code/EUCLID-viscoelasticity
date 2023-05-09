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
% function [target_mat] =  get_true_material(matselect)

function [target_mat] =  get_true_material(matselect)


% Warning: Multiplication times 10 is here needed for the moduli to transform
% N/mm^2 (MPa) to kN/dm^2 that are the units used in the code 

switch matselect

    case 00
        G     = [5.290E+2]'*1e1;
        tauG  = [3.915E+0]'; % Devono essere messi sempre in ordine crescente
        Ginf = 5E2*1e1; % assume the same Einf for each element of the bar.
        NMeG  = length(G);
        G0 = Ginf + sum(G);
        g = G/G0; % Abaqus form

        % Volumetric
        K     = [1.097E+3]'*1e1;
        tauK  = [4.197E+0]';
        Kinf = 2E3*1e1;
        NMeK  = length(K);
        K0 = Kinf + sum(K);
        k = K/K0; % Abaqus form



    case 0
        G     = [1.019E+3 5.290E+2]'*1e1;
        tauG  = [4.824E-1 3.915E+0]'; % Devono essere messi sempre in ordine crescente
        Ginf = 5E2*1e1; % assume the same Einf for each element of the bar.
        NMeG  = length(G);
        G0 = Ginf + sum(G);
        g = G/G0; % Abaqus form

        % Volumetric
        K     = [2.366E+3 1.097E+3]'*1e1;
        tauK  = [4.570E-1 4.197E+0]';
        Kinf = 2E3*1e1;
        NMeK  = length(K);
        K0 = Kinf + sum(K);
        k = K/K0; % Abaqus form

    case 1
        G     = [1.019E+3 5.290E+2 2.011E+2]'*1e1;
        tauG  = [4.824E-1 3.915E+0 3.021E+1]'; % Devono essere messi sempre in ordine crescente
        Ginf = 5E2*1e1; % assume the same Einf for each element of the bar.
        NMeG  = length(G);
        G0 = Ginf + sum(G);
        g = G/G0; % Abaqus form

        % Volumetric
        K     = [2.366E+3 1.097E+3]'*1e1;
        tauK  = [4.570E-1 4.197E+0]';
        Kinf = 2E3*1e1;
        NMeK  = length(K);
        K0 = Kinf + sum(K);
        k = K/K0; % Abaqus form
    case 2
        G     = [1.019E+3 5.290E+2]'*1e1;
        tauG  = [4.824E-1 3.915E+0]'; % Devono essere messi sempre in ordine crescente
        Ginf = 5E2*1e1; % assume the same Einf for each element of the bar.
        NMeG  = length(G);
        G0 = Ginf + sum(G);
        g = G/G0; % Abaqus form

        % Volumetric
        K     = [2.366E+3]'*1e1;
        tauK  = [4.570E-1]';
        Kinf = 2E3*1e1;
        NMeK  = length(K);
        K0 = Kinf + sum(K);
        k = K/K0; % Abaqus form
    case 3
        G     = [7.790E+02 1.019E+3 5.290E+2 2.011E+2 9.600E+1]'*1e1;
        tauG  = [7.280E-02 4.824E-1 3.915E+0 3.021E+1 6.294E+2]'; % Must be in ascending order
        Ginf = 5E2*1e1; 
        NMeG  = length(G);
        G0 = Ginf + sum(G);
        g = G/G0; % Abaqus form

        % Volumetric
        K     = [2.242E+03 2.712E+03 2.366E+3 1.097E+3 4.601E+02]'*1e1;
        tauK  = [7.693E-03 6.344E-02 4.570E-1 4.197E+0 3.512E+01]';
        Kinf = 2E3*1e1;
        NMeK  = length(K);
        K0 = Kinf + sum(K);
        k = K/K0; % Abaqus form


    case 4
        G     = [7.790E+02 1.019E+3 5.290E+2 2.011E+2 9.600E+1]'*1e1;
        tauG  = [7.280E-02 4.824E-1 3.915E+0 3.021E+1 6.294E+2]'; % Devono essere messi sempre in ordine crescente
        Ginf = 5E2*1e1; % assume the same Einf for each element of the bar.
        NMeG  = length(G);
        G0 = Ginf + sum(G);
        g = G/G0; % Abaqus form

        % Volumetric
        K     = [2.712E+03 2.366E+3 1.097E+3 4.601E+02 1.277E-03]'*1e1;
        tauK  = [6.344E-02 4.570E-1 4.197E+0 3.512E+01 6.016E+02]';
        Kinf = 2E3*1e1;
        NMeK  = length(K);
        K0 = Kinf + sum(K);
        k = K/K0; % Abaqus form

    case 5
        G     = [7.790E+02 1.019E+3 5.290E+2 2.011E+2]'*1e1;
        tauG  = [7.280E-02 4.824E-1 3.915E+0 3.021E+1]'; % Devono essere messi sempre in ordine crescente
        Ginf = 5E2*1e1; % assume the same Einf for each element of the bar.
        NMeG  = length(G);
        G0 = Ginf + sum(G);
        g = G/G0; % Abaqus form

        % Volumetric
        K     = [2.712E+03 2.366E+3 1.097E+3 4.601E+02]'*1e1;
        tauK  = [6.344E-02 4.570E-1 4.197E+0 3.512E+01]';
        Kinf = 2E3*1e1;
        NMeK  = length(K);
        K0 = Kinf + sum(K);
        k = K/K0; % Abaqus form

end
% SET TARGET MATERIAL
target_mat.G = G;
target_mat.tauG = tauG;
target_mat.Ginf = Ginf;%

target_mat.K = K;
target_mat.tauK = tauK;
target_mat.Kinf = Kinf;
