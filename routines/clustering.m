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
% function [S,cost_vec] = clustering(tau_G_activated,G_activated,tau_K_activated,K_activated,Ginf_dsc,Kinf_dsc,...
%     omevec,NN,coord,conne,bottomnodes,hatuxc,hatuyc,intnodes,hatrxc,hatryc)

function [S,cost_vec] = clustering(tau_G_activated,G_activated,tau_K_activated,K_activated,Ginf_dsc,Kinf_dsc,...
    omevec,NN,coord,conne,bottomnodes,hatuxc,hatuyc,intnodes,hatrxc,hatryc)

indxGact = find(tau_G_activated>0);
tau_Gact = tau_G_activated(indxGact);
Gact = G_activated(indxGact);

indxKact = find(tau_K_activated>0);
tau_Kact = tau_K_activated(indxKact);
Kact = K_activated(indxKact);


cost_vec = [];
numclust = 1; 

sumdist = ones(numclust,1);
err = 1;
while abs(err) > 1e-7
    [idxG,CG,sumdistG,distG] = kmeans(log10(tau_Gact),numclust);
    tau_G = 10.^CG;
    [idxK,CK,sumdistK,distK] = kmeans(log10(tau_Kact),numclust);
    tau_K = 10.^CK;
    G = zeros(numclust,1);
    for i = 1:length(Gact)
        G(idxG(i),1) = G(idxG(i),1) + Gact(i,1);
    end
    [tau_G,iG] = sort(tau_G);
    G = G(iG);
    Gvec = [Ginf_dsc;G];

    K = zeros(numclust,1);
    for i = 1:length(Kact)
        K(idxK(i),1) = K(idxK(i),1) + Kact(i,1);
    end
    [tau_K,iK] = sort(tau_K);
    K = K(iK);
    Kvec = [Kinf_dsc;K];

    S(numclust).Ginf = Ginf_dsc;
    S(numclust).G = G;
    S(numclust).tau_G = tau_G;

    S(numclust).Kinf = Kinf_dsc;
    S(numclust).K = K;
    S(numclust).tau_K = tau_K;

    S(numclust).param = [[Gvec;Kvec],[NaN;tau_G;NaN;tau_K]];


    [A_total,b_total,Acmplx,bcmplx] = get_system_enzo_v02(omevec,NN,coord,conne,bottomnodes,hatuxc,hatuyc,intnodes,tau_G,tau_K,hatrxc,hatryc);

    cost = A_total*[Gvec;Kvec]- b_total;
    err = (norm(cost))^2/length(cost); %OK!
    cost_vec(numclust) = err;
    if numclust == min(length(tau_Gact),length(tau_Kact))
        break
    end
    numclust = numclust + 1;
end
end