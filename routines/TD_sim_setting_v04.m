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
% function [N,Tsim,uc_th,time,dt,omevec,freqvec,indvec,fmin,fmax,df,phase,Amp] = TD_sim_setting_v04(target_mat)

function [N,Tsim,uc_th,time,dt,omevec,freqvec,indvec,fmin,fmax,df,phase,Amp] = TD_sim_setting_v04(target_mat)

tauGuK = union(target_mat.tauG,target_mat.tauK,'sorted');
fmin = .1*1/tauGuK(end); % actually we do not know these values. Optimal value 0.1 or even 0.01*(1/tau(end)) must be used.
fmax = 1/tauGuK(1);      % actually we do not know these values. Optimal value 5*(1/tau(1));  or even 10* The larger the better you get E_inf
Tmax = 1/fmin;
Tmin = 1/fmax;
numT  = 1.1; 
Tsim = numT*Tmax;
dt = Tmin/5; 
fs = 1/dt;

fNyq = fs/2;



N = Tsim/dt;
N = N - mod(N,2);

time = (0:N-1)*dt;
Tsim = N/fs;
df = 1/Tsim;
dome = 2*pi*df;

freqvec  = [0:N-1]*df;
freqvec = freqvec(2:N/2+1);

% Nyquist check
if freqvec(end) > fNyq
    disp('Abort since maximum bandlimit for generating the signal is larger than Nyquist fr.')
    return
end


S0 = .5;
Amp0 = sqrt(2*S0*dome); 
Amp = Amp0*ones(1,N/2);

num_exc_ome_aux = 15;
fgen_gen = df*unique(round(logspace(0,ceil(log10(ceil(fmax/df))),num_exc_ome_aux)));
num_exc_ome = length(fgen_gen);
for i = 1:num_exc_ome
    [minval,imin] = min(abs(fgen_gen(i) - freqvec));
    fgen_gen2(i) = freqvec(imin);
    indvec(i) = imin;
end
omevec = 2*pi*fgen_gen2;

rng(2)
phase = 2*pi*rand(1,N/2)-pi; 


for i = 1:num_exc_ome
    iex = indvec(i);
    Uc_th(:,i) =  Amp(iex)*cos(2*pi*freqvec(iex)*time + phase(iex));
end
uc_th = sum(Uc_th,2);

figure 
plot(time,uc_th)
