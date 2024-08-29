clear all
close all
clc

%% CONSTANTS
R = 8314;

%% GRID DETAILS

z_left = -0.5;
z_right = 0.5;
ngrid = 100;
z = linspace(z_left, z_right, ngrid);

%% INITIAL GUESS

species=["H2" "O2" "N2" "H2O"];
% reducing and parsing thermo model
ckm = CKM;
mechanism_folder = "GRI30";

ckm.thRed("./" + mechanism_folder + "/thermo.txt",species)
ckm.trRed("./" + mechanism_folder + "/trans.txt",species)

% models
load('./transport_models/red.mat');
load('./thermo_models/red.mat');

data_thermo = dataTh;
data_trans = dataTr;

for k = 1:length(species)

    sp(k).name = cell2mat(data_thermo((k-1)*4+1,1));
    species(k) = sp(k).name;
    W_k(k) = cell2mat(data_thermo((k-1)*4+1, 5));
    sp(k).M = W_k(k);
    sp(k).T_range = [300.000  1000.000  5000.000];
    sp(k).Coes = [];
    temp = [];

    for i = 1:3
        temp = [temp data_thermo((k-1)*4+1+i, :)];
    end
    for i = 1:13
        sp(k).Coes(i) = str2num(cell2mat(temp(i)));
    end


end

% Initialize a cell array to hold the reordered data_trans
reordered_data_trans = cell(size(data_trans, 1), size(data_trans, 2));

% Create a map for quick lookup
[~, idx_map] = ismember(species, data_trans(:,1));

% Reorder data_trans
for i = 1:length(idx_map)
    if idx_map(i) > 0
        reordered_data_trans(i, :) = data_trans(idx_map(i), :);
    else
        % Handle case where molecule is not found in data_trans
        reordered_data_trans(i, :) = repmat({''}, 1, size(data_trans, 2));
    end
end

data_trans = reordered_data_trans;

for k = 1:length(species)
    sp(k).TpCoes(1) = str2num(cell2mat(data_trans(k, 3)));
    sp(k).TpCoes(2) = str2num(cell2mat(data_trans(k, 4)));
end

% The format corresponds to the CHEMKIN-II file formatting.
% Temperature range and Polynomial coefficients from Chemkin therm.dat file
% *.Coes[1:7]: a1 to a7 for upper temperature interval (*.T_range(2) to *.T_range(3))
% *.Coes[8:14]: a1 to a7 for lower temperature interval (*.T_range(1) to *.T_range(2))
% CH4
% CH4.M = 0.01604; % Kg/mol
% CH4.T_range = [300,1000,5000];
% CH4.Coes = [1.68347900e+00,1.02372400e-02,-3.87512900e-06,6.78558500e-10,-4.50342300e-14,...
%              -1.00807900e+04,9.62339500e+00,7.78741500e-01,1.74766800e-02,-2.78340900e-05,...
%              3.04970800e-08,-1.22393100e-11,-9.82522900e+03,1.37221900e+01];
% CH4.TpCoes = [141.400,3.746]; % [epsilon/kB,sigma]


u0 = 30;
T0 = 300;
Y0 = [0.3 0.0 0.0 0.7];
p0 = 101325;

u = ones(1, ngrid);
rho = ones(1, ngrid);
T = ones(1, ngrid);
for k = 1:length(W_k)
    Y(k, :) = ones(1, ngrid);
end
F_u = ones(1, ngrid);
F_T = ones(1, ngrid);
for k = 1:length(W_k)
    F_Y(k,:) = ones(1,ngrid);
end

u = u0*u;
T = T0*T;
for k = 1:length(W_k)
    Y(k, :) = Y0(k)*ones(1,ngrid);
end
for i = 1:length(z)
    W(i) = 1/sum(Y(:,i)./W_k');
    for k = 1:length(W_k)
        X(k, i) = W(i)*Y(k,i)/W_k(k);
    end
end
rho = p0./(R.*T).*W;

%% PROPERTIES



for i = 1:length(z)
    
    cp_mixture(i) = faux.MixCp_CK(sp, X(:,i), Y(:,i),T(i));

    lambda_mixture(i) = faux.MixLambda_CK(sp,X(:,i),T(i));

    for k = 1:length(species)
        h(k, i) = faux.Ent_CK(sp(k),T(i));
    end

end

[diff_species] = importfile_diffspecies("Djk.xlsx", "Sheet1", [1, Inf]);

[~, idx] = ismember(species, diff_species);


Djk = importfile("Djk.xlsx", "Sheet1", [1, Inf]);

Djk = Djk(idx, idx);


%% FLUXES

for k = 1:length(species)
    for i = 1:length(z)
    
        if i == 1

            DX(k, i) = (X(k, i+1) - X(k, i))/(z(i+1)-z(i));

        elseif i == length(z)
            DX(k, i) = (X(k, i) - X(k, i-1))/(z(i)-z(i-1));
            
        else
            DX(k, i) = (X(k, i+1) - X(k, i-1))/(z(i+1)-z(i-1));
            
        end
    
    end
end

for i = 1:length(z)
    
    for k = 1:length(species)

        indexes = [1:k-1 k+1:length(species)];

        Dkm(k, i) = (1 - Y(k, i))/sum(X(indexes, i)./Djk(indexes, k));

        jkstar(k, i) = -rho(i)*W_k(k)/W(i)*Dkm(k, i)*DX(k, i);

        jk(k, i) = jkstar(k, i) - Y(k, i)*sum(jkstar(:, i));

    end
end


for k = 1:length(species)
    for i = 1:length(z)
    
        if i == 1

            jplus05(k, i) = (jk(k, i) + jk(k, i+1))/2;
            jminus05(k, i) = jk(k, i);

        elseif i == length(z)
            jplus05(k, i) = jk(k, i);
            jminus05(k, i) = (jk(k, i) + jk(k, i-1))/2;
            
        else
            jplus05(k, i) = (jk(k, i) + jk(k, i+1))/2;
            jminus05(k, i) = (jk(k, i) + jk(k, i-1))/2;
            
        end
    
    end
end

%% ARRHENIUS

for i = 1:length(z)
    for k = 1:length(species)
        omega(k, i) = 0;
    end
end

%% RESIDUALS COMPUTATION



[F_u, F_Y, F_T] = faux.compute_residuals(z, u,...
                rho, Y, T, F_u, F_Y, F_T, jk, jplus05, jminus05, omega, W_k, ...
                lambda_mixture, cp_mixture, h, T0, Y0);




