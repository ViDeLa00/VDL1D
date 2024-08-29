classdef faux
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Property1
    end

    methods (Static)

        function H_ck = Ent_CK(Sp,T)
            % A small function to evaluate thermodynamic data from chemkin database
            % INPUT:
            %     -- Sp: Name of the species
            %     -- T: Temperature at which the enthalpy will be evaluated (K)
            % OUTPUT:
            %     -- H_ck: The specific enthalpy (options for different units below)

            % Universal gas constant
            R = 8.314510e7 ; % Unit: Chemkin unit, ergs/(mole*K), = 8.314510 J/(mol*K)
            % 1 erg = 1e-10 kilojoule = 1e-7 joule

            tH = 0.0;
            if (T > Sp.T_range(2))
                for i = 1:5
                    tH = tH + Sp.Coes(i).*T.^i./i;  % For upper temperature interval
                end
                tH = tH + Sp.Coes(6);
            else
                for i = 1:5
                    tH = tH + Sp.Coes(i+7).*T.^i./i; % Fot lower temperature interval
                end
                tH = tH + Sp.Coes(13);
            end

            H_ergs = tH .* R; % Unit: ergs/(mole)

            H_J = H_ergs .* 1e-07; % Unit: J/(mole)

            H_Jkg = H_J ./(1000.*Sp.M); % Unit: kJ/(kg)

            H_cal = H_ergs .* 1e-010*238.846; % Unit: cal/(mole), 1 kj/mol = 238.846 cal/mol;

            % H can be returned in different units,
            % for example, change below to H_ck = H_J; to return H in J/(mole)
            % H_ck = H_J; % Return H in J/(mole)
            H_ck = H_Jkg; % Return H in J/kg
        end

        function lambda_ck = Lambda_CK(Sp,T)
            % A small function to evaluate thermodynamic data from chemkin database
            % INPUT:
            %     -- Sp: Name of the species
            %     -- T: Temperature at which the Cp will be evaluated
            % OUTPUT:
            %     -- Cp_ck: The specific heat capacity at constant pressure for Sp at
            %               temperature T.

            % Universal gas constant
            R = 8.314510 ; % Unit: J/(mol*K)
            Cp = faux.CP_CK(Sp,T);  % Unit: J/(mol*K)
            Cv = Cp - R; % Specific heat under constant volume, % Unit: J/(mol*K)
            gama = Cp./Cv;

            Cv = Cv./Sp.M; % Specific heat under constant volume, % Unit: Kg/(mol*K)
            mu = faux.Mu_CK(Sp,T);

            lambda_ck = 0.25.*(9.*gama-5).*mu.*Cv; % Return lambda in W/(m*K)
        end

        function mixCp = MixCp_CK(Mix,X, Y,T)
            % A small function to evaluate thermodynamic data from chemkin database
            % INPUT:
            %     -- Mix: Names of the species in the mixture, e.g. [CO2,H2O,O2,N2]
            %     -- X: mole fractions of each species, e.g. [0.2,0.3,0.1,0.4]
            %     -- T: Temperature at which the viscosity will be evaluated, in K
            % OUTPUT:
            %     -- mixLambda: The specific heat for Mix at temperature T.

            tmixCp = 0.0;

            N = length(Mix);
            mixM = 0.0;
            for a = 1 : N
                tmixCp  = tmixCp + Y(a).*faux.CP_CK(Mix(a),T);
                mixM = mixM + X(a).*Mix(a).M;
            end

            mixCp = tmixCp./(0.001*mixM); % Return Cp of mixture in J/(kg*K)
            % mixCp = tmixCp; % Return Cp of mixture in J/(mole*K)

        end

        function mixLambda = MixLambda_CK(Mix,X,T)
            % A small function to evaluate thermodynamic data from chemkin database
            % INPUT:
            %     -- Mix: Names of the species in the mixture, e.g. [CO2,H2O,O2,N2]
            %     -- X: mole fractions of each species, e.g. [0.2,0.3,0.1,0.4]
            %     -- T: Temperature at which the viscosity will be evaluated, in K
            % OUTPUT:
            %     -- mixLambda: The viscosity for Mix at temperature T.

            tLambda_times = 0.0;
            tLambda_divide = 0.0;
            N = length(Mix);
            for a = 1 : N

                tLambda_times  = tLambda_times + X(a).*faux.Lambda_CK(Mix(a),T);
                tLambda_divide  = tLambda_divide + X(a)./faux.Lambda_CK(Mix(a),T);
            end

            mixLambda = 0.5.*(tLambda_times + 1./tLambda_divide); % Return viscosity in Pa.S, Kg/(s m)


        end

        function Mu_ck = Mu_CK(Sp,T)
            % A small function to evaluate thermodynamic data from chemkin database
            % INPUT:
            %     -- Sp: Name of the species
            %     -- T: Temperature at which the viscosity will be evaluated
            % OUTPUT:
            %     -- Mu_ck: The viscosity at constant pressure for Sp at
            %               temperature T.
            Tep = Sp.TpCoes(1);
            Sigma = Sp.TpCoes(2);
            Omega = 1.147.*(T./Tep).^(-0.145) + (T./Tep + 0.5)^(-2);

            tMu_ck = 2.6693e-05.*(Sp.M.*1000.*T)^0.5/(Sigma.^(2).*Omega);

            % Cp can be returned in different units,
            % for example, change below to Cp_ck = Cp_J; to return Cp in J/(mole*K)
            Mu_ck = 0.1*tMu_ck; % Return Cp in Pa.S, Kg/(s m)


        end

        function Cp_ck = CP_CK(Sp,T)
            % A small function to evaluate thermodynamic data from chemkin database
            % INPUT:
            %     -- Sp: Name of the species
            %     -- T: Temperature at which the Cp will be evaluated
            % OUTPUT:
            %     -- Cp_ck: The specific heat capacity at constant pressure for Sp at
            %               temperature T.

            % Universal gas constant
            % R = 8.314510e7 ; % Unit: Chemkin unit, ergs/(mole*K), = 8.314510 J/(mol*K)
            % 1 erg = 1e-10 kilojoule = 1e-7 joule
            R = 8.314;

            tCp = 0.0;
            if (T > Sp.T_range(2))
                for i = 1:5
                    tCp = tCp + Sp.Coes(i).*T.^(i-1);  % For upper temperature interval
                end
            else
                for i = 1:5
                    tCp = tCp + Sp.Coes(i+7).*T.^(i-1); % Fot lower temperature interval
                end
            end

            % Cp_ergs = tCp .* R; % Unit: ergs/(mole*K)

            % Cp_J = Cp_ergs .* 1e-07; % Unit: J/(mole*K)
            Cp_J = R*tCp;
            % Cp_kJkg = Cp_J ./(0.001.*Sp.M); % Unit: kJ/(kg*K)

            % Cp_cal = Cp_ergs .* 1e-010*238.846; % Unit: cal/(mole*K), 1 kj/mol = 238.846 cal/mol;

            % Cp can be returned in different units,
            % for example, change below to Cp_ck = Cp_J; to return Cp in J/(mole*K)
            Cp_ck = Cp_J; % Return Cp in J/(mole K)
            % Cp_ck = Cp_Jkg; % Return Cp in J/(kg K)
        end

        function [F_u, F_Y, F_T] = compute_residuals(z, u,...
                rho, Y, T, F_u, F_Y, F_T, jk, jplus05, jminus05, omega, W_k, ...
                lambda, cp, h, T0, Y0)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here



            for i = 1:length(z)

                % CONTINUITY RESIDUAL

                % if z(i) > z_fixed
                %     F_u(i) = -(rho(i)*u(i) - rho(i-1)*u(i-1))/(z(i) - z(i-1));
                % elseif z(i) == z_fixed
                %     F_u(i) = T(i) - T_fixed;
                % else
                %     F_u(i) = -(rho(i+1)*u(i+1) - rho(i)*u(i))/(z(i+1) - z(i));
                % end
                if i == length(z)
                    F_u(i) = (rho(i)*u(i) - rho(i-1)*u(i-1))/(z(i) - z(i-1));
                elseif i == 1
                    F_u(i) = (rho(i)*u(i) - rho(i+1)*u(i+1))/(z(i) - z(i+1));
                else
                    F_u(i) = (rho(i)*u(i) - rho(i+1)*u(i+1))/(z(i) - z(i+1));
                end


                % SPECIES EQUATION RESIDUAL



                for k = 1:length(W_k)

                    if i == length(z)
                        F_Y(k, i) = rho(i)*u(i)*Y(k, i) - jminus05(k, i);
                    elseif i == 1
                        F_Y(k, i) = Y(k, i) - Y0(k);
                    else
                        zplus05 = (z(i) + z(i+1))/2;
                        zminus05 = (z(i) + z(i-1))/2;
                        F_Y(k, i) = -rho(i)*u(i)*((Y(k,i) - Y(k,i-1))/(z(i) - z(i-1))) - ...
                            (jplus05(k, i) - jminus05(k, i))/(zplus05 - zminus05) + omega(k,i)*W_k(k);
                    end
                end

                % ENERGY EQUATION RESIDUAL



                if i == length(z)
                    F_T(i) = T(i) - T(i-1);
                elseif i ==1
                    F_T(i) = T(i) - T0;
                else
                    lambda_iplus05 = (lambda(i) + lambda(i+1))/2;

                    lambda_iminus05 = (lambda(i) + lambda(i-1))/2;

                    F_T(i) = -rho(i)*cp(i)*u(i)*((T(i)-T(i-1))/(z(i) - z(i-1))) + ...
                        (lambda_iplus05*((T(i+1) - T(i))/(z(i+1) - z(i))) - lambda_iminus05*...
                        (T(i) - T(i-1))/(z(i) - z(i-1)))/((z(i+1) - z(i-1))/2) - ...
                        sum(jk(:, i).*(h(:, i) - h(:, i-1))/(z(i) - z(i-1)));
                end

            end
        end

    end
end