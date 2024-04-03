classdef Conversion
    %CONVERSION Utility class for unit conversions
    properties(Constant)
        q = 1.60217662e-19; % elementary charge
        c = 299792458; % light speed in vacuum
        h = 6.62607004e-34; % Planck's constant
        k_B = 1.3806485e-23; % Boltzmann's constant
    end
    
    methods(Static)
        function Y = lin2dBm(X) % convert power level from Watt to dBm
            Y = 10*log10(X/1e-3);
        end
        
        function Y = dBm2lin(X)  % convert power level from dBm to Watt
            Y = 10.^(X/10)*1e-3;
        end
        
        function Y = coh_time2length(X) % convert coherence time to coherence lenght
            Y = Conversion.c*X;
        end
        
        function Y = pwr2Np(pwr,frek,BW) % convert power to number of photons
            Y = pwr/(Conversion.h*frek*BW);
        end
        
        function Y = eta2R(eta,frek) % convert quantum efficiency to reposnsivity
            Y = (eta*Conversion.q)/(Conversion.h*frek);
        end
        
        function Y = lambda2frek(X) % convert wavelength to frequency
            Y = Conversion.c/X;
        end
        
        function Y = frek2lambda(X) % convert frequency to wavelength
            Y = Conversion.c/X;
        end
        
        function Y = pwr2photon_flux(pwr,frek) % convert power to mean photon flux [photons/s]
            Y = pwr/(Conversion.h*frek);
        end
        
        function Y = LW2var(LW,Fs) % convert laser linewidth to the corresponding varaince
            Y =  2*pi*LW*(1/Fs);
        end
        
        function Y = dB2lin(X) % conversion dB to linear
            Y = 10^(X./10);
        end
        
        function Y = lin2dB(X) % conversion lin to dB
            Y = 10*log10(X);
        end
    end
end

