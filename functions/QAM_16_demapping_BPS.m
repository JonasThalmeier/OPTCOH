function [demappedConfig_Xpol] = QAM_16_demapping_BPS(downsampledSig_Xpol)
% downsampledSig_Xpol = downsampledSig_Xpol(1:2:end); %downsample to 1SpS
% downsampledSig_Ypol = downsampledSig_Ypol(1:2:end);
division = 2; %round((abs(max(real(downsampledSig_Xpol)))+abs(min(real(downsampledSig_Xpol))))/2,2); 

% Symbol Demapping
demappedConfig_Xpol = zeros(length(downsampledSig_Xpol),1);
%demappedConfig_Ypol = zeros(length(downsampledSig_Ypol),1);

% Assuming QPSK and not accounting for noise, just a direct mapping
for i = 1:length(downsampledSig_Xpol)
    % This is a simplistic approach; real demapping would consider noise, etc.
    if real(downsampledSig_Xpol(i)) > 0
        if imag(downsampledSig_Xpol(i)) > 0
            if real(downsampledSig_Xpol(i)) > division
                if imag(downsampledSig_Xpol(i)) > division
                    demappedConfig_Xpol(i,:) = 3+1i*3;
                else
                    
                    demappedConfig_Xpol(i,:) = 3+1i*1;
                end
            else
                if imag(downsampledSig_Xpol(i)) > division
                    
                    demappedConfig_Xpol(i,:) = 1+1i*3;
                else
                    
                    demappedConfig_Xpol(i,:) = 1+1i*1;
                end
            end
        else
            if real(downsampledSig_Xpol(i)) > division
                if imag(downsampledSig_Xpol(i)) > -division
                    
                    demappedConfig_Xpol(i,:) = 3-1i*1;
                else
                    
                    demappedConfig_Xpol(i,:) = 3-1i*3;
                end
            else
                if imag(downsampledSig_Xpol(i)) > -division
                    
                    demappedConfig_Xpol(i,:) = 1-1i;
                else
                    
                    demappedConfig_Xpol(i,:) = 1-1i*3;
                end
            end
        end
    else
         if imag(downsampledSig_Xpol(i)) > 0
            if real(downsampledSig_Xpol(i)) > -division
                if imag(downsampledSig_Xpol(i)) > division
                    
                    demappedConfig_Xpol(i,:) = -1+1i*3;
                else
                    
                    demappedConfig_Xpol(i,:) = -1+1i*1;
                end
            else
                if imag(downsampledSig_Xpol(i)) > division
                    
                    demappedConfig_Xpol(i,:) = -3+1i*3;
                else
                    
                    demappedConfig_Xpol(i,:) = -3+1i*1;
                end
            end
        else
            if real(downsampledSig_Xpol(i)) > -division
                if imag(downsampledSig_Xpol(i)) > -division
                   
                    demappedConfig_Xpol(i,:) = -1-1i*1;
                else
                    
                    demappedConfig_Xpol(i,:) = -1-1i*3;
                end
            else
                if imag(downsampledSig_Xpol(i)) > -division
                    
                    demappedConfig_Xpol(i,:) = -3-1i*1;
                else
                    
                    demappedConfig_Xpol(i,:) = -3-1i*3;
                end
            end
         end 
    end
end

end