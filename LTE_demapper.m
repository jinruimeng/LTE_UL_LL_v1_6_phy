function LLR_SD_C = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement,varargin)
% Soft Sphere decoder.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at

k = 2.^M;
for layerIdx = 1:nLayers
    rx_layer           = rx_layer_x(layerIdx,:);
    rx_layer_repmat    = rx_layer(ones(1,k(layerIdx)),:);
    const_layer        = symbol_alphabet(layerIdx,1:k(layerIdx)).';
    const_layer_repmat = const_layer(:,ones(1,size(rx_layer_x,2)));
    [C,symbols_ZF(layerIdx,:)] = min(abs(rx_layer_repmat-const_layer_repmat).^2,[],1);      %MINIMUM DISTANCE
end

s_alph = [];
for mm=1:nLayers
    s_alph = [s_alph; symbol_alphabet(mm,symbols_ZF(mm,:))];
end
% dist_ZF = abs(rx_layer_x-s_alph).^2;   % distance to the ZF solution initial value for SS Decoder


if isempty(varargin) || ~strcmp(varargin{1},'ZF')% && ~strcmp(varargin{1},'MMSE'))
    % Soft Sphere Decoder
    if imag(Hg) == 0
        Hg = complex(Hg);
    end    
    %MEX FILE  output: [  M (bits)   x  subcarrier]
    
%     LLR_SD_C = LTE_rx_soft_sd2(Hg,rx_layer_x,dist_ZF,int32(symbols_ZF),int32(M),symbol_alphabet.',bittable)./noise_enhancement;
    
    [Q,R] = qr(Hg);
    siz = size(R,2);
    if (siz < size(R,1)) % chop off unnecessary data
         R = R(1:siz,:);
         Q = Q(:,1:siz);
    end
    stemp = Q'*rx_layer_x;
    
    dist_ZF = sum(abs(stemp-R*s_alph).^2,1);   % distance to the ZF solution initial value for SS Decoder
    if (imag(R) == 0) % The SSD needs a complex matrix, or else the MEX version of it will crash
        R = complex(R);
    end
    
    LLR_SD_C = LTE_rx_soft_sd2(R,stemp,dist_ZF,int32(symbols_ZF),int32(M),symbol_alphabet.',bittable)./noise_enhancement;
    
%     "[LLR] = soft_sd(R,s,dist_ZF,symbols_ZF,symbol_alphabet,bittable) ... soft Sphere Decoder\n\n"
% 				 "  R ... upper triangular matrix obtained from the QR decomposition of the channel H (complex)\n"
% 				 "  s ... received symbol vector, s=Q^H*y (nR x nSym) (complex)\n"
%                  "  dist_ZF ... Distance of the zero forcing solution (real)\n"
%                  "  symbols_ZF_i ... indices to symbols of the ZF solution (nT x nSym) (real integer)\n"
%                  "  M ... number of bits in the corresponding layer (1 x nR) (real)\n"
% 				 "  symbol_alphabet ... for the demapping (2^M_max x nT) (complex)\n"
% 				 "  bittable ... matrix containing the bits according to the symbol_alphabet (M x 2^M) (logical)\n"
% 				 "  LLR  ... max-log-MAP approximation of the LLR values (M*nR) (real)\n\n");
    

% elseif strcmp(varargin{1},'MMSE') % independent detection of layers for
% MMSE; notice: this requires consideration of inter-layer interference,
% which is not yet done
%     
%     dist_ZF = abs(rx_layer_x-s_alph).^2;   % distance to the ZF solution initial value for SS Decoder
%     
%     
%     LLR_SD_C = zeros(size(noise_enhancement));
%     summe = 1;
%     if imag(rx_layer_x) == 0
%         rx_layer_x = complex(rx_layer_x);
% %         rx_layer_x
%     end    
%     for i = 1:nLayers
%         Ht = Hg(i,i);
%         if imag(Ht) == 0
%             Ht = complex(Ht);
%         end
%         LLR_SD_C(summe:summe+M(i)-1,:) = LTE_rx_soft_sd2(Ht,rx_layer_x(i,:),dist_ZF(i,:),int32(symbols_ZF(i,:)),int32(M(i)),symbol_alphabet(i,:).',bittable(summe:summe+M(i)-1,1:2^M(i)))./noise_enhancement(summe:summe+M(i)-1,:);
%         summe = summe+M(i);
%     end
else
%% If you use ZF and there is no interference (inter carrier, imperfect channel knowledge,...) use this demapper to increase the speed (it demaps the layers independently)
    % Soft Sphere Decoder
    % if imag(Hg) == 0
    %     Hg = complex(Hg);
    % end
    
    dist_ZF = abs(rx_layer_x-s_alph).^2;   % distance to the ZF solution initial value for SS Decoder

    LLR_SD_C = zeros(size(noise_enhancement));
    summe = 1;
    if imag(rx_layer_x) == 0
        rx_layer_x = complex(rx_layer_x);
%         rx_layer_x
    end
    for i = 1:nLayers
        LLR_SD_C(summe:summe+M(i)-1,:) = LTE_rx_soft_sd2(1+eps*1i,rx_layer_x(i,:),dist_ZF(i,:),int32(symbols_ZF(i,:)),int32(M(i)),symbol_alphabet(i,:).',bittable(summe:summe+M(i)-1,1:2^M(i)))./noise_enhancement(summe:summe+M(i)-1,:);
        summe = summe+M(i);
    end
end
