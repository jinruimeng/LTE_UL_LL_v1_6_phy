function layer_x = LTE_common_layer_mapping(tx_mode,tx_user_symbols,nLayers,layer_x,nCodewords,cw,nAtPort)

% layer mapping for UPLINK
% 3GPP TS 36.211 V11.4.0 Section 5.3.2A
% author: Stefan Pratschner
% www.nt.tuwien.ac.at

switch tx_mode  % transmission mode used for the UE
    case 1      % single antenna transmission, 3GPP TS 36.211 V11.4.0 Section 5.3.2A.1
        layer_x = tx_user_symbols;
    
    case 4      % spatial multiplexing, 3GPP TS 36.211 V11.4.0 section 5.3.2A.2)
        switch nLayers
            case 1
                layer_x = tx_user_symbols;
            case 2
                switch nCodewords
                    case 1
                        layer_x(1,:)=tx_user_symbols(1:2:end);
                        layer_x(2,:)=tx_user_symbols(2:2:end);
                    case 2
                        layer_x(cw,:)=tx_user_symbols; 
                    otherwise
                        error('number of codewords %d not supported for % layers', nCodewords, nLayers);
                end
            case 3
                if nCodewords == 2
                    if(cw == 1)
                        layer_x(1,:)=tx_user_symbols;
                    else
                        layer_x(2,:)=tx_user_symbols(1:2:end);
                        layer_x(3,:)=tx_user_symbols(2:2:end);
                    end
                else
                    error('number of codewords %d not supported for % layers', nCodewords, nLayers);
                end
            case 4
                if nCodewords == 2
                    layer_x(2*cw-1,:)=tx_user_symbols(1:2:end);
                    layer_x(2*cw,:)  =tx_user_symbols(2:2:end);
                else
                    error('number of codewords %d not supported for % layers', nCodewords, nLayers);
                end
            otherwise
                error('number of layers not supported');
        end
    case 5      % MU MIMO
        switch nLayers
            case 1
                layer_x = tx_user_symbols;
            otherwise
                error('not implemented yet');
        end
    otherwise
        error('TX mode not supported');
end

end