classdef TR36873_Fading_3D_Channel_LL < channel_gain_wrappers.TR36873_Fading_3D_Channel
    methods
        function obj = TR36873_Fading_3D_Channel_LL(LTE_config)
            % call superclass (system level) constructor
            obj@channel_gain_wrappers.TR36873_Fading_3D_Channel(LTE_config);

            obj.bandwidth = LTE_config.bandwidth;
            obj.Nsc = obj.resourceBlock/obj.subcarrierSpacing;
            
            if(obj.bandwidth == 1.4e6)
                obj.Nrb = 6*6;
            else
                obj.Nrb = (obj.bandwidth*0.9) / obj.subcarrierSpacing / 2;
            end
            obj.Ntot = obj.Nsc*obj.Nrb;
            obj.Nfft =  2^ceil(log2(obj.Ntot));
            obj.Tb = 1/obj.subcarrierSpacing;
            obj.fs = obj.subcarrierSpacing*obj.Nfft;
        end
    end
end