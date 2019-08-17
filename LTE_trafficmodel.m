function model = LTE_trafficmodel(trafficmodel_params,initialize,data_traffic_type,N_user)
% This function chooses the actual traffic model for the user
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2016 by ITC

if initialize  % That means that is called by LTE_UL_sim_main_single and the data model for each user is already assigned, it is only needed to initialize
    aPrioriPdf = [0.1,0.2,0.2,0.3,0.2,0]; % a priori pdf according to which a traffic model is picked (defined in RAN R1-070674); just use single digit numbers for this, otherwise set will not work

    if trafficmodel_params.usetraffic_model % if the traffic models shall be used
        
        
        set = [];
        for i = 1:length(aPrioriPdf)
            set = [set,i*ones(1,aPrioriPdf(i)*10)];
        end
        
        %% Random user assigment
%         index = randi([1,length(set)]); % Coment this line and discoment next one if you want "deterministic" traffic model assignemnt
        index=mod(N_user,10);
        
        
        if index==0
            index=10;
        end
        switch set(index)   % randomly pack one of the traffic models according to the aPrioriPdf
            case 1
                model = traffic_models.ftp;
            case 2
                model = traffic_models.http;
            case 3
                model = traffic_models.video;
            case 4
                model = traffic_models.voip;
            case 5
                model = traffic_models.gaming;
            case 6
                model = traffic_models.fullbuffer;
        end
    else
        model = traffic_models.fullbuffer;
    end
else
    switch data_traffic_type   % randomly pack one of the traffic models according to the aPrioriPdf
        case 'ftp'
            model = traffic_models.ftp;
        case 'http'
            model = traffic_models.http;
        case 'video'
            model = traffic_models.video;
        case 'voip'
            model = traffic_models.voip;
        case 'gaming'
            model = traffic_models.gaming;
        case 'fullbuffer'
            model = traffic_models.fullbuffer;
    end
end
end
