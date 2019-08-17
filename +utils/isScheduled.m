function [scheduled] = isScheduled(UE_MCS_and_scheduling_info, bb, uu, connection_table)
% checks if uu is scheduled at bb

scheduled = false;

for bbb = bb
    if  utils.isUEattached( connection_table, bbb, uu )
        u_local = utils.globalToLocalUser( connection_table, bbb, uu );
        scheduled = UE_MCS_and_scheduling_info(bbb, u_local).assigned_RBs > 0;
        
        if scheduled
            return;
        end
    end
end

end