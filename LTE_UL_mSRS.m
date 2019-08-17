% see TS 36.211 v11.1.0 chapter 5.5.3.3, page 39
function [mSRS,Nb,mSRS0] = LTE_UL_mSRS(Csrs,Bsrs,Nrb)
    if Nrb >= 6 && Nrb <= 40
        N = [1 3 3 1;
            1 2 2 2;
            1 6 1 1;
            1 5 1 1;
            1 4 1 1;
            1 3 1 1;
            1 2 1 1;
            1 1 1 1];
        M = [36 12 4 4;
            32 16 8 4;
            24 4 4 4;
            20 4 4 4;
            16 4 4 4;
            12 4 4 4;
            8 4 4 4;
            4 4 4 4];
    elseif Nrb <= 60
        N = [1 2 2 3;
            1 3 2 2;
            1 2 5 1;
            1 3 3 1;
            1 2 2 2;
            1 6 1 1;
            1 5 1 1;
            1 4 1 1];
        M = [48 24 12 4;
            48 16 8 4;
            40 20 4 4;
            36 12 4 4;
            32 16 8 4;
            24 4 4 4;
            20 4 4 4;
            16 4 4 4];
    elseif Nrb <= 80
        N = [1 3 2 3;
            1 2 2 4;
            1 3 5 1;
            1 2 2 3;
            1 3 2 2;
            1 2 5 1;
            1 3 3 1;
            1 2 2 2];
        M = [72 24 12 4;
            64 32 16 4;
            60 20 4 4;
            48 24 12 4;
            48 16 8 4;
            40 20 4 4;
            36 12 4 4;
            32 16 8 4];
    elseif Nrb <= 110
        N = [1 2 2 6;
            1 3 2 4;
            1 2 2 5;
            1 3 2 3;
            1 2 2 4;
            1 3 5 1;
            1 2 2 3;
            1 3 2 2];
        M = [96 48 24 4;
            96 32 16 4;
            80 40 20 4;
            72 24 12 4;
            64 32 16 4;
            60 20 4 4;
            48 24 12 4;
            48 16 8 4];
    else
        error('number of resource blocks not supported');
    end
    
    if Csrs >= 0 && Csrs <= 7 && Bsrs >= 0 && Bsrs <= 3
        if M(Csrs,1) > Nrb
            Csrs = sum(M(:,1)>Nrb);
        end
        Nb = N(Csrs + 1,1:Bsrs + 1);
        mSRS = M(Csrs + 1,1:Bsrs + 1);
        mSRS0 = M(:,1);
    else
        error('invalid higher layer SRS configuration');
    end
end