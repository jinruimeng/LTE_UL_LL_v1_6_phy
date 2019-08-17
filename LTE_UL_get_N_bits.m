function N = LTE_UL_get_N_bits(CQI)
           switch CQI   % randomly pack one of the traffic models according to the aPrioriPdf
        case 1,
            N = 24;
        case 2,
            N = 24;
        case 3,
            N = 24;
        case 4,
            N = 56;
        case 5,
            N = 88;
        case 6,
            N = 128;
        case 7,
            N = 168;
        case 8,
            N = 232;
        case 9,
            N = 296;
        case 10,
            N = 336;
               case 11,
            N = 416;
               case 12,
            N = 488;
               case 13,
            N = 576;
               case 14,
            N = 648;
               otherwise
                   N = 696;
           end
       end