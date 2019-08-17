function [w nDRMS2] = LTE_UL_cyclicShift(CSFiDCIFormat)
    switch CSFiDCIFormat            % tab. 5.5.2.1.1-1
        case 000      
            nDRMS2 = [0 6 3 9];
            w = [1 1;
                1 1;
                1 -1;
                1 -1];
        case 001
            nDRMS2 = [6 0 9 3];
            w = [1 -1;
                1 -1;
                1 1;
                1 1];
        case 010
            nDRMS2 = [3 9 6 0];
            w = [1 -1;
                1 -1;
                1 1;
                1 1];
        case 011
            nDRMS2 = [4 10 7 1];
            w = [1 1;
                1 1;
                1 1;
                1 1];
        case 100
            nDRMS2 = [2 8 5 11];
            w = [1 1;
                1 1;
                1 1;
                1 1];
        case 101
            nDRMS2 = [8 2 11 5];
            w = [1 -1;
                1 -1;
                1 -1;
                1 -1];
        case 110
            nDRMS2 = [10 4 1 7];
            w = [1 -1;
                1 -1;
                1 -1;
                1 -1];
        case 111
            nDRMS2 = [9 3 0 6];
            w = [1 1;
                1 1;
                1 -1;
                1 -1];
        otherwise
            error('DCI Format invalid');
    end
end
