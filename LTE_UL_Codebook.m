% see TS 36.211 v11.1.0 chapter 5.3.3A.2, pages 18-19
function [U] = LTE_UL_Codebook(rank, index, Nt)
    % check if the codebook index is valid!!!
    U = 0;
    switch Nt
        case 1
            U = 1;
        case 2
            switch rank
                case 1
                    U = [1 1;
                        1 -1;
                        1 1i;
                        1 -1i;
                        1 0;
                        0 1];
                    
                    U = 1/sqrt(2)*transpose(U(mod(index,6)+1,:));
                case 2
                    U = 1/sqrt(2)*eye(2);
%                     U = eye(2);
            end
        case 4
            switch rank
                case 1
                    U = [1 1 1 -1;
                        1 1 1i 1i;                    
                        1 1 -1 1;                    
                        1 1 -1i -1i;                    
                        1 1i 1 1i;                    
                        1 1i 1i 1;                    
                        1 1i -1 -1i;                    
                        1 1i -1i -1;                    
                        1 -1 1 1;                    
                        1 -1 1i -1i;
                        1 -1 -1 -1;
                        1 -1 -1i -1i;
                        1 -1i 1 -1i;
                        1 -1i 1i -1;
                        1 -1i -1 1i;
                        1 -1i -1i 1;
                        1 0 1 0;
                        1 0 -1 0;
                        1 0 1i 0;
                        1 0 -1i 0;
                        0 1 0 1;
                        0 1 0 -1;
                        0 1 0 1i;
                        0 1 0 -1i];

                    U = 1/2*transpose(U(mod(index,24)+1,:));
                case 2
                    U = [1 1 0 0; 0 0 1 -1i;
                        1 1 0 0; 0 0 1 1i;
                        1 -1i 0 0; 0 0 1 1;
                        1 -1i 0 0; 0 0 1 -1;
                        1 -1 0 0; 0 0 1 -1i;
                        1 -1 0 0; 0 0 1 1i;
                        1 1i 0 0; 0 0 1 1;
                        1 1i 0 0; 0 0 1 -1;
                        1 0 1 0; 0 1 0 1;
                        1 0 1 0; 0 1 0 -1;
                        1 0 -1 0; 0 1 0 1;
                        1 0 -1 0; 0 1 0 -1;
                        1 0 0 1; 0 1 1 0;
                        1 0 0 1; 0 1 -1 0;
                        1 0 0 -1; 0 1 1 0;
                        1 0 0 -1; 0 1 -1 0];

                    index = 2*mod(index,16);
                    U = 1/2*transpose(U(index+1:index+2,:)); 
                case 3
                    U = [1 1 0 0; 0 0 1 0; 0 0 0 1;
                        1 -1 0 0; 0 0 1 0; 0 0 0 1;
                        1 0 1 0; 0 1 0 0; 0 0 0 1;
                        1 0 -1 0; 0 1 0 0; 0 0 0 1;
                        1 0 0 1; 0 1 0 0; 0 0 1 0;
                        1 0 0 -1; 0 1 0 0; 0 0 1 0;
                        0 1 1 0; 1 0 0 0; 0 0 0 1;
                        0 1 -1 0; 1 0 0 0; 0 0 0 1;
                        0 1 0 1; 1 0 0 0; 0 0 1 0;
                        0 1 0 -1; 1 0 0 0; 0 0 1 0;
                        0 0 1 1; 1 0 0 0; 0 1 0 0;
                        0 0 1 -1; 1 0 0 0; 0 1 0 0];

                    index = 3*mod(index,12);
                    U = 1/2*transpose(U(index+1:index+3,:)); 
                case 4
                    U = 1/2*eye(4);
%                     U = eye(4);
            end
    end            
end
