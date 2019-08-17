% LTE system simulator main simulator file. Check the LTE_sim_batch files
% to check how to launch the simulator.
% [] = LTE_sim_main()
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2016 by ITC
% www.nt.tuwien.ac.at
%
% By using this simulator, you agree to the license terms stated in the license agreement included with this work.
% If you are using the simulator for your scientific work, please reference:
%
% BibTeX:
% @InProceedings{EUSIPCO2009,
%   author =        {Christian Mehlf\"uhrer and Martin Wrulich and Josep Colom Ikuno and Dagmar Bosanska and Markus Rupp},
%   title =         {Simulating the Long Term Evolution Physical Layer},
%   booktitle =     {Proc. of the 17th European Signal Processing Conference (EUSIPCO 2009)},
%   month =         aug,
%   year =          2009,
%   address =       {Glasgow, Scotland},
%   note =          {accepted for publication},
% }
% 
% ASCII
% C. Mehlf�hrer, M. Wrulich, J. C. Ikuno, D. Bosanska and M. Rupp, "Simulating the Long Term Evolution Physical Layer,"
% in Proc. of the 17th European Signal Processing Conference (EUSIPCO 2009), Aug. 2008, Glasgow, Scotland


if DEBUG_LEVEL > 0
    fprintf('LTE-A Uplink Link Level simulator v1.6\n');
    fprintf('(c) 2016, ITC, TU Wien\n');
	fprintf(' This work has been funded by Telekom Austria AG and the Christian Doppler Laboratory for Design Methodology of Signal Processing Algorithms.\n\n');
    fprintf('  By using this simulator, you agree to the license terms stated in the license agreement included with this work\n');
    fprintf('  Contains code from:\n');
    fprintf('    - pycrc (CRC checking)\n');
    fprintf('    - The Coded Modulation Library (convolutional coding & SISO decoding)\n');
    fprintf('  Convolutional coding & SISO decoding MEX files under the GNU lesser GPL license\n\n');
end


% SIMULATION TYPE
switch LTE_params.simulation_method
    case 'normal'
        LTE_UL_sim_main_single;         % normal
    case 'parallel'
        LTE_UL_sim_main_par;            % using parallel toolbox
    otherwise
        error('simulation type not supported');
end


clear tmp_results

return
