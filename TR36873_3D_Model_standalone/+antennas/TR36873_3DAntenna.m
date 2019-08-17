classdef TR36873_3DAntenna < antennas.antenna
    % 3D antenna pattern according to TR 36.873 Table 7.1-1
    % (c) Fjolla Ademaj, Martin Taranetz, ITC 2016
   
    methods
        function obj = TR36873_3DAntenna(LTE_config)
            obj.antenna_type       = 'TR36.873 3D antenna';
            obj.max_antenna_gain   = LTE_config.antenna.max_antenna_gain;
            obj.pattern_is_3D      = true;       
        end
        
        % Print some information of the antenna
        function print(obj)
              fprintf('%s Antenna with maximum gain %f \n',obj.antenna_type,obj.max_antenna_gain);
        end
        
        % Returns antenna gain in [dB] as a function of theta and phi
        function antenna_gain = gain(obj,theta, phi) 
            A_m   = 30;
            SLA_v = 30;
            vert_deg_lobe = 65;
            horz_deg_lobe = 65;
            A_V   = -min(12*((theta-90)/vert_deg_lobe).^2,SLA_v);
            A_H   = -min(12*(phi/horz_deg_lobe).^2,A_m);

            A_HV  = -min(-(A_V+A_H),A_m);
            
            antenna_gain = obj.max_antenna_gain + A_HV; 
%             antenna_gain = zeros(size(antenna_gain));
        end
        
        function minmaxgain = min_max_gain(obj)
            minmaxgain(1) = obj.gain(180,180,0);
            minmaxgain(2) = obj.max_antenna_gain;
        end
        
       %Antenna polarization based on TR 36.873 Model 2 for polarization
        function [field_pattern_theta, field_pattern_phi] = polarization(obj, theta, phi, slant_angle)
            field_pattern_theta = sqrt(10.^(obj.gain(theta,phi)./10)).*cosd(slant_angle);
            field_pattern_phi   = sqrt(10.^(obj.gain(theta,phi)./10)).*sind(slant_angle);
        end
        
        
        %Added from 'TR36873_3DAntenna' lines 43 to 74
        % Calculate the gain pattern of antenna array single column with
        % M elements in vertical dimension
        function array_gain = calculate_array_single_column_field_pattern(obj, config, theta, phi, attached_site_pos)
%            PHI =0;
%            THETA = 0:pi/100:pi; %(-90:1.8:90)*pi/180; %
%        
%            wavelength = 299792458./2e9;
%            sin_cos_p = shiftdim(sin(THETA).*cos(PHI),-1);
%            sin_sin_p = shiftdim(sin(THETA).*sin(PHI),-1);
%            cos_p = shiftdim(cos(THETA),-1);
% 
%            spherical_unit_vector_r_p = [sin_cos_p; sin_sin_p; cos_p];
%            size_vec_= size(spherical_unit_vector_r_p);
%            spherical_unit_vector_p_ = reshape(spherical_unit_vector_r_p,3,size_vec_(2)*size_vec_(3));
%            
%            element_position_ = zeros(config.nr_of_antenna_elements_in_each_column,3);
%            eNodeB_location_ = zeros(config.nr_of_antenna_elements_in_each_column,3);
%            antenna_element_weight_ = zeros(config.nr_of_antenna_elements_in_each_column,1);
% % %            antenna_element_weight_ = zeros(size(THETA,2),LTE_config.nr_of_antenna_elements_in_each_column);
% %            mm= [1:config.nr_of_antenna_elements_in_each_column];
% %             antenna_element_weight_(mm) = 1./sqrt(config.nr_of_antenna_elements_in_each_column).*exp(-1i*2*pi./wavelength*(mm-1)*config.antenna_element_vertical_spacing*cosd(90));
%            % relative_element_position and weight vector
%            for m=1:config.nr_of_antenna_elements_in_each_column
% % %                antenna_element_weight_(:,m) = 1/sqrt(LTE_config.nr_of_antenna_elements_in_each_column)*...
% % %                    exp(1i*2*pi/wavelength*(m-1)*LTE_config.antenna_element_vertical_spacing*(cos(THETA)-cosd(LTE_config.electrical_downtilt)));
%                element_position_(m,:) = [0,0,config.antenna_element_vertical_spacing*(m-1)];
%                antenna_element_weight_(m) = 1/sqrt(config.nr_of_antenna_elements_in_each_column)*exp(-1i*2*pi/wavelength*(m-1)*config.antenna_element_vertical_spacing*cosd(config.electrical_downtilt));
%                eNodeB_location_(m,:) = [attached_site_pos(1),attached_site_pos(2),config.tx_height]+element_position_(m,:);
%                
%                calculate_r_tx_and_location_vector_(:,:,m) = reshape(spherical_unit_vector_p_.'*eNodeB_location_(m,:).',size_vec_(2),size_vec_(3));
%                array_factor_elementwise_terms_(:,:,m) =  antenna_element_weight_(m)*exp(1i*2*pi/wavelength.*calculate_r_tx_and_location_vector_(:,:,m));
%            end  
%            
% % %            sumArr=sum(antenna_element_weight_,2);
%            array_factor_one_column_ = sum(array_factor_elementwise_terms_,3);
%            figure(3); polar(-THETA+pi/2,abs(array_factor_one_column_).*sqrt(10.^(obj.gain(180/pi*THETA,180/pi*PHI)./10)),'b');title('Antenna array radiation pattern in elevation');
%            figure(4); plot(THETA*180/pi, abs(array_factor_one_column_));

           %% Antenna pattern in azimuth  
           PHI =-pi:pi/100:pi;
           THETA = pi/2;
           wavelength = 299792458./2e9;

           sin_cos_p_ = shiftdim(sin(THETA).*cos(PHI),-1);
           sin_sin_p_ = shiftdim(sin(THETA).*sin(PHI),-1);
           cos_p_ = shiftdim(cos(THETA),-1);
           cos_p_=ones(1,1,201)*cos_p_;

           spherical_unit_vector_r_p_ = [sin_cos_p_; sin_sin_p_; cos_p_];
           size_vec_p= size(spherical_unit_vector_r_p_);
           spherical_unit_vector_p_p = reshape(spherical_unit_vector_r_p_,3,size_vec_p(2)*size_vec_p(3));
           
           element_position_ = zeros(config.nr_of_antenna_elements_in_each_column,3);
           eNodeB_location_ = zeros(config.nr_of_antenna_elements_in_each_column,3);
           antenna_element_weight_ = zeros(config.nr_of_antenna_elements_in_each_column,1);
           % relative_element_position and weight vector
           for m=1:config.nr_of_antenna_elements_in_each_column
               element_position_(m,:) = [0,0, config.antenna_element_vertical_spacing*(m - 1)];
               antenna_element_weight_(m) = 1/sqrt(config.nr_of_antenna_elements_in_each_column)*exp(-1i*2*pi/wavelength*(m-1)*config.antenna_element_vertical_spacing*cosd(config.electrical_downtilt));
               eNodeB_location_(m,:) = [attached_site_pos(1),attached_site_pos(2),config.tx_height]+element_position_(m,:);
               
               calculate_r_tx_and_location_vector_p(:,:,m) = reshape(spherical_unit_vector_p_p.'*eNodeB_location_(m,:).',size_vec_p(2),size_vec_p(3));
               array_factor_elementwise_terms_p(:,:,m) =  antenna_element_weight_(m)*exp(1i*2*pi/wavelength.*calculate_r_tx_and_location_vector_p(:,:,m));
           end

           array_factor_one_column_p = sum(array_factor_elementwise_terms_p,3);
       
           % radiation pattern in elevation
           % figure(); polar(PHI,abs(array_factor_one_column_p).*sqrt(10.^(obj.gain(180/pi*THETA,180/pi*PHI)./10)),'m');title('Antenna array radiation pattern in azimuth');
            
            wavelength = 299792458./config.frequency;
            %spherical unit vector
            sin_cos = shiftdim(sind(theta).*cosd(phi),-1);
            sin_sin = shiftdim(sind(theta).*sind(phi),-1);
            cos_ = shiftdim(cosd(theta),-1);
            
            spherical_unit_vector_r = [sin_cos; sin_sin; cos_];
            size_vec= size(spherical_unit_vector_r);
            spherical_unit_vector_ = reshape(spherical_unit_vector_r,3,size_vec(2)*size_vec(3));
            
            element_position = zeros(config.nr_of_antenna_elements_in_each_column,3);
            eNodeB_location = zeros(config.nr_of_antenna_elements_in_each_column,3);
            antenna_element_weight = zeros(config.nr_of_antenna_elements_in_each_column,1);
            
            % relative_element_position and weight vector
            for m=1:config.nr_of_antenna_elements_in_each_column
                element_position(m,:) = [0, 0, config.antenna_element_vertical_spacing*(m - 1)];
                antenna_element_weight(m) = 1/sqrt(config.nr_of_antenna_elements_in_each_column)*...
                                             exp(-1i*2*pi/wavelength*(m-1)*config.antenna_element_vertical_spacing*cosd(config.electrical_downtilt));
                eNodeB_location(m,:) = [attached_site_pos(1),attached_site_pos(2),config.tx_height]+element_position(m,:);
                calculate_r_tx_and_location_vector(:,:,m) = reshape(spherical_unit_vector_.'*eNodeB_location(m,:).',size_vec(2),size_vec(3));
                array_factor_elementwise_terms(:,:,m) =  antenna_element_weight(m)*exp(1i*2*pi/wavelength.*calculate_r_tx_and_location_vector(:,:,m));
            end
            
            array_factor_one_column = sum(array_factor_elementwise_terms,3);
            %Antenna array gain
            array_gain = obj.gain(theta, phi) +20*log10(abs(array_factor_one_column));
            
        end
       
        %Changed from 'TR36873_3DAntenna' lines 79 to 85
        % Returns a horizontal and vertical antenna gain plot for a
        % specific tilt value. Returns antenna gain in dBi
       function [hor_degrees hor_gain ver_degrees ver_gain max_gain] = gain_patterns(obj)
           max_gain    =  obj.max_antenna_gain;
           hor_degrees = -179:180;
           hor_gain    = obj.gain(0, hor_degrees);
           ver_degrees = -89:90;
           ver_gain    = obj.gain(ver_degrees, 0);
       end
    end
    
end