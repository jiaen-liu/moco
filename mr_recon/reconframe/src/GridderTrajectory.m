classdef GridderTrajectory
    % Functions to perform gridding
    
    methods (Static, Hidden)
        function [k, grid_ovs, rot_angles] = calculate_trajectory( data, kx_ovs, kx_range, GridderPars, nus_enc_nrs, sense_factor, echoes, samples )
            rot_angles = [];
            k = [];
            grid_ovs = [];
            
            if ~isempty (data)
                
                % Get number of points to be omitted at the beginning of
                % the spiral trajectory
                origin_shift = 0;
                if isempty( GridderPars.Trajectory)
                    if strcmpi( GridderPars.Preset, 'spiral' )
                        if ~isempty( GridderPars.SpiralLeadingSamples )
                            origin_shift = GridderPars.SpiralLeadingSamples;
                        else
                            slice = ceil(size(data,3)/2);
                            [mi,origin_shift] = min(range(angle(data(1:20,:,slice,:,:)),2));
                            origin_shift = round(median(squeeze(origin_shift)));
                            GridderPars.SpiralLeadingSamples = origin_shift;
                        end
                    end
                end
                
                %---------------------------------------------------------------------
                no_samples      = size(data,1) - origin_shift;
                no_profiles     = size(data,2);
                no_interleaves  = size(data,2);
                no_slices       = size(data,3);
                %---------------------------------------------------------------------
                switch GridderPars.Preset
                    case {'Radial', 'radial', 'RADIAL'}
                        
                        if isempty( GridderPars.GridOvsFactor )
                            grid_ovs = kx_ovs / 2;
                        else
                            grid_ovs = [];
                        end
                        
                        if strcmpi( GridderPars.KooshBall, 'yes' )
                            profiles = linspace(-(no_samples-1)/2,(no_samples-1)/2,no_samples);
                            
                            k = zeros(no_samples,no_profiles,no_slices,3, 'single');
                            
                            %     Darf man py_max / pz_max / radial_3d_koosh_interleaves / radial_nr_angles so berechnen? DEBUG?
                            
                            if ~isempty( GridderPars.RadialAngles ) && (size(GridderPars.RadialAngles, 1) == 2) && (size(GridderPars.RadialAngles, 2) > 2)
                                GridderPars.RadialAngles = GridderPars.RadialAngles';
                            end
                            
                            for pz_number = 1:no_slices
                                for py_number = 1:no_profiles
                                    interleaf = pz_number-1;
                                    
                                    if isempty( GridderPars.RadialAngles )
                                        
                                        z = (py_number-1 + 0.5 - no_profiles) / no_profiles;
                                        % factor of 2. necessary since we cover only half a sphere */
                                        phi = (sqrt( 2. * no_profiles * pi / no_slices ) * asin( z )) + (interleaf * 2.0 * pi / no_slices);
                                        % sign alternation of the radial readout gradient      */
                                        if mod( py_number-1,2 ) == 1  % odd
                                            z = -z;
                                            phi = phi + pi;
                                        end
                                        
                                    else
                                        
                                        if size(GridderPars.RadialAngles, 2) ~= 2
                                            error( 'Error in Gridder: Please specify 2 angles (theta, phi) for the kooshball trajectory. Parameter.Gridder.RadialAngles has to be a vector of size %d x %d', no_profiles*no_slices,2);
                                        end
                                        if length(GridderPars.RadialAngles) < no_profiles*no_slices
                                            error( 'Error in Gridder: the number of user defined rotation angles is smaller than the number of measured profiles. Parameter.Gridder.RadialAngles has to be a vector of size %d x %d', no_profiles*no_slices,2);
                                        end
                                        theta = GridderPars.RadialAngles( (pz_number-1)*no_profiles + py_number, 1 );
                                        phi = GridderPars.RadialAngles( (pz_number-1)*no_profiles + py_number, 2 );
                                        z = cos(theta);
                                        
                                    end
                                    
                                    sin_theta = sqrt( 1.0 - z*z );
                                    rsinphi = profiles .* sin(phi);
                                    rcosphi = profiles .* cos(phi);
                                    
                                    k(:,py_number,pz_number,1) =  rsinphi .* sin_theta;
                                    k(:,py_number,pz_number,2) =  rcosphi .* sin_theta;
                                    k(:,py_number,pz_number,3) =  profiles .* z;
                                end
                            end
                        elseif strcmpi( GridderPars.UTE, 'yes' ) && ismember(0, echoes)
                            % TODO: define the UTE trajectory here
                            error('UTE is not yet implemented in MReconImpl. You will need to define the trajectory yourself');
                            
                        else
                            k = zeros(no_samples,3,no_profiles,'single');
                            rot_angles = zeros( no_profiles, 1);
                            
                            %---------------------------------------------------------------------
                            % Calculate radial trajectory
                            if no_samples/2~=floor(no_samples/2)
                                k0 = single([zeros(1,no_samples);linspace(-floor(no_samples/2),floor(no_samples/2),no_samples)]);
                            else
                                k0 = single([zeros(1,no_samples);linspace(-floor(no_samples/2),ceil(no_samples/2-1),no_samples)]);
                            end
                            
                            for i = 0:no_profiles-1
                                if isempty( GridderPars.RadialAngles )
                                    if strcmpi( GridderPars.UTE, 'yes' )
                                        rot_angle = -i*2*pi/no_profiles;
                                    elseif strcmpi( GridderPars.GoldenAngle, 'yes' )
                                        rot_angle = -i*GridderPars.GoldenAngleValue/180*pi;
                                    elseif strcmpi( GridderPars.PseudoGoldenAngle, 'yes' )
                                        rot_angle = -i*2*pi/no_profiles;
                                    else
                                        rot_angle = -i*pi/no_profiles;
                                    end
                                    
                                    if strcmpi( GridderPars.AlternatingRadial, 'yes') && mod( i, 2)
                                        rot_angle = rot_angle + pi;
                                    end
                                else
                                    if length(GridderPars.RadialAngles) < no_profiles
                                        error( 'Error in Gridder: the number of user defined rotation angles is smaller than the number of measured profiles. Parameter.Gridder.RadialAngles has to be a vector of size %d x %d', no_profiles,1);
                                    end
                                    rot_angle = GridderPars.RadialAngles(i+1);
                                end
                                rot_angles(i+1) = rot_angle;
                                R = [cos(rot_angle),-sin(rot_angle);sin(rot_angle),cos(rot_angle)];
                                k(:,1:2,i+1) = (R*k0)';
                            end
                            k = permute(k,[1,3,4,2]);
                            clear k0 R;
                        end
                    case {'Spiral', 'spiral', 'SPIRAL'}
                        if isempty( GridderPars.GridOvsFactor )
                            grid_ovs       = kx_ovs / 1;
                        else
                            grid_ovs = [];
                        end
                        
                        %---------------------------------------------------------------------
                        % Spiral fine-tuning (scanner specific)
                        channel_delay  = 0;
                        phase_offset   = 0;
                        lambda         = 3;
                        
                        %---------------------------------------------------------------------
                        % Definition of k-space trajectory is found in mpuspiral__g.c in
                        % Philips pulse programming environment
                        acq_samples    = round(length(kx_range(1):kx_range(2) ) / kx_ovs);
                        no_turns       = acq_samples /(2*no_interleaves);                       % ok
                        lambda         = min(no_turns,lambda);                                  % /Ta
                        phi_top        = pi*acq_samples*sqrt(1+lambda)/(no_interleaves*no_samples);
                        A              = no_interleaves/(2*pi);
                        k              = zeros(no_samples,no_interleaves, 3,'single');
                        
                        %--------------------------------------------------
                        %-------------------
                        % Calculate spiral trajectory   
                        if ~isempty( GridderPars.Trajectory)
                            for interleave=0:no_interleaves-1
                                phi = -2*pi*interleave/no_interleaves+phase_offset;
                                R = [cos(phi), -sin(phi), 0; sin(phi), cos(phi), 0; 0, 0, 1];
                                cur_k = R*GridderPars.Trajectory;
                                k(:,interleave+1,:) = permute(cur_k, [2,1]);
                            end
                            k = k .* acq_samples;
                        else
                            for interleave=0:no_interleaves-1
                                phi_0 = 2*pi*interleave/no_interleaves+phase_offset;
                                for sample=fix(channel_delay):no_samples-1
                                    samp = channel_delay - fix(channel_delay) + sample;
                                    phi_t = phi_top*samp/sqrt(1+lambda*samp/no_samples);
                                    k(sample+1,interleave+1,1) = A*phi_t*cos(phi_t-phi_0);
                                    k(sample+1,interleave+1,2) = A*phi_t*sin(phi_t-phi_0);
                                end
                            end
                        end
                        k = reshape(k, size(k,1), size(k,2), 1, size(k,3) );
                    case {'Epi', 'epi', 'EPI'}
                        k = [];
                        if ~isempty( nus_enc_nrs )
                            if length( nus_enc_nrs ) ~= size(data,1)
                                error( 'Error in GridData: The length of the EPI NUS samples and nr of dsamples in Data is different --> cannot grid the EPI data' );
                            end
                            k = zeros( size(data,1), 1, 1, 3);
                            kx = nus_enc_nrs;
                            k(:,:,:,1) = kx;
                            k(:,:,:,2) = 0;
                            k(:,:,:,3) = 0;
                            grid_ovs = 1;
                        end
                    case {'Cartesian', 'cartesian', 'CARTESIAN'}
                        r_res = [size(data, 1),size(data, 2), size(data, 3)];
                        res = round(sense_factor(1,:).*r_res);
                        
                        for i = 1:3
                            pos = 0:sense_factor(1,i):ceil(res(i)/2-1);
                            neg = -sense_factor(1,i):-sense_factor(1,i):-floor(res(i)/2);
                            kxyz{i} = [neg(end:-1:1), pos];
                        end
                        [kx, ky, kz] = ndgrid( kxyz{1}, kxyz{2}, kxyz{3} );
                        k = zeros( length(kx(:)), 3, 'single');
                        k(:,1) = kx(:);
                        k(:,2) = ky(:);
                        k(:,3) = kz(:);
                        k = reshape(k, r_res(1),r_res(2),r_res(3), 3);
                        grid_ovs = kx_ovs / 1;
                    otherwise
                        error( 'Error in Gridder: Gridder Preset Unknown');
                end
            end
        end
        function weight = calculate_weights( k, GridderPars )
            weight = [];
            if ~isempty(k)
                switch GridderPars.Preset
                    case {'Radial', 'radial', 'RADIAL'}
                        if strcmpi( GridderPars.KooshBall, 'yes' )
                            weight = k(:,:,:,1).^2 + k(:,:,:,2).^2 + k(:,:,:,3).^2;
                            weight = weight./(size(k,1)/2)^2;
                        else
                            %---------------------------------------------------------------------
                            % Calculate simple k-space weights
                            gk = k(2:end,:,:,:)-k(1:end-1,:,:,:);
                            gk(end+1,:,:,:) = gk(end,:,:,:);
                            weight = abs(k(:,:,:,1).*gk(:,:,:,1)+k(:,:,:,2).*gk(:,:,:,2));
                            weight(end,:,:,:) = weight(end-1,:,:,:);
                            weight = weight./( max([size(k,1), size(k,2)]) /2 );
                            clear gk;
                        end
                    case {'Spiral', 'spiral', 'SPIRAL'}
                        %---------------------------------------------------------------------
                        % Calculate simple k-space weights
                        gk = abs( abs( k(2:end,:,:,1) + 1i*k(2:end,:,:,2) ) - abs( k(1:end-1,:,:,1) + 1i*k(1:end-1,:,:,2) ) );
                        gk(end+1,:,:,:) = gk(end,:,:,:);
                        weight = abs( k(:,:,:,1) + 1i*k(:,:,:,2) ) .* gk;
                        weight(end,:,:,:) = weight(end-1,:,:,:);
                        weight( 1,:,:,:) = 0;
                        clear gk;
                    case {'Epi', 'epi', 'EPI'}
                        weight = ones( size(k,1), 1, 1 );
                    case {'Cartesian', 'cartesian', 'CARTESIAN'}
                        weight = ones( size(k,1), size(k,2), size(k,3) );
                    otherwise
                        weight = ones( size(k,1), size(k,2), size(k,3) );
                end
            end
        end
    end
end