classdef MRecon < handle
    % MRecon is an extensive object-oriented Matlab library, implementing many
    % common image and spectrum reconstruction task. The object-oriented design
    % makes it easy to expand or alter the existing functionality. Code examples
    % of complete reconstructions are provided. They can serve as a starting point
    % for custom development.
    
    properties (Dependent)
        Parameter;  % Scan and reconstruction parameter
        Data;       % The current data
        DataProperties;  % Holds information about the Data dimensions and labels
    end
    properties (Dependent, Hidden)
         DataClass;
    end
    properties (Hidden)
        MRImpl;       
    end
    properties (Hidden, Constant)
        dim = struct('kx', 		1, ...
            'ky',               2, ...
            'kz',               3, ...
            'coil',             4, ...
            'dyn',              5, ...
            'hp',               6, ...
            'echo',             7, ...
            'loc',              8, ...
            'mix',              9, ...
            'extr1',            10, ...
            'extr2',            11, ...
            'meas',             12);
    end
    methods
        % ---------------------------------------------------------------%
        % Constructor
        % ---------------------------------------------------------------%
        function MR = MRecon( varargin )
            p = inputParser;   
            filenameValidation = @(x) ischar(x) || isempty(x);
            retroHoleValidation = @(x) ischar(x) && (strcmpi(x, 'no') || strcmpi(x, 'nearest') || strcmpi(x, 'average') || strcmpi(x, 'linear') || strcmpi(x, 'cubic'));
            addOptional(p,'Filename', [], filenameValidation);
            addOptional(p,'Datafile', [], filenameValidation);
            addParameter(p,'Sinfile', [], filenameValidation);  
            addParameter(p,'Dicomfile', [], filenameValidation);  
            addParameter(p,'Rcfile', [], filenameValidation);   
            addParameter(p,'FileDialogTitle', 'Pick a file', @ischar);    
            addParameter(p,'RetroHoleInterpolation', [], retroHoleValidation);    
            parse(p, varargin{:});  
            
            sinfile = p.Results.Sinfile;
            if isempty(p.Results.Sinfile) && ~isempty(p.Results.Rcfile)
                sinfile = p.Results.Rcfile;
            end
            
            MRecon.AddPath;
            MR.MRImpl = MReconImpl(p.Results.Filename, p.Results.Datafile, sinfile, p.Results.Dicomfile, p.Results.FileDialogTitle, p.Results.RetroHoleInterpolation);                                    
        end
        
        % ---------------------------------------------------------------%
        % Perform Reconstruction
        % ---------------------------------------------------------------%
        function Perform( MR )
            % Perform: Performs a complete reconstruction of the selected MR data. Please see the
            % source code below for meore details.
            %
            % Syntax:     r.Perform;
            
            switch MR.Parameter.DataFormat
                case{'ExportedRaw', 'Raw', 'Bruker' }
                    
                    %Reconstruct only standard (imaging) data
                    MR.Parameter.Parameter2Read.typ = 1;
                    MR.Parameter.Parameter2Read.Update;
                    
                    % Check if enough memory is available to reconstruct
                    % the whole file at once. Otherwise reconstruct the
                    % data in chunks
                    if strcmpi(MR.Parameter.Recon.AutoChunkHandling, 'yes')
                        [MemoryNeeded, MemoryAvailable] = MR.GetMemoryInformation;
                        if MemoryNeeded > MemoryAvailable
                            if strcmpi( MR.Parameter.Recon.ImmediateAveraging, 'yes' ) || strcmpi( MR.Parameter.Recon.Average, 'yes' )
                                MR.Parameter.Chunk.Def = {'kx', 'ky', 'kz', 'chan', 'aver'};
                            else
                                MR.Parameter.Chunk.Def = {'kx', 'ky', 'kz', 'chan'};
                            end
                        end
                    end
                    
                    % Switch off for performance reasons (after recon it is
                    % switched back to original state again)
                    AutoUpdateStatus = MR.Parameter.Recon.AutoUpdateInfoPars;
                    MR.Parameter.Recon.AutoUpdateInfoPars = 'no';
                    
                    % Define new counter
                    counter = Counter( 'Performing Recon --> Chunk %d/%d\n');
                    
                    % Loop over all chunks
                    for cur_loop = 1:MR.Parameter.Chunk.NrLoops
                        
                        % Update Counter
                        if strcmpi( MR.Parameter.Recon.StatusMessage, 'yes')
                            counter.Update( {cur_loop ,MR.Parameter.Chunk.NrLoops} );
                        end
                        
                        % Set the chunk-loop which automatically determines the
                        % image parameter to read for the current chunk
                        MR.Parameter.Chunk.CurLoop = cur_loop;
                        
                        % --------------------------------------------------------
                        % Perform the Reconstruction for the Current Chunk (Start)
                        % --------------------------------------------------------
                        
                        % spectro begin ----------------------------
                        if MR.isSpectro
                            MR.ReadData;  
                            MR.NonLinearityCorrection;
                            MR.RandomPhaseCorrection;
                            MR.PDACorrection;
                            MR.DcOffsetCorrection;
                            MR.SortData;
                            MR.Average;
                            MR.RemoveOversampling;
                            MR.RingingFilter;
                            MR.ZeroFill;
                            MR.SENSEUnfold;
                            MR.EddyCurrentCorrection;
                            MR.CombineCoils;
                            MR.GeometryCorrection;
                            MR.K2I;
                            %MR.RotateImage;
                        else
                            % spectro end ------------------------------
                            MR.ReadData;
                            MR.NonLinearityCorrection;
                            MR.RandomPhaseCorrection;
                            MR.RemoveOversampling;
                            MR.PDACorrection;
                            MR.DcOffsetCorrection;
                            MR.MeasPhaseCorrection;
                            MR.SortData;
                            MR.GridData;
                            MR.RingingFilter;
                            MR.ZeroFill;
                            MR.K2IM;
                            MR.EPIPhaseCorrection;
                            MR.K2IP;
                            MR.GridderNormalization;
                            MR.SENSEUnfold;
                            MR.PartialFourier;
                            MR.ConcomitantFieldCorrection;
                            MR.DivideFlowSegments;
                            MR.CombineCoils;
                            MR.Average;
                            MR.GeometryCorrection;
                            MR.RemoveOversampling;                            
                            MR.ReconTKE;
                            MR.FlowPhaseCorrection;
                            MR.ZeroFill;
                            MR.RotateImage;
                        end
                        
                        % The current chunk is now reconstructed. If the data is
                        % reconstructed in more than one chunk write the result to
                        % a temporary file on the disk.
                        if MR.Parameter.Chunk.NrLoops > 1
                            [exported_datafile, exported_listfile] = MR.WriteExportedRaw( [MR.Parameter.Filename.Data, '_temp.data'], MR.Parameter.Parameter2Read );
                        end
                        
                        % --------------------------------------------------------
                        % Perform the Reconstruction for the Current Chunk (End)
                        % --------------------------------------------------------
                    end
                    
                    % If data has been written to a temporary file read it
                    % again
                    if MR.Parameter.Chunk.NrLoops > 1
                        r_temp = MRecon(exported_datafile);
                        r_temp.ReadData;
                        r_temp.Parameter.Recon.ImmediateAveraging = 'no';
                        r_temp.SortData;
                        MR.Data = r_temp.Data;                       
                        delete(exported_datafile);
                        delete(exported_listfile);
                        clear r_temp;
                    end
                    if strcmpi( MR.Parameter.Recon.StatusMessage, 'yes')
                        fprintf('\n');
                    end
                    MR.Parameter.Recon.AutoUpdateInfoPars = AutoUpdateStatus;
                    MR.Parameter.Reset;
                case 'ExportedCpx'
                    MR.ReadData;
                    MR.SortData;
                    MR.CombineCoils;
                case 'Cpx'
                    MR.ReadData;
                    MR.CombineCoils;
                case 'Rec'
                    MR.ReadData;
                    MR.RescaleREC;
                    MR.CreateComplexREC;
                otherwise
                    error( 'Error in Perform: Unknown data format' );
            end
        end
        
        % ---------------------------------------------------------------%
        % Data Reading
        % ---------------------------------------------------------------%
        function ReadData( MR, typ )
            % ReadData: Reads all supported data formats into Matlab.
            %
            % Syntax:     r.ReadData;
            %
            % Parameters used:
            %
            %           All formats:
            %           - <a href="matlab:helpwin('MRparameterDoc.DataFormat')">Parameter.DataFormat</a>:
            %             Used to check which data format is passed to the reader
            %           - <a href="matlab:helpwin('Parameter2ReadParsDoc')">Parameter.Parameter2Read</a>:
            %             Specifies what data is read from the file.
            %           - <a href="matlab:helpwin('MRparameterDoc.Filename')">Parameter.Filename.ijk</a>:
            %             The filename of the data and parameter file to be read.
            %           - <a href="matlab:helpwin('MRparameterDoc.Labels')">Parameter.Labels</a>:
            %             The data labels which hold the image attributes for every profile.
            %
            %           Raw and ExportedRaw data
            %           - <a href="matlab:helpwin('EncodingParsDoc.KxRange')">Parameter.Encoding.KxRange</a>:
            %             The sampled k-space range in readout direction used to place the data
            %             correctly in k-space (e.g. in partial echo acquisitions). Please note that
            %             for data which has to be gridded (e.g. radial, spiral, EPI data), this
            %             parameter specifies the k-space range after gridding.
            %           - <a href="matlab:helpwin('ReconParsDoc.ArrayCompression')">Parameter.Recon.ArrayCompression</a>:
            %             Enables the array compression, according to: Buehrer, M., Pruessmann, K. P.
            %             , Boesiger, P. and Kozerke, S. (2007), Array compression for MRI with large
            %               coil arrays. Magnetic Resonance in Medicine, 57: 1131–1139.
            %             The individual profiles are compressed to the userd defined number of
            %             virtual coils immediately after reading.
            %           - <a href="matlab:helpwin('ReconParsDoc.ACNrVirtualChannels')">Parameter.Recon.ACNrVirtualChannels</a>:
            %             The number of virtual coils used in array compression.
            %           - <a href="matlab:helpwin('ReconParsDoc.ACMatrix')">Parameter.Recon.ACMatrix</a>:
            %             The compression matrix used in array compression.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.isread')">Parameter.ReconFlags.isread</a>
            %
            % Location:   Beginning of recon
            %
            % Formats:    Raw | ExportedRaw | Cpx | Rec | Bruker
            %
            % Description/Algorithm: The ReadData function reads the MR data from the selected file
            %             and stores it into Matlab. The function supports selective reading,
            %             meaning that the user can specify what data should be read via
            %             r.Parameter.Parameter2Read. After calling the ReadData function, the data
            %             is stored in the order it was saved in the data file. This means that for
            %             raw data we obtain a matrix containing the individual k-space profiles in
            %             the order they have been measured. For Rec and Cpx data however we obtain
            %             (partly) reconstructed images in image-space.
            %             Generally the data is stored in a cell array where the different rows hold
            %             different data types (e.g. imaging data, correction data, noise data etc)
            %             and the columns different data sizes (e.g. kt-undersampled and training
            %             data). The number of rows in the cell array is predefined while the number
            %             of columns varies depending on the data in the file. Generally the cell
            %             array has the following form:
            %
            %                                               Data size 1     Data size 2
            %
            %               1. Accepted imaging data        { STD },        { STD }
            %               2. Rejected imaging data        { REJ },        { REJ }
            %               3. Phase correction data        { PHX },        { PHX }
            %               4. Frequency correction data    { FRX },        { FRX }
            %               5. Noise data                   { NOI },        { NOI }
            %
            %             I a data type does not occur in the file the the selected cell element is
            %             left empty. If only one data type and size is read (cell array with one
            %             element) then a ordinary Matrix is returned instead of a cell.
            %
            %             For convenience the data type can directly be passed to ReadData as argument:
            %                          r.ReadData(1)
            %             is equivalent to:
            %                          r.Parameter.Parameter2Read.typ = 1;
            %                          r.ReadData;
            %
            % Examples:   Read only imaging data from the first coil and the third dynamic:
            %
            %               r = MRecon( 'rawfile.raw' );
            %               r.Parameter.Parameter2Read.typ  = 1;   % Imaging data
            %               r.Parameter.Parameter2Read.chan = 0;   % First coil
            %               r.Parameter.Parameter2Read.dyn  = 2;   % Third dynamic
            %               r.ReadData;
            %
            %             Read the noise samples from a raw file:
            %               r = MRecon( 'rawfile.raw' );           
            %               r.ReadData(5);  % 5 = Noise data
            
            if nargin == 1
                typ = NaN;
            end 
            if ~isnan(typ)
                 MR.Parameter.Parameter2Read.typ = typ;
            end
            MR.Parameter.ReconFlags.Init(MR.Parameter.DataFormat);
            MR.Parameter.InitWorkEncoding( MR.Parameter );
            MR.Parameter.Gridder.InitWorkingPars;
            try
                MR.Parameter.InitCurFOV( 'ReadData');
            catch
            end
            MR.Parameter.ResetImageInformation;
            
            switch MR.Parameter.DataFormat
                case 'Rec'
                    MR.Parameter.Recon.AutoUpdateInfoPars = 'no';
                    [MR.Data, MR.Parameter.LabelLookupTable, MR.Parameter.ImageInformation] = Reader.readrec( MR.Parameter.Filename.Data, MR.Parameter.Parameter2Read, MR.Parameter.Labels );
                    MR.Parameter.Recon.AutoUpdateInfoPars = 'yes';
                case 'Cpx'
                    MR.Data = Reader.read_cpx(MR.Parameter.Filename.Data, 0, 0, 0, MR.Parameter.Parameter2Read, []);
                case {'ExportedRaw', 'Raw', 'ExportedCpx', 'Bruker'}                               
                    MR.Parameter.Scan.ijk = MR.Parameter.Scan.MPS;
                    % Check if it is a radial or spiral scan
                    radial_spiral =  any( strcmpi( MR.Parameter.Scan.AcqMode, {'radial', 'spiral'})) || strcmpi( MR.Parameter.Gridder.Preset, 'radial' );
                    [MR.Data, MR.Parameter.LabelLookupTable] = ...
                        cellfun( @(x1, x2, x3, x4, x5, x6)Reader.read_raw( MR.Parameter, ...                        
                        MR.Parameter.Parameter2Read, x1, x2, x3, x4, ...
                        radial_spiral, strcmpi( MR.Parameter.Recon.ArrayCompression, 'yes' ), ...
                        MR.Parameter.Recon.ACNrVirtualChannels, MR.Parameter.Recon.ACMatrix ), ...
                        MR.Parameter.Encoding.WorkEncoding.Typ, ...
                        MR.Parameter.Encoding.WorkEncoding.Mix, ...
                        MR.Parameter.Encoding.WorkEncoding.Echo, ...
                        MR.Parameter.Encoding.WorkEncoding.KxRange, ...
                        'UniformOutput', false );
            end
            
            MR.Data = Helper.UnconvertCell( MR.Data );
            MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
            MR.Parameter.ReconFlags.isread = 1;
        end
        function ReadChunk(MR, varargin)
            if ~any(strcmp(MR.Parameter.DataFormat, {'ExportedRaw', 'Raw', 'ExportedCpx', 'Bruker'}))
                error('This function only works on raw data');
            end
            
            p = inputParser;                        
            sizeValidator = @(x) length(x) == 1 && isnumeric(x);
            resetValidator = @(x) islogical(x) || (isnumeric(x) && all(x <= 1) && all(x>=0));
            addParameter(p,'reset', false, resetValidator);
            addParameter(p,'size', 10, sizeValidator);
            addParameter(p,'typ', 1, sizeValidator);
            parse(p,varargin{:});
                                                
            if p.Results.reset
                MR.MRImpl.ChunkInfo = [];
            end
            
            MR.Parameter.Parameter2Read.typ = p.Results.typ;
            
            if( isempty(MR.MRImpl.ChunkInfo))
                MR.MRImpl.ChunkInfo.cur_ind = 1;
                MR.MRImpl.ChunkInfo.chunk_ind = [];
            end
            if( ~isfield(MR.MRImpl.ChunkInfo, 'size'))
                MR.MRImpl.ChunkInfo.size = p.Results.size;
            end            
            
            MR.Parameter.Scan.ijk = MR.Parameter.Scan.MPS;
            % Check if it is a radial or spiral scan
            radial_spiral =  any( strcmpi( MR.Parameter.Scan.AcqMode, {'radial', 'spiral'})) || strcmpi( MR.Parameter.Gridder.Preset, 'radial' );
            [MR.Data, MR.Parameter.LabelLookupTable, MR.MRImpl.ChunkInfo] = ...
                cellfun( @(x1, x2, x3, x4, x5, x6)Reader.read_raw( MR.Parameter, ...                
                MR.Parameter.Parameter2Read, x1, x2, x3, x4, ...
                radial_spiral, strcmpi( MR.Parameter.Recon.ArrayCompression, 'yes' ), ...
                MR.Parameter.Recon.ACNrVirtualChannels, MR.Parameter.Recon.ACMatrix, MR.MRImpl.ChunkInfo ), ...
                MR.Parameter.Encoding.WorkEncoding.Typ, ...
                MR.Parameter.Encoding.WorkEncoding.Mix, ...
                MR.Parameter.Encoding.WorkEncoding.Echo, ...
                MR.Parameter.Encoding.WorkEncoding.KxRange, ...
                'UniformOutput', false );
            
            MR.Data = Helper.UnconvertCell( MR.Data );
            MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
            MR.MRImpl.ChunkInfo = Helper.UnconvertCell( MR.MRImpl.ChunkInfo );
            MR.Parameter.ReconFlags.isread = 1;
        end
        
        % ---------------------------------------------------------------%
        % Data Writing
        % ---------------------------------------------------------------%
        function [datafile, listfile] = WriteExportedRaw( MR, Filename, Parameter2Write )
            % WriteExportedRaw: Exports the data to Philips exported raw format.
            %
            % Syntax:     [datafile, listfile] = r.WriteExportedRaw( Filename, Parameter2Write );
            %
            % Parameters used:
            %
            % Location:   k-space | image-space.
            %
            % Formats:    Raw | ExportedRaw
            %
            % Description/Algorithm: Exports the images in the data array to Philips exported format
            %             (data/list file pair). The Parameter2Write which have to be given as input
            %             is usually the current Parameter2Read struct (r.Parameter.Parameter2Read).
            
            if ~MR.Parameter.ReconFlags.issorted
                error( 'Please sort the data first' );
            end
            [datafile, listfile] = Writer.write_eported_raw( MR, Filename, MR.Data, Parameter2Write );
        end
        function WriteRec( MR, Filename )
            % WriteRec: Exports the data to Philips rec format.
            %
            % Syntax:     r.WriteRec( filename );
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.ExportRECImgTypes')">Parameter.Recon.ExportRECImgTypes</a>: {'M', 'P', 'R', 'I'}
            %             Specifies the images which should be exported (Magnitude, Phase, Real
            %             part, Imagenary part).
            %           - <a href="matlab:helpwin('InfoParsDoc.RescaleIntercept')">Parameter.ImageInformation.RescaleIntercept</a>, <a href="matlab:helpwin('InfoParsDoc.RescaleSlope')">Parameter.ImageInformation.RescaleSlope</a>:
            %             Defines the scaling of the images in the rec file. If these values are
            %             empty (default) then the scaling is calculated automatically in the
            %             WriteRec function. For details see below.
            %
            % Location:   image-space.
            %
            % Formats:    Raw | ExportedRaw | Cpx | Rec
            %
            % Description/Algorithm: Exports the images in the data array to Philips rec format.
            %             During the export the data is rescaled according to:
            %
            %             Magnitude data, Real part, Imagenary part:
            %                Value_in_Rec File = round( (Value_in_r.Data - RescaleIntercept) / RescaleSlope )
            %
            %             Phase data:
            %                Value_in_Rec_File = round( (1000 * Phase_in_r.Data - RescaleIntercept) / RescaleSlope )
            %
            % Notes:    - When the exported rec data is reimported into the scanner database the
            %             values are rescaled such that they correspond the the original values
            %             in r.Data before exporting.
                                                
            % define the data type of the recfile
            data_type = 'uint16';
            
            % define default filename
            if nargin == 1 || isempty( Filename )
                Filename = 'recon.rec';
            end
                                    
            if isempty(MR.Parameter.Scan.Multivenc) || strcmpi(MR.Parameter.Scan.Multivenc,'no')
                if (~isempty(MR.Parameter.Scan.Venc) && sum(MR.Parameter.Scan.Venc ~= 0)>1)
                    nr_files = size(MR.Data,10);
                    fid = cell(1, nr_files);
                    for ffile = 1:nr_files
                        if contains(Filename, '.rec' )
                            fid{ffile} = fopen( [strrep(Filename, '.rec', ''),sprintf('_%d.rec',ffile)],'w');
                        elseif  contains(Filename, '.REC' ) 
                            fid{ffile} = fopen( [strrep(Filename, '.REC', ''),sprintf('_%d.REC',ffile)],'w');
                        else
                            fid{ffile} = fopen( [Filename,sprintf('_%d.REC',ffile)],'w');
                        end
                    end
                else
                    fid{1} = fopen(Filename,'w');
                end
            else %multivenc
                venc_dir = sum(MR.Parameter.Scan.Venc>0,1);
                if sum(venc_dir)>0
                    if size(MR.Parameter.Scan.Venc, 1) > 1
                        nr_files = size(MR.Parameter.Scan.Venc, 1);
                    else
                        nr_files = sum(MR.Parameter.Scan.Venc > 0);
                    end
                    
                    if iscell(MR.Data)
                        nr_segments = size(MR.Data{1}, 10);                    
                    else
                        nr_segments = size(MR.Data, 10);
                    end
                    
                    if nr_files ~= nr_segments
                        error('The number of flow segments in r.Data does not match the Venc in r.Parameter.Scan');
                    end
                    
                    fid = cell(1, nr_files);
                    for nr=1:nr_files %export one magnitude/phase per flow encoded direction
                        if  contains(Filename, '.rec' ) 
                            fid{nr} = fopen( [strrep(Filename, '.rec', ''),sprintf('%d.rec',nr)],'w');
                        elseif  contains(Filename, '.REC' ) 
                            fid{nr} = fopen( [strrep(Filename, '.REC', ''),sprintf('%d.rec',nr)],'w');
                        else
                            fid{nr} = fopen( [Filename,sprintf('%d.rec',nr)],'w');
                        end
                    end
                    % if TKE maps available, then copy them to an extra dataset:
                    nr_images2 = 0;
                    if (strcmpi(MR.Parameter.Scan.Multivenc,'Yes')&& (size(MR.Data,2)==2))
                        nr_images2 = size( MR.Data{1, 2}(:,:,:), 3 );
                    end
                    
                    if nr_images2>0 %tke
                        if  contains(Filename, '.rec' ) 
                            fid{nr+1} = fopen( [strrep(Filename, '.rec', ''),'TKE.rec'],'w');
                        elseif  contains(Filename, '.REC' ) 
                            fid{nr+1} = fopen( [strrep(Filename, '.REC', ''),'TKE.rec'],'w');
                        else
                            fid{nr+1} = fopen( [Filename,'TKE.rec'],'w');
                        end
                    end
                else
                    if  contains(Filename, '.rec' ) 
                        fid{1} = fopen( [strrep(Filename, '.rec', ''),sprintf('%d.rec',1)],'w');
                    elseif  contains(Filename, '.REC' ) 
                        fid{1} = fopen( [strrep(Filename, '.REC', ''),sprintf('%d.rec',1)],'w');
                    else
                        fid{1} = fopen( [Filename,sprintf('%d.rec',1)],'w');
                    end
                end
            end
            
            rec_images = MR.GetRecImages;
            if length(rec_images) ~= length(fid)
                error('mismatch between number of recfiles and number of file contents');
            end
                        
            for i = 1:length(fid)                
                fwrite(fid{i},reshape( rec_images{i}, [], 1), data_type );
                fclose(fid{i});
            end                                                            
        end
        function WritePar( MR, Filename )
            % WritePar: Exports the parameter to a Philips par file.
            %
            % Syntax:     r.WritePar( filename );
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ScanParsDoc')">Parameter.Scan</a>:
            %             Holds the general parameters (uper part of the par file)
            %           - <a href="matlab:helpwin('InfoParsDoc')">Parameter.ImageInformation</a>:
            %             Holds the individual parameters for every image (lower part of the par
            %             file).
            %
            % Location:   image-space.
            %
            % Formats:    Raw | ExportedRaw | Cpx | Rec
            %
            % Description/Algorithm: Exports the parameters in to a Philips par file.                                   
            if nargin == 1 || isempty( Filename )
                Filename = 'recon.par';
            end                      
            
            if isempty(MR.Parameter.Scan.Multivenc) || (strcmpi(MR.Parameter.Scan.Multivenc,'no'))
                if (~isempty(MR.Parameter.Scan.Venc) && sum(MR.Parameter.Scan.Venc ~= 0)>1)
                    nr_files = size(MR.Data,10);
                    fid = cell(1, nr_files);
                    for ffile = 1:nr_files
                        if  contains(Filename, '.par' ) 
                            fid{ffile} = fopen( [strrep(Filename, '.par', ''),sprintf('_%d.par',ffile)],'w');
                        elseif  contains(Filename, '.PAR' ) 
                            fid{ffile} = fopen( [strrep(Filename, '.PAR', ''),sprintf('_%d.PAR',ffile)],'w');
                        else
                            fid{ffile} = fopen( [Filename,sprintf('_%d.PAR',ffile)],'w');
                        end
                    end
                else
                    fid{1} = fopen(Filename,'w');
                end
            else  %multivenc                                     
                venc_dir = sum(MR.Parameter.Scan.Venc>0,1);
                if sum(venc_dir)>0
                    if size(MR.Parameter.Scan.Venc, 1) > 1
                        nr_files = size(MR.Parameter.Scan.Venc, 1);
                    else
                        nr_files = sum(MR.Parameter.Scan.Venc > 0);
                    end
                    
                    if iscell(MR.Data)
                        nr_segments = size(MR.Data{1}, 10);                    
                    else
                        nr_segments = size(MR.Data, 10);
                    end
                        
                    if nr_files ~= nr_segments
                        error('The number of flow segments in r.Data does not match the Venc in r.Parameter.Scan');
                    end
                    
                    fid = cell(1, nr_files);
                    for nr=1:sum(nr_files) %export one magnitude/phase per flow encoded direction
                        if  contains(Filename, '.par' ) 
                            fid{nr} = fopen( [strrep(Filename, '.par', ''),sprintf('%d.par',nr)],'w');
                        elseif  contains(Filename, '.PAR' ) 
                            fid{nr} = fopen( [strrep(Filename, '.PAR', ''),sprintf('%d.par',nr)],'w');
                        else
                            fid{nr} = fopen( [Filename,sprintf('%d.par',nr)],'w');
                        end
                    end
                    % if TKE maps available, then copy them to an extra dataset:
                    nr_images2 = 0;
                    if (strcmpi(MR.Parameter.Scan.Multivenc,'Yes')&& (size(MR.Data,2)==2))
                        nr_images2 = size( MR.Data{1, 2}(:,:,:), 3 );
                    end
            
                    if nr_images2>0
                        if  contains(Filename, '.par' ) 
                            fid{nr+1} = fopen( [strrep(Filename, '.par', ''),'TKE.par'],'w');
                        elseif  contains(Filename, '.PAR' ) 
                            fid{nr+1} = fopen( [strrep(Filename, '.PAR', ''),'TKE.par'],'w');
                        else
                            fid{nr+1} = fopen( [Filename,'TKE.par'],'w');
                        end
                    end
                else
                    if  contains(Filename, '.par' ) 
                        fid{1} = fopen( [strrep(Filename, '.par', ''),sprintf('%d.par',1)],'w');
                    elseif  contains(Filename, '.PAR' ) 
                        fid{1} = fopen( [strrep(Filename, '.PAR', ''),sprintf('%d.par',1)],'w');
                    else
                        fid{1} = fopen( [Filename,sprintf('%d.par',1)],'w');
                    end
                end
            end
            
            ParString = MR.WritePar2String(Filename);
            if length(fid) ~= length(ParString)
                error('mismatch between number of files and number of file contents');
            end
            for i = 1:length(fid)
                fprintf( fid{i}, '%s', ParString{i});
                fclose(fid{i});
            end                                                         
        end
        function ParString = WritePar2String( MR, Filename )
            if isempty( MR.Parameter.ImageInformation )
                error( 'Parameter.ImageInformation struct is empty. Please set Parameter.Recon.AutoUpdateInfoPars to yes' );
            end
            
            if any( any( MR.Parameter.Scan.Venc ~= 0))
                MR.Parameter.Recon.ExportRECImgTypes = {'M', 'P'};
            end
            
            auto_update_status = MR.Parameter.Recon.AutoUpdateInfoPars;
            MR.Parameter.Recon.AutoUpdateInfoPars = 'no';
            
            MR.DataClass.Convert2Cell;
            MR.Parameter.UpdateScalingPars;
            
            NY = 'NY';
            
            nr_images = 0;
            for i = 1:size( MR.Data, 2)
                if ~isempty( MR.Data{1,i} )
                    nr_images = nr_images + size( MR.Data{1, i}(:,:,:), 3 );
                end
            end
            % if TKE maps available, then copy them to an extra dataset:
            nr_images2 = 0;
            if (strcmpi(MR.Parameter.Recon.TKE,'Yes')&& (size(MR.Data,2)==2))
                nr_images2 = size( MR.Data{1, 2}(:,:,:), 3 );
            end
            I = InfoPars( nr_images );
            cur_img = 1;
            for i = 1:size( MR.Data, 2)
                if ~isempty( MR.Data{1,i} )
                    % Get the parameters for each image (used for scaling)
                    if iscell( MR.Parameter.ImageInformation )
                        I( cur_img:cur_img+size(MR.Data{1,i}(:,:,:),3 )-1 ) = MR.Parameter.ImageInformation{1,i}(:);
                    else
                        I( cur_img:cur_img+size(MR.Data{1,i}(:,:,:),3 )-1 ) = MR.Parameter.ImageInformation(:);
                    end
                    cur_img = cur_img+size(MR.Data{1,i}(:,:,:),3 );
                end
            end
            
            if isempty(MR.Parameter.Scan.Multivenc) || (strcmpi(MR.Parameter.Scan.Multivenc,'no'))
                if (~isempty(MR.Parameter.Scan.Venc) && sum(MR.Parameter.Scan.Venc ~= 0)>1)
                    nr_files = size(MR.Data{1},10);
                    ParString = cell(1, nr_files);
                    for ffile = 1:nr_files
                        ParString{ffile} = '';
                    end
                else
                    ParString{1} = '';
                end
            else  %multivenc
                venc_dir = sum(MR.Parameter.Scan.Venc>0,1);
                if sum(venc_dir)>0
                    if size(MR.Parameter.Scan.Venc, 1) > 1
                        nr_files = size(MR.Parameter.Scan.Venc, 1);
                    else
                        nr_files = sum(MR.Parameter.Scan.Venc > 0);
                    end
                    
                    ParString = cell(1, nr_files);
                    for nr=1:nr_files %export one magnitude/phase per flow encoded direction
                        ParString{nr} = '';
                    end
                    if nr_images2>0
                        ParString{nr+1} = '';
                    end
                else
                    ParString{1} = '';
                end
            end
            
            v = MRparameter.convert_parameter2output_struct( MR.Parameter  );   
            if nargin == 1 || isempty(Filename)
                Filename = MR.Parameter.Filename.Data;
            end
            [~,dataset_name] = fileparts(Filename); 
            
            venc_ind = abs(MRparameter.coord2num(MRparameter.unformat_coord_str( MR.Parameter.Scan.MPS(1,:) )));
            for i = 1:length(ParString)
                
                if isempty( MR.Parameter.Scan.Venc )
                    MR.Parameter.Scan.Venc = MR.Parameter.Scan.Venc;
                end
                if size(MR.Parameter.Scan.Venc, 2) < 3
                    MR.Parameter.Scan.Venc(:, end+1:3) = 0;
                end
                
                if strcmpi(MR.Parameter.Scan.Multivenc,'yes') && size(MR.Parameter.Scan.Venc, 1) >= length(ParString)
                    venc = MR.Parameter.Scan.Venc(i,:);                    
                else
                    if ~isempty(MR.Parameter.Scan.Venc)
                        venc = v.PhaseEncodingVelocity.*0;
                        if sum(MR.Parameter.Scan.Venc ~= 0)>1
                            try
                                venc(venc_ind(i)) = MR.Parameter.Scan.Venc(i);
                            catch
                                venc = [0,0,0];
                            end
                        else
                            venc = MR.Parameter.Scan.Venc;
                        end
                    else
                        venc = [0,0,0];
                    end
                end
                
                
                ParString{i} = [ParString{i}, sprintf( '# === DATA DESCRIPTION FILE ======================================================\n')];
                ParString{i} = [ParString{i}, sprintf(  '#\n')];
                ParString{i} = [ParString{i}, sprintf(  '# CAUTION - Investigational device.\n')];
                ParString{i} = [ParString{i}, sprintf(  '# Limited by Federal Law to investigational use.\n')];
                ParString{i} = [ParString{i}, sprintf(  '#\n')];
                ParString{i} = [ParString{i}, sprintf(  '# Dataset name: %s\n', dataset_name)];
                ParString{i} = [ParString{i}, sprintf(  '# Exported by MRecon (c)Gyrotools GmbH, Zuerich Switzerland (http://www.gyrotools.ch/)\n')];
                ParString{i} = [ParString{i}, sprintf(  '# CLINICAL TRYOUT             Research image export tool     V4.1\n')];
                ParString{i} = [ParString{i}, sprintf(  '#\n')];
                ParString{i} = [ParString{i}, sprintf(  '# === GENERAL INFORMATION ========================================================\n')];
                ParString{i} = [ParString{i}, sprintf(  '#\n')];
                
                ParString{i} = [ParString{i}, sprintf(  '.    Patient name                       :   %s\n', v.PatientName)];
                ParString{i} = [ParString{i}, sprintf(  '.    Examination name                   :   %s\n', v.ExaminationName)];
                ParString{i} = [ParString{i}, sprintf(  '.    Protocol name                      :   %s\n', v.ProtocolName)];
                ParString{i} = [ParString{i}, sprintf(  '.    Examination date/time              :   %s / %s\n', v.ExaminationDate, v.ExaminationTime )];
                ParString{i} = [ParString{i}, sprintf(  '.    Series_data_type                   :   %d\n',0)];
                ParString{i} = [ParString{i}, sprintf(  '.    Acquisition nr                     :   %d\n',v.AquisitionNumber)];
                ParString{i} = [ParString{i}, sprintf(  '.    Reconstruction nr                  :   %d\n',v.ReconstructionNumber)];
                ParString{i} = [ParString{i}, sprintf(  '.    Scan Duration [sec]                :   %0.2f\n',0)];
                ParString{i} = [ParString{i}, sprintf(  '.    Max. number of cardiac phases      :   %d\n', v.MaxNoPhases )];
                ParString{i} = [ParString{i}, sprintf(  '.    Max. number of echoes              :   %d\n', v.MaxNoEchoes)];
                ParString{i} = [ParString{i}, sprintf(  '.    Max. number of slices/locations    :   %d\n', v.MaxNoSlices)];
                ParString{i} = [ParString{i}, sprintf(  '.    Max. number of dynamics            :   %d\n', v.MaxNoDynamics)];
                ParString{i} = [ParString{i}, sprintf(  '.    Max. number of mixes               :   %d\n', v.MaxNoMixes)];
                ParString{i} = [ParString{i}, sprintf(  '.    Patient Position                   :   %s\n', v.PatientPosition)];
                ParString{i} = [ParString{i}, sprintf(  '.    Preparation direction              :   %s\n', v.PreparationDirection(1,:))];
                ParString{i} = [ParString{i}, sprintf(  '.    Technique                          :   %s\n', v.Technique)];
                ParString{i} = [ParString{i}, sprintf(  '.    Scan resolution  (x, y)            :   %-3d  %3d\n', v.ScanResolutionX, v.ScanResolutionY)];
                ParString{i} = [ParString{i}, sprintf(  '.    Scan mode                          :   %s\n', v.ScanMode)];
                ParString{i} = [ParString{i}, sprintf(  '.    Repetition time [msec]             :   %0.2f\n', v.RepetitionTimes(1) )];
                ParString{i} = [ParString{i}, sprintf(  '.    FOV (ap,fh,rl) [mm]                :   %-5.2f %5.2f %5.2f\n', v.FOVAP, v.FOVFH, v.FOVRL )];
                ParString{i} = [ParString{i}, sprintf(  '.    Water Fat shift [pixels]           :   %0.2f\n',v.WaterFatShift)];
                ParString{i} = [ParString{i}, sprintf(  '.    Angulation midslice(ap,fh,rl)[degr]:   %-6.2f %-6.2f %-6.2f\n', v.AngulationAP, v.AngulationFH, v.AngulationRL )];
                ParString{i} = [ParString{i}, sprintf(  '.    Off Centre midslice(ap,fh,rl) [mm] :   %-6.2f %-6.2f %-6.2f\n', v.OffCenterAP, v.OffCenterFH, v.OffCenterRL)];
                ParString{i} = [ParString{i}, sprintf(  '.    Flow compensation <0=no 1=yes> ?   :   %d\n',strfind( NY, v.FlowCompensation )-1)];
                ParString{i} = [ParString{i}, sprintf(  '.    Presaturation     <0=no 1=yes> ?   :   %d\n',strfind( NY, v.Presaturation )-1 )];
                ParString{i} = [ParString{i}, sprintf(  '.    Phase encoding velocity [cm/sec]   :   %-5.2f %5.2f %5.2f\n',venc(1),venc(2),venc(3))];
                ParString{i} = [ParString{i}, sprintf(  '.    MTC               <0=no 1=yes> ?   :   %d\n',strfind( NY, v.MTC )-1)];
                ParString{i} = [ParString{i}, sprintf(  '.    SPIR              <0=no 1=yes> ?   :   %d\n',strfind( NY, v.SPIR )-1)];
                ParString{i} = [ParString{i}, sprintf(  '.    EPI factor        <0,1=no EPI>     :   %d\n',v.EPIfactor)];
                ParString{i} = [ParString{i}, sprintf(  '.    Dynamic scan      <0=no 1=yes> ?   :   %d\n',strfind( NY, v.DynamicScan )-1)];
                ParString{i} = [ParString{i}, sprintf(  '.    Diffusion         <0=no 1=yes> ?   :   %d\n',strfind( NY, v.Diffusion )-1)];
                ParString{i} = [ParString{i}, sprintf(  '.    Diffusion echo time [msec]         :   %0.2f\n',v.DiffusionEchoTime)];
                ParString{i} = [ParString{i}, sprintf(  '.    Max. number of diffusion values    :   %d\n',v.DiffusionValues)];
                ParString{i} = [ParString{i}, sprintf(  '.    Max. number of gradient orients    :   %d\n',v.MaxNoGradientOrients)];
                ParString{i} = [ParString{i}, sprintf(  '.    Number of label types   <0=no ASL> :   %d\n',v.ASLNolabelTypes)];
                ParString{i} = [ParString{i}, sprintf(  '#\n')];
                
                ParString{i} = [ParString{i}, sprintf(  '# === PIXEL VALUES =============================================================\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  PV = pixel value in REC file, FP = floating point value, DV = displayed value on console\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  RS = rescale slope,           RI = rescale intercept,    SS = scale slope\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  DV = PV * RS + RI             FP = PV /  SS\n')];
                ParString{i} = [ParString{i}, sprintf(  '#\n')];
                ParString{i} = [ParString{i}, sprintf(  '# === IMAGE INFORMATION DEFINITION =============================================\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  The rest of this file contains ONE line per image, this line contains the following information:\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  \n')];
                ParString{i} = [ParString{i}, sprintf(  '#  slice number                             (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  echo number                              (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  dynamic scan number                      (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  cardiac phase number                     (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  image_type_mr                            (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  scanning sequence                        (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  index in REC file (in images)            (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  image pixel size (in bits)               (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  scan percentage                          (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  recon resolution (x,y)                   (2*integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  rescale intercept                        (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  rescale slope                            (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  scale slope                              (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  window center                            (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  window width                             (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  image angulation (ap,fh,rl in degrees )  (3*float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  image offcentre (ap,fh,rl in mm )        (3*float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  slice thickness                          (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  slice gap                                (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  image_display_orientation                (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  slice orientation ( TRA/SAG/COR )        (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  fmri_status_indication                   (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  image_type_ed_es  (end diast/end syst)   (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  pixel spacing (x,y) (in mm)              (2*float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  echo_time                                (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  dyn_scan_begin_time                      (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  trigger_time                             (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  diffusion_b_factor                       (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  number of averages                       (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  image_flip_angle (in degrees)            (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  cardiac frequency                        (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  min. RR. interval                        (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  max. RR. interval                        (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  turbo factor                             (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  inversion delay                          (float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  diffusion b value number    (imagekey!)  (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  gradient orientation number (imagekey!)  (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  contrast type                            (string)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  diffusion anisotropy type                (string)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  diffusion (ap, fh, rl)                   (3*float)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#  label type (ASL)            (imagekey!)  (integer)\n')];
                ParString{i} = [ParString{i}, sprintf(  '#\n')];
                ParString{i} = [ParString{i}, sprintf(  '# === IMAGE INFORMATION ==========================================================\n')];
                ParString{i} = [ParString{i}, sprintf(  '#sl ec dyn ph ty  idx pix %% rec size (re)scale     window       angulation      offcentre         thick  gap   info   spacing   echo  dtime ttime    diff avg  flip  freq RR_int  turbo  delay b grad cont anis diffusion\n\n')];
                
                format_str = '%-4d %-3d %-3d %-3d %-3d %-3d %-5d %-2d %-4d %-4d %-4d %8f %8f %8f %6d %6d %6.2f %6.2f %6.2f %7.3f %7.3f %7.3f %5.2f %5.2f %d %d %d %d %5.3f %5.3f %5.3f %9.3f %5.1f %5d %5d %4.1f %5d %5d %5d %3d %4.2f %3d %3d %6d %6d %8f %8f %8f %3d\n';
            end
            
            total_nr_images = length( v.ImageInformation.Slice );
            
            index = zeros( length( ParString ), 1);
            for i = 1:total_nr_images - nr_images2
                
                if v.ImageInformation.NoData(i)
                    continue;
                end
                
                % Write data to file
                if length( ParString ) > 1
                    cur_fid = v.Extr1(i);
                else
                    cur_fid = 1;
                end
                
                switch v.ImageInformation.Type(i,:)
                    case 'M'
                        type = 0;
                    case 'R'
                        type = 1;
                    case 'I'
                        type = 2;
                    case 'P'
                        type = 3;
                end
                
                ri = v.ImageInformation.RescaleIntercept(i);
                rs = v.ImageInformation.RescaleSlope(i);
                ss = v.ImageInformation.ScaleSlope(i);
                wc = v.ImageInformation.WindowCenter(i);
                ww = v.ImageInformation.WindowWidth(i);
                
                switch v.ImageInformation.Sequence(i,:)
                    case 'FFE'
                        seq = 0;
                    case 'SE'
                        seq = 1;
                end
                switch strtrim(v.ImageInformation.SliceOrientation(i,:))
                    case 'Transversal'
                        ori = 1;
                    case 'Coronal'
                        ori = 3;
                    case 'Sagital'
                        ori = 2;
                end
                
                ParString{cur_fid} = [ParString{cur_fid}, sprintf( format_str, ...
                    v.ImageInformation.Slice(i), ...                slice number
                    v.ImageInformation.Echo(i), ...                 echo number
                    v.ImageInformation.Dynamic(i), ...              dynamic scan number
                    v.ImageInformation.Phase(i), ...                cardiac phase number
                    type, ...                                       image_type_mr
                    seq, ...                                        scanning sequence
                    index( cur_fid ), ...                           index in REC file (in images)
                    16, ...                                         image pixel size (in bits)
                    v.ImageInformation.ScanPercentage(i), ...       scan percentage
                    v.ImageInformation.ResolutionX(i), ...          recon resolution (x,y)
                    v.ImageInformation.ResolutionY(i), ...          recon resolution (x,y)
                    ri, ...                                         rescale intercept
                    rs, ...                                         rescale slope
                    ss, ...                                         scale slope
                    wc,...                                          window center
                    ww,...                                          window width
                    v.ImageInformation.AngulationAP(i),...          image angulation (ap,fh,rl in degrees )
                    v.ImageInformation.AngulationFH(i),...          image angulation (ap,fh,rl in degrees )
                    v.ImageInformation.AngulationRL(i), ...         image angulation (ap,fh,rl in degrees )
                    v.ImageInformation.OffcenterAP(i), ...          image offcentre (ap,fh,rl in mm )
                    v.ImageInformation.OffcenterFH(i),...           image offcentre (ap,fh,rl in mm )
                    v.ImageInformation.OffcenterRL(i), ...          image offcentre (ap,fh,rl in mm )
                    v.ImageInformation.SliceThickness(i), ...       slice thickness
                    v.ImageInformation.SliceGap(i), ...             slice gap
                    0, ...                                          image_display_orientation
                    ori, ...                                        slice orientation ( TRA/SAG/COR )
                    v.ImageInformation.fMRIStatusIndication(i), ... fmri_status_indication
                    0, ...                                          image_type_ed_es  (end diast/end syst)
                    v.ImageInformation.PixelSpacing(i,1), ...       pixel spacing (x,y) (in mm)
                    v.ImageInformation.PixelSpacing(i,1), ...       pixel spacing (x,y) (in mm)
                    v.ImageInformation.EchoTime(i), ...             echo_time
                    v.ImageInformation.DynScanBeginTime(i),...      dyn_scan_begin_time
                    v.ImageInformation.TriggerTime(i),...           trigger_time
                    v.ImageInformation.DiffusionBFactor(i),...      diffusion_b_factor
                    v.ImageInformation.NoAverages(i), ...           number of averages
                    v.ImageInformation.ImageFlipAngle(i), ...       image_flip_angle (in degrees)
                    v.ImageInformation.CardiacFrequency(i),...      cardiac frequency
                    v.ImageInformation.MinRRInterval(i), ...        min. RR. interval
                    v.ImageInformation.MaxRRInterval(i), ...        max. RR. interval
                    v.ImageInformation.TURBOFactor(i), ...          turbo factor
                    v.ImageInformation.InversionDelay(i), ...       inversion delay
                    v.ImageInformation.BValue(i), ...               diffusion b value number    (imagekey!)
                    v.ImageInformation.GradOrient(i), ...           gradient orientation number (imagekey!)
                    0, ...                                          contrast type
                    0, ...                                          diffusion anisotropy type
                    v.ImageInformation.DiffusionAP(i), ...          diffusion (ap, fh, rl)
                    v.ImageInformation.DiffusionFH(i), ...          diffusion (ap, fh, rl)
                    v.ImageInformation.DiffusionRL(i), ...          diffusion (ap, fh, rl)
                    v.ImageInformation.LabelTypeASL(i) )];           %label type
                
                index( cur_fid ) = index( cur_fid )+1;
            end
            %then write TKE:
            if ( (nr_images2 > 2) && strcmpi(MR.Parameter.Recon.TKE,'Yes'))
                cur_fid = cur_fid + 1;
                for i = total_nr_images - nr_images2 + 1:total_nr_images
                    
                    
                    type = 0;
                    
                    ri = v.ImageInformation.RescaleIntercept(i);
                    rs = v.ImageInformation.RescaleSlope(i);
                    ss = v.ImageInformation.ScaleSlope(i);
                    wc = v.ImageInformation.WindowCenter(i);
                    ww = v.ImageInformation.WindowWidth(i);
                    
                    switch v.ImageInformation.Sequence(i,:)
                        case 'FFE'
                            seq = 0;
                        case 'SE'
                            seq = 1;
                    end
                    switch strtrim(v.ImageInformation.SliceOrientation(i,:))
                        case 'Transversal'
                            ori = 1;
                        case 'Coronal'
                            ori = 3;
                        case 'Sagital'
                            ori = 2;
                    end
                    
                    ParString{cur_fid} = [ParString{cur_fid}, sprintf(format_str, ...
                        v.ImageInformation.Slice(i), ...                slice number
                        v.ImageInformation.Echo(i), ...                 echo number
                        v.ImageInformation.Dynamic(i), ...              dynamic scan number
                        v.ImageInformation.Phase(i), ...                cardiac phase number
                        type, ...                                       image_type_mr
                        seq, ...                                        scanning sequence
                        index( cur_fid ), ...                           index in REC file (in images)
                        16, ...                                         image pixel size (in bits)
                        v.ImageInformation.ScanPercentage(i), ...       scan percentage
                        v.ImageInformation.ResolutionX(i), ...          recon resolution (x,y)
                        v.ImageInformation.ResolutionY(i), ...          recon resolution (x,y)
                        ri, ...                                         rescale intercept
                        rs, ...                                         rescale slope
                        ss, ...                                         scale slope
                        wc,...                                          window center
                        ww,...                                          window width
                        v.ImageInformation.AngulationAP(i),...          image angulation (ap,fh,rl in degrees )
                        v.ImageInformation.AngulationFH(i),...          image angulation (ap,fh,rl in degrees )
                        v.ImageInformation.AngulationRL(i), ...         image angulation (ap,fh,rl in degrees )
                        v.ImageInformation.OffcenterAP(i), ...          image offcentre (ap,fh,rl in mm )
                        v.ImageInformation.OffcenterFH(i),...           image offcentre (ap,fh,rl in mm )
                        v.ImageInformation.OffcenterRL(i), ...          image offcentre (ap,fh,rl in mm )
                        v.ImageInformation.SliceThickness(i), ...       slice thickness
                        v.ImageInformation.SliceGap(i), ...             slice gap
                        0, ...                                          image_display_orientation
                        ori, ...                                        slice orientation ( TRA/SAG/COR )
                        v.ImageInformation.fMRIStatusIndication(i), ... fmri_status_indication
                        0, ...                                          image_type_ed_es  (end diast/end syst)
                        v.ImageInformation.PixelSpacing(i,1), ...       pixel spacing (x,y) (in mm)
                        v.ImageInformation.PixelSpacing(i,1), ...       pixel spacing (x,y) (in mm)
                        v.ImageInformation.EchoTime(i), ...             echo_time
                        v.ImageInformation.DynScanBeginTime(i),...      dyn_scan_begin_time
                        v.ImageInformation.TriggerTime(i),...           trigger_time
                        v.ImageInformation.DiffusionBFactor(i),...      diffusion_b_factor
                        v.ImageInformation.NoAverages(i), ...           number of averages
                        v.ImageInformation.ImageFlipAngle(i), ...       image_flip_angle (in degrees)
                        v.ImageInformation.CardiacFrequency(i),...      cardiac frequency
                        v.ImageInformation.MinRRInterval(i), ...        min. RR. interval
                        v.ImageInformation.MaxRRInterval(i), ...        max. RR. interval
                        v.ImageInformation.TURBOFactor(i), ...          turbo factor
                        v.ImageInformation.InversionDelay(i), ...       inversion delay
                        v.ImageInformation.BValue(i), ...               diffusion b value number    (imagekey!)
                        v.ImageInformation.GradOrient(i), ...           gradient orientation number (imagekey!)
                        0, ...                                          contrast type
                        0, ...                                          diffusion anisotropy type
                        v.ImageInformation.DiffusionAP(i), ...          diffusion (ap, fh, rl)
                        v.ImageInformation.DiffusionFH(i), ...          diffusion (ap, fh, rl)
                        v.ImageInformation.DiffusionRL(i), ...          diffusion (ap, fh, rl)
                        v.ImageInformation.LabelTypeASL(i) )];           %label type
                    
                    index( cur_fid ) = index( cur_fid )+1;
                end
            end
            
            for i = 1:length(ParString)
                ParString{i} = [ParString{i}, newline];
                ParString{i} = [ParString{i}, sprintf(  '# === END OF DATA DESCRIPTION FILE ===============================================\n')];
            end
            
            MR.Data = Helper.UnconvertCell( MR.Data );
            MR.Parameter.Recon.AutoUpdateInfoPars = auto_update_status;
        end
        function RecImages = GetRecImages( MR )
            % Check if the data is empty                   
            if isempty( MR.Data )
                error( 'Error in WriteRec: The data matrix is empty. Please reconstruct the data first');
            end
            
            % if the data is not a cell yet convert it to one
            MR.DataClass.Convert2Cell;
            
            % Check if ImageInformation is empty. If it is create one
            if isempty( MR.Parameter.ImageInformation )
                info_dim = cell(size(MR.Data,1), size(MR.Data, 2));
                for ci = 1:size(MR.Data,1)
                    for cj = 1:size(MR.Data, 2)
                        info_dim{ci, cj} = size( MR.Data{ci,cj} );
                        if sum(info_dim{ci,cj}) ~= 0
                            if length( info_dim{ci,cj} ) > 2
                                info_dim{ci,cj} = info_dim{ci,cj}(3:end);
                            else
                                info_dim{ci,cj} = 1;
                            end
                        else
                            info_dim{ci,cj} = [];
                        end
                        if length( info_dim{ci, cj} ) == 1
                            info_dim{ci, cj} = [info_dim{ci, cj}, 1];
                        end
                    end
                end
                info_dim = Helper.UnconvertCell(info_dim);
                MR.Parameter.CreateInfoPars( info_dim );                
            end
            
            if any( any( MR.Parameter.Scan.Venc ~= 0))
                MR.Parameter.Recon.ExportRECImgTypes = {'M', 'P'};
            end                        
            
            % switch off the automatic update of the image information
            % parameter
            MR.Parameter.UpdateImageInfo = 0;
            
            % the images which are written to the recfile have to be
            % square. Furthermore their matrix sizes have to be equal.
            % Therefore here we make it work for cells and different image sizes
            if size( MR.Data, 1) > 1
                warning( 'MATLAB:MRecon', 'Only standard data is written to the .rec file');
            end
            if size( MR.Data, 2) > 1
                if any( cellfun( @(x) size(x,1) ~= size(x,2), MR.Data(1,:) ))
                    warning( 'MATLAB:MRecon', 'There are different sized images in the data. The sizes are made equal before writing');
                end
            end
            
            % calculate the total number of images in r.Data (if data is a
            % cell the look in all cell elements)
            nr_images = 0;
            for i = 1:size( MR.Data, 2)
                if ~isempty( MR.Data{1,i} )
                    nr_images = nr_images + size( MR.Data{1, i}(:,:,:), 3 );
                end
            end
            % if TKE maps available, then copy them to an extra dataset:
            nr_images2 = 0;
            if (strcmpi(MR.Parameter.Scan.Multivenc,'Yes')&& (size(MR.Data,2)==2))
                nr_images2 = size( MR.Data{1, 2}(:,:,:), 3 );
            end
            
            MR.Parameter.UpdateScalingPars;
            
            % define the data matrix which is written to the recfile
            %             data = zeros(MR.Parameter.Encoding.WorkEncoding.XReconRes{1,1}, MR.Parameter.Encoding.WorkEncoding.YReconRes{1,1}, nr_images );
            data = zeros(size(MR.Data{1,1},1), size(MR.Data{1,1},2), nr_images );
            I = InfoPars( size(data,3) );
            cur_img = 1;
            % if the data is a cell then check if the images in the
            % different cell elements are of equal size. If not make them
            % equal by zero-padding.
            % Also make all image sizes equal to:
            % [MR.Parameter.Encoding.WorkEncoding.XReconRes, MR.Parameter.Encoding.WorkEncoding.YReconRes]
            for i = 1:size( MR.Data, 2)
                if ~isempty( MR.Data{1,i} )
                    ox = 1;
                    oy = 1;
                    oz = 1;
                    if ~isempty( MR.Parameter.Encoding.WorkEncoding.XReconRes{1,1}) && ...
                            size(MR.Data{1,i},1) ~= size(data,1)
                        ox = size(MR.Data{1,i},1) / size(data,1);
                    end
                    if ~isempty( MR.Parameter.Encoding.WorkEncoding.YReconRes{1,1}) && ...
                            size(MR.Data{1,i},2) ~= size(data,2)
                        oy = size(MR.Data{1,i},2) / size(data,2);
                    end
                    
                    if any( [ox, oy, oz] ~= [1,1,1] )
                        temp = zeros( size(data,1), size(data,2), size(MR.Data{1,i},3) );
                        for j = 1:size(MR.Data{1,i},3)
                            temp(:,:,j) = complex( imresize( real(MR.Data{1,i}(:,:,j)) , [size(data,1), size(data,2)] ), ...
                                imresize( imag(MR.Data{1,i}(:,:,j)) , [size(data,1), size(data,2)] ) );
                        end
                        data(:,:,cur_img:cur_img+size(MR.Data{1,i}(:,:,:),3 )-1) = temp(:,:,:);
                    else
                        data(:,:,cur_img:cur_img+size(MR.Data{1,i}(:,:,:),3 )-1) = MR.Data{1,i}(:,:,:);
                    end
                    
                    % Get the parameters for each image (used for scaling)
                    if iscell( MR.Parameter.ImageInformation )
                        I( cur_img:cur_img+size(MR.Data{1,i}(:,:,:),3 )-1 ) = MR.Parameter.ImageInformation{1,i}(:);
                    else
                        I( cur_img:cur_img+size(MR.Data{1,i}(:,:,:),3 )-1 ) = MR.Parameter.ImageInformation(:);
                    end
                    cur_img = cur_img+size(MR.Data{1,i}(:,:,:),3 );
                end
            end
            
            if isempty(MR.Parameter.Scan.Multivenc) || strcmpi(MR.Parameter.Scan.Multivenc,'no')
                if (~isempty(MR.Parameter.Scan.Venc) && sum(MR.Parameter.Scan.Venc ~= 0)>1)
                    nr_files = size(MR.Data{1},10);
                    RecImages = cell(1, nr_files);
                    for ffile = 1:nr_files
                        RecImages{ffile} = [];
                    end
                else
                    RecImages{1} = [];
                end
            else %multivenc
                venc_dir = sum(MR.Parameter.Scan.Venc>0,1);
                if sum(venc_dir)>0
                    if size(MR.Parameter.Scan.Venc, 1) > 1
                        nr_files = size(MR.Parameter.Scan.Venc, 1);
                    else
                        nr_files = sum(MR.Parameter.Scan.Venc > 0);
                    end
                    
                    RecImages = cell(1, nr_files);
                    for nr=1:nr_files %export one magnitude/phase per flow encoded direction
                        RecImages{nr} = [];
                    end
                    if nr_images2>0 %tke
                        RecImages{nr+1} = [];
                    end
                else
                    RecImages{1} = [];
                end
            end
            
            
            nr_loops = length( MR.Parameter.Recon.ExportRECImgTypes );
            
            % determine the number of images per parfile
            nr_images_per_par = zeros(length( RecImages ), 1);
            for phase_loop = 1:nr_loops
                for k = 1:(size( data, 3) - nr_images2)
                    if I(k).NoData
                        continue;
                    end
                    cur_fid = min( [ length( RecImages ), I(k).Extra1] );
                    nr_images_per_par(cur_fid) = nr_images_per_par(cur_fid)+1;
                end
            end
            for i = 1:length(RecImages)
                RecImages{i} = zeros(size(data,2), size(data, 1), nr_images_per_par(i), 'uint16');
            end
            
            
            cur_image_per_par = zeros(length( RecImages ), 1);
            for phase_loop = 1:nr_loops                                
                for k = 1:(size( data, 3) - nr_images2)
                    if I(k).NoData
                        continue;
                    end
                    
                    cur_fid = min( [ length( RecImages ), I(k).Extra1] );
                    cur_image_per_par(cur_fid) = cur_image_per_par(cur_fid) + 1;
                    
                    % Scale the images according to the Rescale values in
                    % the image information struct.
                    % Only write the specified image type to the recfile
                    % (abs, angle, real, imag)
                    switch MR.Parameter.Recon.ExportRECImgTypes{ phase_loop }
                        case 'M'
                            try
                                RecImages{cur_fid}(:,:,cur_image_per_par(cur_fid)) = uint16(round( (abs( data(:,:,k) ) - I(k).RescaleIntercept.M) ./ I(k).RescaleSlope.M)');
                            catch
                                RecImages{cur_fid}(:,:,cur_image_per_par(cur_fid)) = uint16(round( (abs( data(:,:,k) ) - I(k).RescaleIntercept(1)) ./ I(k).RescaleSlope(1))');
                            end
                        case 'P'
                            try
                                RecImages{cur_fid}(:,:,cur_image_per_par(cur_fid)) = uint16(round( (floor(1000.*angle( data(:,:,k) )) - I(k).RescaleIntercept.P) ./ I(k).RescaleSlope.P)');
                            catch
                                RecImages{cur_fid}(:,:,cur_image_per_par(cur_fid)) = uint16(round( (floor(1000.*angle( data(:,:,k) )) - I(k).RescaleIntercept(2)) ./ I(k).RescaleSlope(2))');
                            end
                        case 'R'
                            try
                                RecImages{cur_fid}(:,:,cur_image_per_par(cur_fid)) = uint16(round( (real( data(:,:,k) ) - I(k).RescaleIntercept.R) ./ I(k).RescaleSlope.R)');
                            catch
                                RecImages{cur_fid}(:,:,cur_image_per_par(cur_fid)) = uint16(round( (real( data(:,:,k) ) - I(k).RescaleIntercept(3)) ./ I(k).RescaleSlope(3))');
                            end
                        case 'I'
                            try
                                RecImages{cur_fid}(:,:,cur_image_per_par(cur_fid)) = uint16(round( (imag( data(:,:,k) ) - I(k).RescaleIntercept.I) ./ I(k).RescaleSlope.I)');
                            catch
                                RecImages{cur_fid}(:,:,cur_image_per_par(cur_fid)) = uint16(round( (imag( data(:,:,k) ) - I(k).RescaleIntercept(4)) ./ I(k).RescaleSlope.(4))');
                            end
                    end                    
                end
            end
            if (nr_images2 > 0)                               
                loop = 1;
                for k = (size( data, 3) - nr_images2 + 1):size(data,3)
                    
                    % Scale the images according to the Rescale values in
                    % the image information struct.
                    % Only write abs for TKE
                    try
                        RecImages{cur_fid+1}(:,:,loop) = uint16(round( (abs( data(:,:,k) ) - I(k).RescaleIntercept.M) ./ I(k).RescaleSlope.M)');
                    catch
                        RecImages{cur_fid+1}(:,:,loop) = uint16(round( (abs( data(:,:,k) ) - I(k).RescaleIntercept(1)) ./ I(k).RescaleSlope(1))');
                    end                    
                    loop = loop+1;
                end
            end
            
            % Restore the original data matrix
            MR.Data = Helper.UnconvertCell( MR.Data );
            MR.Parameter.UpdateImageInfo = 1;
            
        end
        function WriteXMLPar( MR, Filename )
            % WriteXMLPar: Exports the parameter to a Philips XML par file.
            %
            % Syntax:     r.WriteXMLPar( filename );
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ScanParsDoc')">Parameter.Scan</a>:
            %             Holds the general parameters (upper part of the par file)
            %           - <a href="matlab:helpwin('InfoParsDoc')">Parameter.ImageInformation</a>:
            %             Holds the individual parameters for every image (lower part of the par
            %             file).
            %
            % Location:   image-space.
            %
            % Formats:    Raw | ExportedRaw | Cpx | Rec
            %
            % Description/Algorithm: Exports the parameters in to a Philips XML par file.
            
            auto_update_status = MR.Parameter.Recon.AutoUpdateInfoPars;
            MR.Parameter.Recon.AutoUpdateInfoPars = 'no';
            
            if nargin == 1 || isempty( Filename )
                Filename = 'recon.xml';
            end
            
            MR.DataClass.Convert2Cell;
            MR.Parameter.UpdateScalingPars;
            
            nr_images = 0;
            for i = 1:size( MR.Data, 2)
                if ~isempty( MR.Data{1,i} )
                    nr_images = nr_images + size( MR.Data{1, i}(:,:,:), 3 );
                end
            end
            I = InfoPars( nr_images );
            cur_img = 1;
            for i = 1:size( MR.Data, 2)
                if ~isempty( MR.Data{1,i} )
                    % Get the parameters for each image (used for scaling)
                    if iscell( MR.Parameter.ImageInformation )
                        I( cur_img:cur_img+size(MR.Data{1,i}(:,:,:),3 )-1 ) = MR.Parameter.ImageInformation{1,i}(:);
                    else
                        I( cur_img:cur_img+size(MR.Data{1,i}(:,:,:),3 )-1 ) = MR.Parameter.ImageInformation(:);
                    end
                    cur_img = cur_img+size(MR.Data{1,i}(:,:,:),3 );
                end
            end
            
            v = MRparameter.convert_parameter2output_struct( MR.Parameter  );
            if all(MR.Parameter.Scan.Venc ~= 0) && length(unique([I.Extra1])) == 3
                if  contains(Filename, '.xml' ) 
                    fid{1} = fopen( [strrep(Filename, '.xml', ''),'M.xml'],'w');
                    fid{2} = fopen( [strrep(Filename, '.xml', ''),'P.xml'],'w');
                    fid{3} = fopen( [strrep(Filename, '.xml', ''),'S.xml'],'w');
                elseif  contains(Filename, '.XML' ) 
                    fid{1} = fopen( [strrep(Filename, '.XML', ''),'M.XML'],'w');
                    fid{2} = fopen( [strrep(Filename, '.XML', ''),'P.XML'],'w');
                    fid{3} = fopen( [strrep(Filename, '.XML', ''),'S.XML'],'w');
                else
                    fid{1} = fopen( [Filename,'M.XML'],'w');
                    fid{2} = fopen( [Filename,'P.XML'],'w');
                    fid{3} = fopen( [Filename,'S.XML'],'w');
                end
            else
                fid = fopen(Filename,'w');
            end
            
            fprintf( fid, '<PRIDE_V5>\r\n' );
            fprintf( fid, '<Series_Info>\r\n' );
            
            fprintf( fid, '<Attribute Name="Patient Name" Tag="0x00100010" Level="Patient" Type="String">%s</Attribute>\r\n', strtrim(v.PatientName) );
            fprintf( fid, '<Attribute Name="Examination Name" Tag="0x00400254" Level="Examination" Type="String">%s</Attribute>\r\n', strtrim(v.ExaminationName) );
            fprintf( fid, '<Attribute Name="Protocol Name" Tag="0x00181030" Level="MRSeries" Type="String">%s</Attribute>\r\n', strtrim(v.ProtocolName) );
            fprintf( fid, '<Attribute Name="Examination Date" Tag="0x00400244" Level="Examination" Type="Date">%s</Attribute>\r\n', strtrim(v.ExaminationDate) );
            fprintf( fid, '<Attribute Name="Examination Time" Tag="0x00400245" Level="Examination" Type="Time">%s</Attribute>\r\n', strtrim(v.ExaminationTime) );
            fprintf( fid, '<Attribute Name="Series Data Type" Tag="0x20051035" Level="MRSeries" Type="String">%s</Attribute>\r\n',  strtrim(v.SeriesDataType));
            fprintf( fid, '<Attribute Name="Aquisition Number" Tag="0x2001107B" Level="MRSeries" Type="Int32">%d</Attribute>\r\n', v.AquisitionNumber);
            fprintf( fid, '<Attribute Name="Reconstruction Number" Tag="0x2001101D" Level="MRSeries" Type="Int32">%d</Attribute>\r\n', v.ReconstructionNumber);
            fprintf( fid, '<Attribute Name="Scan Duration" Tag="0x2001101D" Level="MRSeries" Type="Float">%f</Attribute>\r\n', v.ScanDuration);
            fprintf( fid, '<Attribute Name="Max No Phases" Tag="0x20011017" Level="MRSeries" Type="Int32">%d</Attribute>\r\n', v.MaxNoPhases);
            fprintf( fid, '<Attribute Name="Max No Echoes" Tag="0x20011014" Level="MRSeries" Type="Int32">%d</Attribute>\r\n', v.MaxNoEchoes );
            fprintf( fid, '<Attribute Name="Max No Slices" Tag="0x20011018" Level="MRSeries" Type="Int32">%d</Attribute>\r\n', v.MaxNoSlices );
            fprintf( fid, '<Attribute Name="Max No Dynamics" Tag="0x20011081" Level="MRSeries" Type="Int32">%d</Attribute>\r\n', v.MaxNoDynamics );
            fprintf( fid, '<Attribute Name="Max No Mixes" Tag="0x20051021" Level="MRSeries" Type="Int16">%d</Attribute>\r\n', v.MaxNoMixes );
            fprintf( fid, '<Attribute Name="Max No B Values" Tag="0x20051414" Level="MRSeries" Type="Int32">%d</Attribute>\r\n', v.MaxNoBValues );
            fprintf( fid, '<Attribute Name="Max No Gradient Orients" Tag="0x20051415" Level="MRSeries" Type="Int32">%d</Attribute>\r\n', v.MaxNoGradientOrients );
            fprintf( fid, '<Attribute Name="No Label Types" Tag="0x20051428" Level="MRSeries" Type="Int32">%d</Attribute>\r\n', v.NoLabelTypes );
            fprintf( fid, '<Attribute Name="Patient Position" Tag="0x00185100" Level="MRSeries" Type="String">%s</Attribute>\r\n', strtrim(v.PatientPosition) );
            fprintf( fid, '<Attribute Name="Preparation Direction" Tag="0x2005107B" Level="MRStack" Type="String">%s</Attribute>\r\n', strtrim(v.PreparationDirection) );
            fprintf( fid, '<Attribute Name="Technique" Tag="0x20011020" Level="MRSeries" Type="String">%s</Attribute>\r\n', strtrim(v.Technique) );
            fprintf( fid, '<Attribute Name="Scan Resolution X" Tag="0x2005101D" Level="MRSeries" Type="Int16">%d</Attribute>\r\n', v.ScanResolutionX );
            fprintf( fid, '<Attribute Name="Scan Resolution Y" Tag="0x00180089" Level="MRImage" Type="Int32">%d</Attribute>\r\n', v.ScanResolutionY );
            fprintf( fid, '<Attribute Name="Scan Mode" Tag="0x2005106F" Level="MRSeries" Type="String">%s</Attribute>\r\n', strtrim(v.ScanMode) );
            if length( v.RepetitionTimes ) > 1
                temp1 = sprintf( '%2.4E', v.RepetitionTimes(1) );
                temp2 = sprintf( '%2.4E', v.RepetitionTimes(2) );
                temp1( end-2 ) = [];
                temp2( end-2 ) = [];
                fprintf( fid, '<Attribute Name="Repetition Times" Tag="0x20051030" Level="MRSeries" Type="Float" ArraySize="2">%s %s</Attribute>\r\n', temp1, temp2  );
            else
                temp1 = sprintf( '%2.4E', v.RepetitionTimes(1) );
                temp1( end-2 ) = [];
                fprintf( fid, '<Attribute Name="Repetition Times" Tag="0x20051030" Level="MRSeries" Type="Float" ArraySize="2">%s</Attribute>\r\n', temp1 );
            end
            temp1 = sprintf( '%2.4E', v.FOVAP );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="FOV AP" Tag="0x20051074" Level="MRStack" Type="Float">%s</Attribute>\r\n', temp1);
            temp1 = sprintf( '%2.4E', v.FOVFH );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="FOV FH" Tag="0x20051075" Level="MRStack" Type="Float">%s</Attribute>\r\n', temp1 );
            temp1 = sprintf( '%2.4E', v.FOVRL );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="FOV RL" Tag="0x20051076" Level="MRStack" Type="Float">%s</Attribute>\r\n', temp1 );
            temp1 = sprintf( '%2.4E', v.WaterFatShift );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="Water Fat Shift" Tag="0x20011022" Level="MRSeries" Type="Float">%s</Attribute>\r\n', temp1 );
            temp1 = sprintf( '%2.4E', v.AngulationAP );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="Angulation AP" Tag="0x20051071" Level="MRStack" Type="Float">%s</Attribute>\r\n', temp1 );
            temp1 = sprintf( '%2.4E', v.AngulationFH );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="Angulation FH" Tag="0x20051072" Level="MRStack" Type="Float">%s</Attribute>\r\n', temp1 );
            temp1 = sprintf( '%2.4E', v.AngulationRL );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="Angulation RL" Tag="0x20051073" Level="MRStack" Type="Float">%s</Attribute>\r\n', temp1 );
            temp1 = sprintf( '%2.4E', v.OffCenterAP );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="Off Center AP" Tag="0x20051078" Level="MRStack" Type="Float">%s</Attribute>\r\n', temp1 );
            temp1 = sprintf( '%2.4E', v.OffCenterFH );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="Off Center FH" Tag="0x20051079" Level="MRStack" Type="Float">%s</Attribute>\r\n', temp1 );
            temp1 = sprintf( '%2.4E', v.OffCenterRL );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="Off Center RL" Tag="0x2005107A" Level="MRStack" Type="Float">%s</Attribute>\r\n', temp1 );
            fprintf( fid, '<Attribute Name="Flow Compensation" Tag="0x20051016" Level="MRSeries" Type="Boolean">%s</Attribute>\r\n', strtrim(v.FlowCompensation) );
            fprintf( fid, '<Attribute Name="Presaturation" Tag="0x2005102F" Level="MRSeries" Type="Boolean">%s</Attribute>\r\n', strtrim(v.Presaturation) );
            temp1 = sprintf( '%2.4E', v.PhaseEncodingVelocity(1) );
            temp2 = sprintf( '%2.4E', v.PhaseEncodingVelocity(2) );
            temp3 = sprintf( '%2.4E', v.PhaseEncodingVelocity(3) );
            temp1( end-2 ) = [];
            temp2( end-2 ) = [];
            temp3( end-2 ) = [];
            fprintf( fid, '<Attribute Name="Phase Encoding Velocity" Tag="0x2001101A" Level="MRSeries" Type="Float" ArraySize="3">%s %s %s</Attribute>\r\n', temp1, temp2, temp3 );
            fprintf( fid, '<Attribute Name="MTC" Tag="0x2005101C" Level="MRSeries" Type="Boolean">%s</Attribute>\r\n', strtrim(v.MTC) );
            fprintf( fid, '<Attribute Name="SPIR" Tag="0x20011021" Level="MRSeries" Type="Boolean">%s</Attribute>\r\n', strtrim(v.SPIR) );
            fprintf( fid, '<Attribute Name="EPI factor" Tag="0x20011013" Level="MRSeries" Type="Int32">%d</Attribute>\r\n', v.EPIfactor );
            fprintf( fid, '<Attribute Name="Dynamic Scan" Tag="0x20011012" Level="MRSeries" Type="Boolean">%s</Attribute>\r\n', strtrim(v.DynamicScan) );
            fprintf( fid, '<Attribute Name="Diffusion" Tag="0x20051014" Level="MRSeries" Type="Boolean">%s</Attribute>\r\n', strtrim(v.Diffusion) );
            temp1 = sprintf( '%2.4E', v.DiffusionEchoTime );
            temp1( end-2 ) = [];
            fprintf( fid, '<Attribute Name="Diffusion Echo Time" Tag="0x20011011" Level="MRSeries" Type="Float">%s</Attribute>\r\n', temp1 );
            fprintf( fid, '</Series_Info>\r\n' );
            
            fprintf( fid, '<Image_Array>' );
            
            for i = 1:length( v.ImageInformation.Slice )
                if v.ImageInformation.NoData(i)
                    continue;
                end
                
                fprintf( fid,'<Image_Info>');
                fprintf( fid,'<Key>');
                fprintf( fid,'<Attribute Name="Slice" Tag="0x2001100A" Type="Int32">%d</Attribute>\r\n', v.ImageInformation.Slice(i) );
                fprintf( fid,'<Attribute Name="Echo" Tag="0x00180086" Type="Int32">%d</Attribute>\r\n', v.ImageInformation.Echo(i) );
                fprintf( fid,'<Attribute Name="Dynamic" Tag="0x00200100" Type="Int32">%d</Attribute>\r\n', v.ImageInformation.Dynamic(i) );
                fprintf( fid,'<Attribute Name="Phase" Tag="0x20011008" Type="Int32">%d</Attribute>\r\n', v.ImageInformation.Phase(i) );
                fprintf( fid,'<Attribute Name="BValue" Tag="0x20051412" Type="Int32">%d</Attribute>\r\n', v.ImageInformation.BValue(i) );
                fprintf( fid,'<Attribute Name="Grad Orient" Tag="0x20051413" Type="Int32">%d</Attribute>\r\n', v.ImageInformation.GradOrient(i) );
                fprintf( fid,'<Attribute Name="Label Type" Tag="0x20051429" Type="Enumeration" EnumType="Label_Type">%s</Attribute>\r\n', strtrim(v.ImageInformation.LabelType(i,:)) );
                fprintf( fid,'<Attribute Name="Type" Tag="0x20051011" Type="Enumeration" EnumType="Image_Type">%s</Attribute>\r\n', strtrim(v.ImageInformation.Type(i,:)));
                fprintf( fid,'<Attribute Name="Sequence" Tag="0x2005106E" Type="Enumeration" EnumType="Image_Sequence">%s</Attribute>\r\n', strtrim(v.ImageInformation.Sequence(i,:)) );
                fprintf( fid,'<Attribute Name="Index" Type="Int32" Calc="Index">%d</Attribute>\r\n', v.ImageInformation.Index(i) );
                fprintf( fid,'</Key>\r\n' );
                fprintf( fid,'<Attribute Name="Pixel Size" Tag="0x00280100" Type="UInt16">%d</Attribute>\r\n', v.ImageInformation.PixelSize(i) );
                fprintf( fid,'<Attribute Name="Scan Percentage" Tag="0x00180093" Type="Double">%2.6E</Attribute>\r\n', v.ImageInformation.ScanPercentage(i) );
                fprintf( fid,'<Attribute Name="Resolution X" Tag="0x00280011" Type="UInt16">%d</Attribute>\r\n', v.ImageInformation.ResolutionX(i) );
                fprintf( fid,'<Attribute Name="Resolution Y" Tag="0x00280010" Type="UInt16">%d</Attribute>\r\n', v.ImageInformation.ResolutionY(i) );
                
                ri = v.ImageInformation.RescaleIntercept(i);
                rs = v.ImageInformation.RescaleSlope(i);
                ss = v.ImageInformation.ScaleSlope(i);
                wc = v.ImageInformation.WindowCenter(i);
                ww = v.ImageInformation.WindowWidth(i);
                
                fprintf( fid,'<Attribute Name="Rescale Intercept" Tag="0x00281052" Type="Double">%2.6E</Attribute>\r\n', ri );
                fprintf( fid,'<Attribute Name="Rescale Slope" Tag="0x00281053" Type="Double">%2.6E</Attribute>\r\n', rs );
                temp1 = sprintf( '%2.4E', ss );
                temp1( end-2 ) = [];
                fprintf( fid,'<Attribute Name="Scale Slope" Tag="0x2005100E" Type="Float">%s</Attribute>\r\n', temp1 );
                fprintf( fid,'<Attribute Name="Window Center" Tag="0x00281050" Type="Double">%2.6E</Attribute>\r\n', wc );
                fprintf( fid,'<Attribute Name="Window Width" Tag="0x00281051" Type="Double">%2.6E</Attribute>\r\n', ww );
                fprintf( fid,'<Attribute Name="Slice Thickness" Tag="0x00180050" Type="Double">%2.6E</Attribute>\r\n', v.ImageInformation.SliceThickness(i) );
                fprintf( fid,'<Attribute Name="Slice Gap" Tag="0x00180088" Type="Double">%2.6E</Attribute>\r\n', v.ImageInformation.SliceGap(i) );
                fprintf( fid,'<Attribute Name="Display Orientation" Tag="0x20051004" Type="Enumeration" EnumType="Display_Orientation">%s</Attribute>\r\n', strtrim(v.ImageInformation.DisplayOrientation(i,:)) );
                fprintf( fid,'<Attribute Name="fMRI Status Indication" Tag="0x20051063" Type="Int16">%d</Attribute>\r\n', v.ImageInformation.fMRIStatusIndication(i) );
                fprintf( fid,'<Attribute Name="Image Type Ed Es" Tag="0x20011007" Type="Enumeration" EnumType="Type_ed_es">%s</Attribute>\r\n', strtrim(v.ImageInformation.ImageTypeEdEs(i,:)) );
                fprintf( fid,'<Attribute Name="Pixel Spacing" Tag="0x00280030" Type="Double" ArraySize="2">%4.6E %4.6E</Attribute>\r\n', v.ImageInformation.PixelSpacing(i,1), v.ImageInformation.PixelSpacing(i,2) );
                fprintf( fid,'<Attribute Name="Echo Time" Tag="0x00180081" Type="Double">%4.6E</Attribute>\r\n', v.ImageInformation.EchoTime(i) );
                temp1 = sprintf( '%2.4E', v.ImageInformation.DynScanBeginTime(i) );
                temp1( end-2 ) = [];
                fprintf( fid,'<Attribute Name="Dyn Scan Begin Time" Tag="0x200510A0" Type="Float">%s</Attribute>\r\n', temp1 );
                fprintf( fid,'<Attribute Name="Trigger Time" Tag="0x00181060" Type="Double">%4.6E</Attribute>\r\n', v.ImageInformation.TriggerTime(i) );
                temp1 = sprintf( '%2.4E', v.ImageInformation.DiffusionBFactor(i) );
                temp1( end-2 ) = [];
                fprintf( fid,'<Attribute Name="Diffusion B Factor" Tag="0x20011003" Type="Float">%s</Attribute>\r\n', temp1 );
                fprintf( fid,'<Attribute Name="No Averages" Tag="0x00180083" Type="Double">%4.6E</Attribute>\r\n', v.ImageInformation.NoAverages(i) );
                fprintf( fid,'<Attribute Name="Image Flip Angle" Tag="0x00181314" Type="Double">%4.6E</Attribute>\r\n', v.ImageInformation.ImageFlipAngle(i) );
                fprintf( fid,'<Attribute Name="Cardiac Frequency" Tag="0x00181088" Type="Int32">%d</Attribute>\r\n', v.ImageInformation.CardiacFrequency(i) );
                fprintf( fid,'<Attribute Name="Min RR Interval" Tag="0x00181081" Type="Int32">%d</Attribute>\r\n', v.ImageInformation.MinRRInterval(i) );
                fprintf( fid,'<Attribute Name="Max RR Interval" Tag="0x00181082" Type="Int32">%d</Attribute>\r\n', v.ImageInformation.MaxRRInterval(i) );
                fprintf( fid,'<Attribute Name="TURBO Factor" Tag="0x00180091" Type="Int32">%d</Attribute>\r\n', v.ImageInformation.TURBOFactor(i) );
                fprintf( fid,'<Attribute Name="Inversion Delay" Tag="0x00180082" Type="Double">%4.6E</Attribute>\r\n', v.ImageInformation.InversionDelay(i) );
                fprintf( fid,'<Attribute Name="Contrast Type" Tag="0x00089209" Type="String">%s</Attribute>\r\n', strtrim(v.ImageInformation.ContrastType(i,:)) );
                fprintf( fid,'<Attribute Name="Diffusion Anisotropy Type" Tag="0x00189147" Type="String">%s</Attribute>\r\n', strtrim(v.ImageInformation.DiffusionAnisotropyType(i,:)) );
                temp1 = sprintf( '%2.4E', v.ImageInformation.DiffusionAP(i) );
                temp1( end-2 ) = [];
                fprintf( fid,'<Attribute Name="Diffusion AP" Tag="0x200510B1" Type="Float">%s</Attribute>\r\n', temp1 );
                temp1 = sprintf( '%2.4E', v.ImageInformation.DiffusionFH(i) );
                temp1( end-2 ) = [];
                fprintf( fid,'<Attribute Name="Diffusion FH" Tag="0x200510B2" Type="Float">%s</Attribute>\r\n', temp1 );
                temp1 = sprintf( '%2.4E', v.ImageInformation.DiffusionRL(i) );
                temp1( end-2 ) = [];
                fprintf( fid,'<Attribute Name="Diffusion RL" Tag="0x200510B0" Type="Float">%s</Attribute>\r\n', temp1 );
                fprintf( fid,'<Attribute Name="Angulation AP" Type="Double" Calc="AngulationAP">%4.6E</Attribute>\r\n', v.ImageInformation.AngulationAP(i) );
                fprintf( fid,'<Attribute Name="Angulation FH" Type="Double" Calc="AngulationFH">%4.6E</Attribute>\r\n', v.ImageInformation.AngulationFH(i) );
                fprintf( fid,'<Attribute Name="Angulation RL" Type="Double" Calc="AngulationRL">%4.6E</Attribute>\r\n', v.ImageInformation.AngulationRL(i) );
                fprintf( fid,'<Attribute Name="Offcenter AP" Type="Double" Calc="OffcenterAP">%4.6E</Attribute>\r\n', v.ImageInformation.OffcenterAP(i) );
                fprintf( fid,'<Attribute Name="Offcenter FH" Type="Double" Calc="OffcenterFH">%4.6E</Attribute>\r\n', v.ImageInformation.OffcenterFH(i) );
                fprintf( fid,'<Attribute Name="Offcenter RL" Type="Double" Calc="OffcenterRL">%4.6E</Attribute>\r\n', v.ImageInformation.OffcenterRL(i) );
                fprintf( fid,'<Attribute Name="Slice Orientation" Type="Enumeration" Calc="SliceOrient" EnumType="Slice_Orientation">%s</Attribute>\r\n', strtrim(v.ImageInformation.SliceOrientation(i,:)) );
                fprintf( fid,'</Image_Info>\r\n' );
                
            end
            
            fprintf( fid, '</Image_Array>\r\n' );
            fprintf( fid, '</PRIDE_V5>\r\n' );
            
            MR.Data = Helper.UnconvertCell( MR.Data );
            MR.Parameter.Recon.AutoUpdateInfoPars = auto_update_status;
            
            fclose(fid);
        end
        function WriteCpx( MR, Filename )
            
            % Data: 1:resolution y, 2:resolution x 3:Stack,
            %       4:Slice, 5:Coil, 6:Heart phase, 7:Echo,
            %       8:Dynamics, 9:Segments, 10:Segments2
            
            
            %    h1(2);                  % Stack
            %    h1(3);                  % Slice
            %    h2(2);                  % Coil
            %    h1(6);                  % Heart phase
            %    h1(5);                  % Echo
            %    h1(7);                  % Dynamics
            %    h1(8);                  % Segments
            %    h1(10);                 % Data offset
            %    h1(11);                 % Resolution x
            %    h1(12);                % Resolution y
            %    h1(14);                % Compression
            %    h2(111);               % Flip
            %    factor(1);             % Scaling factor 1
            %    factor(2);             % Scaling factor 2
            %    h1(1);                 % mix
            %    h1(4);                 % Prep Dir
            %    h1(15);                % Sequence Number
            %    h2(1);                 % Segment2
            %    h2(3);                 % Syncho number
            
            if nargin == 1
                Filename = 'recon.cpx';
            end
            
            MR.Parameter.UpdateImageInfo = 0;
            MR.DataClass.Convert2Cell;
            MR.Parameter.UpdateImageInfo = 1;
            
            switch class(MR.Data{1})
                case('double')
                    compression = 1;
                case('single')
                    compression = 1;
                otherwise
                    MR.Data = single(MR.Data);
                    compression = 1;
            end
            
            fid = fopen(Filename,'w');
            
            h1 = zeros(15,1);
            h2 = zeros(15,1);
            skip = zeros(384,1);
            
            for m = 1:size(MR.Data,2)
                siz = ones(1,12);
                siz(1:ndims(MR.Data{1,m})) = size(MR.Data{1,m} );
                
                res_x = siz(1);
                res_y = siz(2);
                data_slice = zeros(ceil(res_x*res_y*8/512)*512/8*2,1);
                
                siz = siz(3:end);
                total_nr_images = size(MR.Data{1,m}(:,:,:),3 );
                
                for i = 1:total_nr_images
                    sub = Helper.ind2sub(siz,i);
                    % [kz, chan, dyn, card, echo, loca, mix, extr1,extr2, aver]
                    
                    if size(MR.Data,2) == 1
                        real_mix = sub(7);
                    else
                        real_mix = m;
                    end
                    
                    slice = sub(1);
                    chan  = sub(2);
                    dyn   = sub(3);
                    card  = sub(4);
                    echo  = sub(5);
                    loca  = siz(6)*(real_mix-1) + sub(6);
                    extr1 = sub(8);
                    extr2 = siz(9)*(sub(10)-1) + sub(9);
                    
                    h1(2) = loca ;
                    h1(3) = slice ;
                    h2(2) = chan ;
                    h1(6) = card ;
                    h1(5) = echo ;
                    h1(7) = dyn ;
                    h1(8) = extr1 ;
                    h2(1) = extr2 ;
                    
                    h1(11)= res_x;
                    h1(12)= res_y;
                    h1(14)= compression;
                    h1(9) = 1;
                    h1(13)= ceil(res_x*res_y*8/512);
                    h1(10)= ftell(fid)+512;
                    
                    fwrite(fid,h1,'long');
                    fwrite(fid,0,'float');
                    fwrite(fid,0,'float');
                    fwrite(fid,h2,'long');
                    fwrite(fid,skip,'char');
                    
                    data_slice1 = MR.Data{1,m}(:,:,i);
                    data_slice1 = reshape(data_slice1,res_x*res_y,1);
                    data_slice(1:2:2*length(data_slice1))=real(data_slice1);
                    data_slice(2:2:2*length(data_slice1))=imag(data_slice1);
                    fwrite(fid,data_slice,'float');
                    
                end
            end
            
            
            h1(9) = 0;
            fwrite(fid,h1,'long');
            fwrite(fid,0,'float');
            fwrite(fid,0,'float');
            fwrite(fid,h2,'long');
            fwrite(fid,skip,'char');
            
            fclose( fid );
        end
        function ExportLabels( MR, Filename )
            % ExportLabels: Exports the imaging labels, which are stored in the .lab file to a
            % textfile
            %
            % Syntax:     r.ExportLabels( filename );
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('MRparameterDoc.Labels.Index')">Parameter.Labels.Index</a>:
            %             The labels to be exported.
            %
            % Location:   k-space, image-space.
            %
            % Formats:    Raw | ExportedRaw | Bruker
            %
            % Description/Algorithm: This function exports the imaging labels, which were sent from
            %             the spectrometer to the reconstructor, to a textfile. Every row in the
            %             file corresponds to an measured k-space line in the order they have been
            %             acquired. Exporting the labels is useful to get an overview over the
            %             sampling order for example.
            %
            % Notes:    - The format of the textfile is virtually identical to the one of the .list
            %             file of the exported raw data format.
            
            if ~MR.Parameter.ReconFlags.isreadparameter
                error( 'Please read the parameter file first' );
            end
            
            if nargin == 1
                Filename = 'labels.list';
            end
            Writer.export_labels( MR.Parameter.Labels.Index, Filename );
        end
        function WriteDICOM( MR, OutputDirectory )
            % Check if the data is empty
            if isempty( MR.Data )
                error( 'Error in WriteDICOM: The data matrix is empty. Please reconstruct the data first');
            end
            
            if strcmpi(MR.Parameter.DataFormat,'raw') && ~isempty( MR.Parameter.Sin )
                warning('The dataset was not measured with the ReconFrame patch. Some DICOM tags will be missing');
            elseif strcmpi(MR.Parameter.DataFormat,'raw') && isempty( MR.Parameter.Par40 )
                warning('The dataset was measured with the old ReconFrame patch. Some DICOM tags will be missing');
            elseif strcmpi(MR.Parameter.DataFormat,'rec')
                warning('Some DICOM tags will be missing when converting a REC file to DICOM');
            else
                if( ~MR.Parameter.IsParameter('RFR_EXAM_OBJECT_OID') && ~MR.Parameter.IsParameter('RFR_SERIES_objOid'))
                    warning('DICOM export on release 3 data is not fully supported. Some DICOM tags will be missing');
                end
                if( MR.Parameter.IsParameter('RFR_SERIES_objOid'))
                    warning('DICOM export on release 11 is an experimental feature');
                end
            end
            
            if nargin == 1
                OutputDirectory = cd;
            end
            
            slashind = [strfind( MR.Parameter.Filename.Data, '\' ), strfind( MR.Parameter.Filename.Data, '/' ) ];
            dotind = strfind( MR.Parameter.Filename.Data, '.' );
            if isempty( slashind )
                slashind = 1;
            end
            if isempty( dotind )
                dotind = length( MR.Parameter.Filename.Data );
            end
            dataset_name = MR.Parameter.Filename.Data( slashind(end)+1:dotind(end)-1 );
            
            par = MR.WritePar2String;
            rec = MR.GetRecImages;
            Parfile = cell(1, length(par));
            for i = 1:length(par)
                Parfile{i} = MR.Parameter.parread_from_string(par{i});
            end
            D = DICOMExporter;
            
            if( ~isfolder( [OutputDirectory, filesep, dataset_name] ))
                try
                    mkdir(OutputDirectory, dataset_name);
                catch
                    error( ['cannot create the directory: ', OutDir, filesep, 'DICOM'] );
                end
            end
            
            OutputDirectory = [OutputDirectory, filesep, dataset_name];
            if( isempty(MR.Parameter.Par40) )
                D.Export(OutputDirectory, rec, Parfile);
            else
                D.Export(OutputDirectory, rec, Parfile, MR.Parameter.Par40);
            end
        end
        function WriteExamCard(MR, Filename)
            if isempty( MR.Parameter.Par40 )
                error('The dataset was measured with the old ReconFrame patch. Cannot export an ExamCard.');
            end
            
            if nargin == 1
                Filename = [];
            end
            
            MR.Parameter.Par40.WriteExamCard(Filename);
        end
        
        % ---------------------------------------------------------------%
        % Data Sorting
        % ---------------------------------------------------------------%
        function SortData( MR, varargin )
            % SortData: Sorts the acquired profiles by their image attributes and generates proper
            % k-space data.
            %
            % Syntax:     r.SortData;
            %
            % Parameters used:
            %
            %           All formats:
            %           - <a href="matlab:helpwin('ReconParsDoc.ImmediateAveraging')">Parameter.Recon.ImmediateAveraging</a>:
            %             When enabled the averages are added immediately during SortData. If you
            %             want to reconstruct the individual averages seperately then disable this
            %             parameter.
            %            - <a href="matlab:helpwin('MRparameterDoc.Labels')">Parameter.Labels</a>:
            %             The data labels which hold the image attributes for every profile.
            %           - <a href="matlab:helpwin('MRparameterDoc.LabelLookupTable')">Parameter.LabelLookupTable</a>:
            %             Links the profiles which are currently stored in the Data array with the
            %             data lables.
            %           - <a href="matlab:helpwin('EncodingParsDoc.KyRange')">Parameter.Encoding.KyRange</a>,  <a href="matlab:helpwin('EncodingParsDoc.KzRange')">Parameter.Encoding.KzRange</a>:
            %             The sampled k-space ranges in phase encoding directions used to place the
            %             data correctly in k-space (e.g. in halfscan acquisitions).
            %           - <a href="matlab:helpwin('EncodingParsDoc.KyOversampling')">Parameter.Encoding.KyOversampling</a>,  <a href="matlab:helpwin('EncodingParsDoc.KyOversampling')">Parameter.Encoding.KyOversampling</a>:
            %             The oversampling factors in phase encoding directions used to place the
            %             data correctly in k-space.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.issorted')">Parameter.ReconFlags.issorted</a>
            %
            % Location:   After ReadData
            %
            % Formats:    Raw | ExportedRaw | ExportedCpx
            %
            % Description/Algorithm: After reading the data into Matlab with ReadData, the acquired
            %             profiles are stored in temporal order (in the order they have been
            %             measured). The SortData function sorts the profiles by their image
            %             attributes phase encoding number, cardiac phase, dynamic scan number etc.
            %             After the sorting the data the data array holds proper k-space data which
            %             can be fourier transformed to generate images. Please note that averages
            %             are summed up immediately during sorting if immediate averaging is enabled
            %             (see above).
            %
            %             After SortData the data array has a predifined form with 12 dimensions,
            %             where every array dimension corresponds to an image attribute:
            %
            %             r.Data = [kx, ky, kz, coils, dynamics, cadiac phases, echoes, locations,
            %                       mixes, extra1, extra2, averages ]
            
            if ~MR.Parameter.ReconFlags.isread
                error( 'Please read the data first' );
            end
            if MR.Parameter.ReconFlags.issorted
                error( 'The data is already sorted' );
            end
            if isfield(MR.Parameter.Labels, 'NumberOfEncodingDimensions') && any(MR.Parameter.Labels.NumberOfEncodingDimensions > 1) && isempty(MR.Parameter.Encoding.KyRange)
                error('r.Parameter.Encoding.KyRange is empty');
            end
            
            defaultAverageIdenticalProfs = [];
            validAverageIdenticalProfs = @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= 2);
            
            p = inputParser;
            addParameter(p,'AverageIdenticalProfs',defaultAverageIdenticalProfs,validAverageIdenticalProfs);
            parse(p,varargin{:});
            
            if isempty(p.Results.AverageIdenticalProfs)
                if any( strcmpi( MR.Parameter.Cardiac.Synchronization, {'retro', 'retrospective'} ))
                    switch MR.Parameter.Cardiac.RetroHoleInterpolation
                        case {'Nearest', 'nearest'}
                            AverageIdenticalProfs = 0;
                        case {'Average', 'average'}
                            AverageIdenticalProfs = 1;
                        case {'Linear', 'linear', 'Cubic', 'cubic'}
                            AverageIdenticalProfs = 2;
                        otherwise
                            AverageIdenticalProfs = 0;
                    end
                else
                    AverageIdenticalProfs = 1;
                end
            else
                AverageIdenticalProfs = p.Results.AverageIdenticalProfs;
            end
                    
            ImmediateAveraging    = strcmpi( MR.Parameter.Recon.ImmediateAveraging, 'yes' ) && strcmpi( MR.Parameter.Recon.Average, 'yes' );
            
            if ~isempty(MR.Data)
                MR.Parameter.UpdateImageInfo = 0;
                MR.Data = Helper.Convert2Cell( MR.Data );
                
                MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                
                NoZeroPad(1) = any( strcmpi( MR.Parameter.Scan.AcqMode, {'Radial', 'Spiral'}));
                % Allow zero padding for a stack of radial or spirals
                % (halfscan along z direction)
                NoZeroPad(2) = (strcmpi( MR.Parameter.Scan.AcqMode, 'Radial') && strcmpi( MR.Parameter.Scan.KooshBall, 'yes') );
                               
                % Turn off the zero filling for non-cartesian spectro (spiral)
                if MR.isSpectro && any( strcmpi( MR.Parameter.Scan.AcqMode, {'Radial', 'Spiral'}))
                    NoZeroPad(:) = 1;
                end                
                
                if isfield( MR.Parameter.Labels, 'Samples' )
                    samples = MR.Parameter.Labels.Samples(1,:);
                    ovs = cellfun(@(x,y,z)[max([x,1]),max([y,1]),max([z,1])], MR.Parameter.Encoding.WorkEncoding.KxOversampling, MR.Parameter.Encoding.WorkEncoding.KyOversampling, MR.Parameter.Encoding.WorkEncoding.KzOversampling, 'UniformOutput', 0);
                    samples = cellfun( @(x) round(x.*samples), ovs, 'UniformOutput', 0);
                else
                    samples = cellfun(@(x)[], MR.Parameter.Encoding.WorkEncoding.KxOversampling, 'UniformOutput', 0);
                end
                FixedKyRange = ~isfield( MR.Parameter.Labels, 'KxRange' );
                
                [MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable] = Helper.check_cell_sizes( MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable, MR.Data );
                
                % Check if the ky and kz labels are inside the KyRange,
                % KzRange
                ky_outside_range = ~NoZeroPad(1) && (any(any(cellfun(@(x,y) max([min(MR.Parameter.Labels.Index.ky(x)) < min(y), 0]) > 0, MR.Parameter.LabelLookupTable, MR.Parameter.Encoding.WorkEncoding.KyRange ))) || ...
                    any(any(cellfun(@(x,y) max([max(MR.Parameter.Labels.Index.ky(x)) > max(y), 0]) > 0, MR.Parameter.LabelLookupTable, MR.Parameter.Encoding.WorkEncoding.KyRange ))));
                if( ky_outside_range )
                    error( 'The ky-labels (Parameter.Labels.Index.ky) are outside the ky-range (r.Parameter.Encoding.KyRange). Please adapt your parameters');
                end
                kz_outside_range = ~NoZeroPad(2) && (any(any(cellfun(@(x,y) max([min(MR.Parameter.Labels.Index.kz(x)) < min(y), 0]) > 0, MR.Parameter.LabelLookupTable, MR.Parameter.Encoding.WorkEncoding.KzRange ))) || ...
                    any(any(cellfun(@(x,y) max([max(MR.Parameter.Labels.Index.kz(x)) > max(y), 0]) > 0, MR.Parameter.LabelLookupTable, MR.Parameter.Encoding.WorkEncoding.KzRange ))));
                if( kz_outside_range )
                    error( 'The kz-labels (Parameter.Labels.Index.kz) are outside the kz-range (r.Parameter.Encoding.KzRange). Please adapt your parameter');
                end
                
                if( ismac || strcmpi(MR.Parameter.DataFormat, 'ExportedRaw') || strcmpi(MR.Parameter.Recon.UseMatlabInternals, 'yes') )
                    if AverageIdenticalProfs > 1
                        AverageIdenticalProfs = 1;
                    end
                    [MR.Data, MR.Parameter.LabelLookupTable] = cellfun( @(x,y,z,u,v,w)Sorter.sort_exported_raw( x, MR.Parameter.Labels, y, z, u, v,w, NoZeroPad, AverageIdenticalProfs, ImmediateAveraging ), ...
                        MR.Data, MR.Parameter.LabelLookupTable, ...
                        MR.Parameter.Encoding.WorkEncoding.KyRange, MR.Parameter.Encoding.WorkEncoding.KzRange, ...
                        MR.Parameter.Encoding.WorkEncoding.KyOversampling, MR.Parameter.Encoding.WorkEncoding.KzOversampling,'UniformOutput', false );
                else                    
                    [MR.Data, MR.Parameter.LabelLookupTable] = cellfun( @(x,y,z,u,v,w,s)sort_data( x, y, MR.Parameter.Labels.Index, z, u, v, w, s, NoZeroPad, AverageIdenticalProfs, ImmediateAveraging,  FixedKyRange), ...
                        MR.Data, MR.Parameter.LabelLookupTable, ...
                        MR.Parameter.Encoding.WorkEncoding.KyRange, MR.Parameter.Encoding.WorkEncoding.KzRange, ...
                        MR.Parameter.Encoding.WorkEncoding.KyOversampling, MR.Parameter.Encoding.WorkEncoding.KzOversampling, samples, 'UniformOutput', false );
                end
                MR.Parameter.UpdateImageInfo = 1;
                MR.Data = Helper.UnconvertCell( MR.Data );
                MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                
            end
            MR.Parameter.ReconFlags.issorted = 1;
            % spectro begin ----------------------------
            if ImmediateAveraging
                MR.Parameter.ReconFlags.isaveraged = 1;
            end
            % spectro end ------------------------------
            
        end
        
        % ---------------------------------------------------------------%
        % Corrections
        % ---------------------------------------------------------------%
        function BasicCorrections(MR)
            MR.NonLinearityCorrection;
            MR.RandomPhaseCorrection;
            MR.RemoveOversampling;
            MR.PDACorrection;
            MR.DcOffsetCorrection;
            MR.MeasPhaseCorrection;
        end
        function NonLinearityCorrection( MR )
            % NonLinearityCorrection: Corrects for receiver channel non-linearities
            %
            % Syntax:     r.NonLinearityCorrection;
            %            
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.israndphasecorr')">Parameter.ReconFlags.isnonlincorr</a>
            %
            % Location:   k-space
            %
            % Formats:    Raw
            %  
            % Description/Algorithm: During the preparation phase the scanner 
            % sends channel non-linearity information to the reconstructor.
            % This function corrects for these non-linearities.
            
            if strcmpi( MR.Parameter.Recon.NonLinearityCorrection, 'yes' ) && ...
                    strcmpi( MR.Parameter.Recon.ArrayCompression, 'no' ) && ...
                    strcmpi( MR.Parameter.DataFormat, 'raw' ) && ...
                    isfield(MR.Parameter.Labels, 'NonLinLabels')                          
                
                if ~MR.Parameter.ReconFlags.isread
                    error( 'Error in Non-Linearity Correction: Please read the data first' );
                end
                if MR.Parameter.ReconFlags.ispartialfourier
                    error( 'Error in Non-Linearity Correction: The random phase correction has to be applied before the partial fourier reconstruction' );
                end
                if MR.Parameter.ReconFlags.isnonlincorr
                    error( 'Error in Non-Linearity Correction: The data is already random-phase corrected' );
                end
                if any( MR.Parameter.ReconFlags.isimspace)
                    error( 'Error in Non-Linearity Correction: The random phase correction has to be applied in k-space' );
                end
                
                if ~isempty(MR.Data)
                    
                    MR.Parameter.UpdateImageInfo = 0;
                    MR.DataClass.Convert2Cell;
                    
                    MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                                                           
                    MR.Data = cellfun( @(x,y)BasicCorrections.nonlin_corr( x, MR.Parameter.Labels, y, MR.Parameter.Filename.Data ), ...
                        MR.Data, MR.Parameter.LabelLookupTable, 'UniformOutput', false );                   
                    
                    [MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable] = Helper.check_cell_sizes( MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable, MR.Data );
                                                         
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    
                    MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                end
                MR.Parameter.ReconFlags.isnonlincorr = 1;
                
            end
        end
        function RandomPhaseCorrection( MR )
            % RandomPhaseCorrection: Corrects for the random phase which is added to the measured
            % profiles
            %
            % Syntax:     r.RandomPhaseCorrection;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.RandomPhaseCorrection')">Parameter.Recon.RandomPhaseCorrection</a>: {'yes'} | 'no'
            %             Enables/disables the random phase correction.
            %           - <a href="matlab:helpwin('MRparameterDoc.DataFormat')">Parameter.DataFormat</a>:
            %             Used to check if raw data is passed to the function
            %           - <a href="matlab:helpwin('MRparameterDoc.LabelLookupTable')">Parameter.LabelLookupTable</a>:
            %             Defines the profiles which are currently stored in the Data array.
            %           - <a href="matlab:helpwin('MRparameterDoc.Labels')">Parameter.Labels.Index.random_phase</a>:
            %             Holds the random phases which are subtracted from the profiles.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.israndphasecorr')">Parameter.ReconFlags.israndphasecorr</a>
            %
            % Location:   k-space
            %
            % Formats:    Raw
            %
            % Description/Algorithm: The receiver adds a random phase to every acquired profile to
            %             reduce artifacts. This random phase has to be subtracted in
            %             reconstruction since the phase is used for spatial encoding. The added
            %             phase is stored in the profile labels and the phase which is subtraced is
            %             given by:
            %
            %             Let ri be the random phase label of the i-th profile:
            %
            %                 ri = Parameter.Labels.Index.meas_phase(i);
            %
            %             then the phase subtracted to that profile is given by:
            %
            %                 phase = 2*pi / double( intmax('uint16') ) * ri;
            
            if strcmpi( MR.Parameter.Recon.RandomPhaseCorrection, 'yes' ) && ...
                    strcmpi( MR.Parameter.Recon.ArrayCompression, 'no' ) && ...
                    strcmpi( MR.Parameter.DataFormat, 'raw' )
                
                if ~MR.Parameter.ReconFlags.isread
                    error( 'Error in Random Phase Correction: Please read the data first' );
                end
                if MR.Parameter.ReconFlags.ispartialfourier
                    error( 'Error in Random Phase Correction: The random phase correction has to be applied before the partial fourier reconstruction' );
                end
                if MR.Parameter.ReconFlags.israndphasecorr
                    error( 'Error in Random Phase Correction: The data is already random-phase corrected' );
                end
                if any( MR.Parameter.ReconFlags.isimspace)
                    error( 'Error in Random Phase Correction: The random phase correction has to be applied in k-space' );
                end
                
                if ~isempty(MR.Data)
                    
                    MR.Parameter.UpdateImageInfo = 0;
                    MR.DataClass.Convert2Cell;
                    
                    MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                    
                    % Decide weather a FEAR or REAR correction should be
                    % performed
                    if( isfield(MR.Parameter.Labels, 'FEARFactor' ) && MR.Parameter.Labels.FEARFactor > 0 )
                        if ~MR.Parameter.ReconFlags.isoversampled(1)
                            error( 'Error in Random Phase Correction: The RandomPhaseCorrection must be performed before RemoveOversampling for release 5 data and higher' );
                        end
                        
                        center_k = MR.Parameter.Encoding.WorkEncoding.FEARCenterK;
                        radial_spiral =  any( strcmpi( MR.Parameter.Scan.AcqMode, {'radial', 'spiral'})) || strcmpi( MR.Parameter.Gridder.Preset, 'radial' );
                        if( radial_spiral )
                            center_k = cellfun( @(x)[], center_k, 'UniformOutput', 0);
                        end
                        
                        MR.Data = cellfun( @(x,y,z, u)BasicCorrections.fear_corr( x, MR.Parameter.Labels, y, z, u, radial_spiral ), ...
                            MR.Data, MR.Parameter.LabelLookupTable, MR.Parameter.Encoding.WorkEncoding.FEARFactor, center_k, 'UniformOutput', false );
                    else
                        MR.Data = cellfun( @(x,y)BasicCorrections.random_phase_corr( x, MR.Parameter.Labels, y ), ...
                            MR.Data, MR.Parameter.LabelLookupTable, 'UniformOutput', false );
                    end
                    
                    [MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable] = Helper.check_cell_sizes( MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable, MR.Data );
                                                         
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    
                    MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                end
                MR.Parameter.ReconFlags.israndphasecorr = 1;
                
            end
        end
        function MeasPhaseCorrection( MR )
            % MeasPhaseCorrection: Corrects for an offset phase which might have been added to
            % certain profiles
            %
            % Syntax:     r.DcOffsetCorrection;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.MeasPhaseCorrection')">Parameter.Recon.MeasPhaseCorrection</a>: {'yes'} | 'no'
            %             Enables/disables the measurement phase correction.
            %           - <a href="matlab:helpwin('MRparameterDoc.DataFormat')">Parameter.DataFormat</a>:
            %             Used to check if raw data is passed to the function
            %           - <a href="matlab:helpwin('MRparameterDoc.LabelLookupTable')">Parameter.LabelLookupTable</a>:
            %             Defines the profiles which are currently stored in the Data array.
            %           - <a href="matlab:helpwin('MRparameterDoc.Labels')">Parameter.Labels.Index.meas_phase</a>:
            %             Holds the measurement phases.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.isdcoffsetcorr')">Parameter.ReconFlags.isdcoffsetcorr</a>
            %
            % Location:   k-space
            %
            % Formats:    Raw
            %
            % Description/Algorithm: Different profiles can sometimes exhibit a phase shift, which
            %             originates from either the RF pulse phase or the acquistion phase
            %             specified in the AQ object. Different averages for example are usually
            %             acquired with a phase difference of pi between them. If we wouldn't
            %             correct for this phase then the signal would effectively cancel out when
            %             summing up the averages. MeasPhaseCorrection therefore corrects for this
            %             possible phase shift using the values in the labels, received by the
            %             scanner.
            %             Let mi be the measurement of the i-th profile:
            %
            %                 mi = Parameter.Labels.Index.meas_phase(i);
            %
            %             then the phase added to that profile is given by:
            %
            %                 phase = -mi * pi/2
            
            if strcmpi( MR.Parameter.Recon.MeasPhaseCorrection, 'yes' ) && ...
                    strcmpi( MR.Parameter.Recon.ArrayCompression, 'no' ) && ...
                    strcmpi( MR.Parameter.DataFormat, 'raw' )
                
                if ~MR.Parameter.ReconFlags.isread
                    error( 'Error in Measurement Phase Correction: Please read the data first' );
                end
                if any( MR.Parameter.ReconFlags.isimspace)
                    error( 'Error in Measurement Phase Correction: The random phase correction has to be applied in k-space' );
                end
                if MR.Parameter.ReconFlags.ismeasphasecorr
                    error( 'Error in Measurement Phase Correction: The data is already measurement-phase corrected' );
                end
                if MR.Parameter.ReconFlags.issorted && strcmpi(MR.Parameter.Recon.ImmediateAveraging, 'yes')
                    error( 'Error in Measurement Phase Correction: Please apply the measurement before averaging the data. Either perform it before SortData or switch immediate averaging off in Parameter.Recon.ImmediateAveraging' );
                end
                
                if( MR.isSpectro == 1 )
                    warning( 'Running the measurement phase correction on spectro data is not recommended. The spectra quality might suffer considerably' );
                end
                
                if ~isempty(MR.Data)
                    MR.Parameter.UpdateImageInfo = 0;
                    MR.DataClass.Convert2Cell;
                    
                    MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                    
                    [MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable] = Helper.check_cell_sizes( MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable, MR.Data );
                    
                    MR.Data = cellfun( @(x,y)BasicCorrections.meas_phase_corr( x, MR.Parameter.Labels, y ), ...
                        MR.Data, MR.Parameter.LabelLookupTable, 'UniformOutput', false );
                    
                    %                     cellfun( @(x,y)meas_phase_correction( x, y, MR.Parameter.Labels.Index ), ...
                    %                         MR.Data, MR.Parameter.LabelLookupTable, 'UniformOutput', false );
                    
                    
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    
                    MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                    
                end
                MR.Parameter.ReconFlags.ismeasphasecorr = 1;
            end
        end
        function DcOffsetCorrection( MR )
            % DcOffsetCorrection: Corrects for a signal offset that might occur in the acquired
            % samples.
            %
            % Syntax:     r.DcOffsetCorrection;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.DcOffsetCorrection')">Parameter.Recon.DcOffsetCorrection</a>: {'yes'} | 'no'
            %             Enables/disables the DC-offset correction.
            %           - <a href="matlab:helpwin('MRparameterDoc.DataFormat')">Parameter.DataFormat</a>:
            %             Used to check if raw data is passed to the function
            %           - <a href="matlab:helpwin('MRparameterDoc.LabelLookupTable')">Parameter.LabelLookupTable</a>:
            %             Defines the profiles which are currently stored in the Data array.
            %           - <a href="matlab:helpwin('EncodingPars.KxRange')">Parameter.Encoding.KxRange</a>:
            %             If no noise samples are present in the raw file, this parameter specifies
            %             where the outer k-space (noise like) samples are located.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.isdcoffsetcorr')">Parameter.ReconFlags.isdcoffsetcorr</a>
            %
            % Location:   k-space (before the partial fourier reconstruction)
            %
            % Formats:    Raw
            %
            % Description/Algorithm: Acquired raw samples can sometimes exhibit a signal offset
            %             originating from hardware inacuracies of the receiver. If not corrected
            %             this offset can cause a artifact which will appear as bright signal peak
            %             in the middle of the image. Since noise in the MR signal is normally
            %             distributed with a mean of zero, the offset can be calculated by
            %             calculating the mean of the noise samples. Therefore if the raw data
            %             contains noise samples (from a noise acquisition prior to the scan), these
            %             samples are read and the mean is calculated. If the data does not contain
            %             any noise, the offset is calculated from the noise like outer k-space
            %             positions. After calculating the offset it is subtracted from the data.
            %
            % Notes:    - The offset is calculated seperately for every channel.
            
            if strcmpi( MR.Parameter.Recon.DcOffsetCorrection, 'yes' ) && ...
                    strcmpi( MR.Parameter.Recon.ArrayCompression, 'no' ) && ...
                    strcmpi( MR.Parameter.DataFormat, 'raw' )
                
                if ~MR.Parameter.ReconFlags.isread
                    error( 'Error in DC Offset Correction: Please read the data first' );
                end
                if MR.Parameter.ReconFlags.ispartialfourier
                    error( 'Error in DC Offset Correction: The DC Offset correction has to be applied before the partial fourier reconstruction' );
                end
                if MR.Parameter.ReconFlags.isdcoffsetcorr
                    error( 'Error in DC Offset Correction: The data is already DC Offset corrected' );
                end
                if any( MR.Parameter.ReconFlags.isimspace)
                    error( 'Error in DC OffsetA Correction: The DC Offset correction has to be applied in k-space' );
                end
                
                if ~isempty(MR.Data)
                    if isempty( MR.MRImpl.MeanNoise )
                        MR.MRImpl.MeanNoise = BasicCorrections.get_mean_noise( MR.Parameter.Filename, ...
                            MR.Parameter.DataType, MR.Parameter.Labels );
                    end
                    
                    MR.Parameter.UpdateImageInfo = 0;
                    MR.DataClass.Convert2Cell;
                    
                    MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                    
                    [MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable] = Helper.check_cell_sizes( MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable, MR.Data );
                    
                    MR.Data = cellfun( @(x,y,z)BasicCorrections.dc_offset_corr( x, MR.MRImpl.MeanNoise, MR.Parameter.Labels, y,  z, MR.Parameter.ReconFlags.ispdacorr ), ...
                        MR.Data, MR.Parameter.Encoding.WorkEncoding.KxRange, MR.Parameter.LabelLookupTable, 'UniformOutput', false );                                      
                    
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    
                    MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                    
                end
                MR.Parameter.ReconFlags.isdcoffsetcorr = 1;
            end
        end
        function PDACorrection( MR )
            % PDACorrection: Corrects the profile dependent amplification (PDA)
            %
            % Syntax:     r.PDACorrection;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.PDACorrection')">Parameter.Recon.PDACorrection</a>: {'yes'} | 'no'
            %             Enables/disables the PDA correction.
            %           - <a href="matlab:helpwin('MRparameterDoc.DataFormat')">Parameter.DataFormat</a>:
            %             Used to check if raw data is passed to the function
            %           - <a href="matlab:helpwin('MRparameterDoc.LabelLookupTable')">Parameter.LabelLookupTable</a>:
            %             Defines the profiles which are currently stored in the Data array.
            %           - <a href="matlab:helpwin('MRparameterDoc.Labels')">MR.Parameter.Labels.PDAFactors</a>:
            %             Holds the PDA correction factors.
            %           - <a href="matlab:helpwin('MRparameterDoc.Labels')">MR.Parameter.Labels.Index.pda_index</a>:
            %             Specifies which correction factor should be taken for a certain profile.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.ispdacorr')">Parameter.ReconFlags.ispdacorr</a>
            %
            % Location:   k-space
            %
            % Formats:    Raw
            %
            % Description/Algorithm: In some cases the receiver weights certain k-space profiles
            %             differently than others. In 3D imaging for example the outer k-space lines
            %             are amplified compared to the central ones, to ensure a sufficiant signal
            %             range. Of course the amplifacation has to be corrected in
            %             reconstruction, if not the edges in the image would be enhanced. The
            %             amplification factors are provided by the scanner and this function simply
            %             divides the profiles by them.
            
            if strcmpi( MR.Parameter.Recon.PDACorrection, 'yes' ) && ...
                    strcmpi( MR.Parameter.Recon.ArrayCompression, 'no' ) && ...
                    strcmpi( MR.Parameter.DataFormat, 'raw' ) && ...
                    isfield( MR.Parameter.Labels, 'PDAFactors' ) && ...
                    ~all( MR.Parameter.Labels.PDAFactors == 0 )
                
                if ~MR.Parameter.ReconFlags.isread
                    error( 'Error in PDA Correction: Please read the data first' );
                end
                if MR.Parameter.ReconFlags.ispartialfourier
                    error( 'Error in PDA Correction: The PDA correction has to be applied before the partial fourier reconstruction' );
                end
                if MR.Parameter.ReconFlags.ispdacorr
                    error( 'Error in PDA Correction: The data is already PDA corrected' );
                end
                if any( MR.Parameter.ReconFlags.isimspace)
                    error( 'Error in PDA Correction: The PDA correction has to be applied in k-space' );
                end
                
                if ~isempty(MR.Data)
                    MR.Parameter.UpdateImageInfo = 0;
                    MR.DataClass.Convert2Cell;
                    
                    MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                    
                    [MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable] = Helper.check_cell_sizes( MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable, MR.Data );
                    
                    MR.Data = cellfun( @(x,y)BasicCorrections.pda_corr( x, MR.Parameter.Labels, y ), ...
                        MR.Data, MR.Parameter.LabelLookupTable, 'UniformOutput', false );                                       
                    
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    
                    MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                    
                end
                MR.Parameter.ReconFlags.ispdacorr = 1;
            end
        end
        function EPIPhaseCorrection( MR )
            % EPIPhaseCorrection: Performes the EPI correction to correct for the FOV/2 ghost
            % originating from eddy current effects.
            %
            % Syntax:     r.EPIPhaseCorrection;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.EPIPhaseCorrection')">Parameter.Recon.EPIPhaseCorrection</a>: {'yes'} | 'no'
            %             Enables/disables the EPI correction.
            %           - <a href="matlab:helpwin('MRparameterDoc.DataFormat')">Parameter.DataFormat</a>:
            %             Used to check if the right data format is passed to the function
            %           - <a href="matlab:helpwin('MRparameterDoc.EPICorrData')">Parameter.EPICorrData</a>:
            %             The linear/nonlinear phase used in the EPI correction. If a linear
            %             correction is used (default) then the phase is described by a slope and
            %             offset value. For the nonlinear phase correction the full correction
            %             profile is given. If this parameter is empty these values are calculated
            %             in the function itself. The parameter can also be filled with user defined
            %             slopes and offsets prior to the correction. Please see the example below
            %             on how to do that.
            %           - <a href="matlab:helpwin('ReconFlagParsDoc.iszerofilled')">Parameter.ReconFlags.iszerofilled</a>, <a href="matlab:helpwin('ReconFlagParsDoc.isoversampled')">Parameter.ReconFlags.isoversampled</a> and <a href="matlab:helpwin('ReconFlagParsDoc.isgridded')">Parameter.ReconFlags.isgridded</a>:
            %             Used to check in what reconstruction state the data currently is. The
            %             correction profiles have to be brought into the same state to match the
            %             imaging data.
            %           - <a href="matlab:helpwin('ReconParsDoc.Recon.EPICorrectionMethod')">Parameter.Recon.EPICorrectionMethod</a>: {linear} | nonlin
            %             Enables/disables the non-linear EPI correction (see below).
            %           - <a href="matlab:helpwin('ReconParsDoc.Recon.EPI2DCorr')">Parameter.Recon.EPI2DCorr</a>: yes | {no}
            %             Enables/disables the image based EPI correction (see below).
            %           - <a href="matlab:helpwin('EncodingPars')">Encoding Parameter</a>:
            %             Used to read and format the correction profiles approprately.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.isdepicorr')">Parameter.ReconFlags.isdepicorr</a>
            %
            % Location:   Between FFT in readout and phase encoding direction.
            %
            % Formats:    Raw | ExportedRaw | Bruker
            %
            % Description/Algorithm: In EPI imaging multiple k-space lines are acquired after a
            %             single excitation of the imaging volume. This is achieved by inverting the
            %             readout gradient after every acquired profile and thus sampling the
            %             k-space in a forward-backward fashion. Due to eddy current effects
            %             however, the positiv readout gradients (forward sampling) is not exactly
            %             the same as the negativ gradient (backward sampling), which results in a
            %             temporal shift between the forward and backward sampled profiles. After
            %             fourier transformation this shift results in a ghost like artifact located
            %             at half of the FOV.
            %
            %             To correct for this effect EPI correction profiles are acquired prior to
            %             every scan. The correction profiles consist of a full EPI shot, but
            %             acquired without phase encoding blips, therefore sampling the central
            %             k-space line multiple times. Since these profiles should be identical
            %             (apart from relaxation effects), the temporal shift, caused by eddy current
            %             effects of the readout gradients, can be calculated. A temporal shift in
            %             k-space corresponds to a linear phase in image space. Therefore the
            %             correction profiles are first transfered into image space and the phase
            %             difference between two subsequent profiles is calculated. The phase
            %             difference is then linearly fitted and the linear phase is added to every
            %             odd (or even) profile, thus realigning them in k-space.
            %
            %             Optionally a nonlinear correction can be enabled via <a href="matlab:helpwin('ReconParsDoc.Recon.EPICorrectionMethod')">Parameter.Recon.EPICorrectionMethod</a>.
            %             In the nonlinear correction the calculated linear phase is subtracted from
            %             the original phase difference. The remaining part is then smoothed and
            %             added to the linear phase.
            %
            %             Optionally an image based correction can be enabled via <a href="matlab:helpwin('ReconParsDoc.Recon.EPI2DCorr')">Parameter.Recon.EPI2DCorr</a>.
            %             This correction tries to improve the images further by minimizing the
            %             amount of ghosts in the image. Thereby the offsets found by the regular
            %             EPI crrection are slighly modified and the amount of ghosts in the
            %             resulting images is quantified, to find the optimal correction value.
            %
            % Notes:    - Every coil has its own EPI correction parameters.
            %           - The x-axis for the linear fit is given by:
            %                  x = -floor(n/2):ceil(n/2)-1;
            %             where n is the number of samples in the correction profiles.
            %           - The linear fit is weighted by the magnitude of the correction profiles
            %             (points with a high magnitude value are weighted higher in the fit)
            %           - The EPI correction profiles are labeled as data type 3 and can be read by
            %             setting: r.Parameter.Parameter2Read.typ = 3, followed by r.ReadData.
            %
            % Examples:   The user has the possibility to specify own correction parameters and pass
            %             them to the function. To do that <a href="matlab:helpwin('MRparameterDoc.EPICorrData')">Parameter.EPICorrData</a> must be filled:
            %
            %             Every coil has its own slope and offset. To identify the coil we must
            %             label them in the correction parameter:
            %
            %                  r.Parameter.EPICorrData.coil = r.Parameter.Parameter2Read.chan;
            %
            %             Assign the same slope for all the coils:
            %
            %                  r.Parameter.EPICorrData.slope = 0.1;
            %
            %             Assign a different offset for every coil (4 coils in this example):
            %
            %                  r.Parameter.EPICorrData.offset = [0, 0.2, 0.4, 0.6];
            %
            %             Call the EPI correction:
            %
            %                  r.EPIPhaseCorrection;
            
            if strcmpi( MR.Parameter.Recon.EPIPhaseCorrection, 'yes' )
                
                if ~isempty(MR.Data)
                    if ~MR.Parameter.ReconFlags.issorted
                        error( 'Please sort the data first' );
                    end
                    if ~MR.Parameter.ReconFlags.isimspace(1)
                        error( 'Please fourier transform the data along the measurement direction first' );
                    end
                    
                    if any( strcmpi( MR.Parameter.DataFormat, {'Raw', 'Bruker'} ) ) && ...
                            any( strcmpi( MR.Parameter.Scan.FastImgMode, {'EPI', 'TFEEPI', 'TFEPI' } ))
                        if MR.Parameter.ReconFlags.isdepicorr
                            error( 'The data is already EPI-phase corrected' );
                        end
                        
                        % Check if the Epi correction parameters should be
                        % recalculated.
                        % They are calculated again if the user didn't set them
                        % manually.
                        recalculate_EPI_corr = ~(MR.Parameter.EPICorrData.coil_changed || ...
                            MR.Parameter.EPICorrData.loca_changed || ...
                            MR.Parameter.EPICorrData.offset_changed || ...
                            MR.Parameter.EPICorrData.slope_changed || ...
                            MR.Parameter.EPICorrData.corr_factors_changed);
                        
                        if recalculate_EPI_corr
                            if isempty( MR.MRImpl.MeanNoise)
                                MR.MRImpl.MeanNoise = BasicCorrections.get_mean_noise( MR.Parameter.Filename, MR.Parameter.DataType, MR.Parameter.Labels );
                            end
                            
                            % Calculate the position off the isocenter. The linear phase has its zero crossing
                            % there.
                            isocenterMPS = MR.Transform( [0,0,0], 'xyz', 'ijk', 1);
                            
                            ACMatrix = [];
                            if strcmpi( MR.Parameter.Recon.ArrayCompression, 'yes' )
                                ACMatrix = MR.Parameter.Recon.ACMatrix;
                            end
                            
                            epi_corr_data = EPI.get_epi_corr_data( MR.Parameter.Filename, MR.Parameter.DataType, ...
                                MR.Parameter.Labels, MR.Parameter.Gridder.WorkingPars, MR.Parameter.Encoding, ...
                                MR.Parameter.LabelLookupTable, MR.Data, ...
                                MR.Parameter.ReconFlags.isgridded, MR.Parameter.ReconFlags.iszerofilled(1), ...
                                ~MR.Parameter.ReconFlags.isoversampled(1), strcmpi( MR.Parameter.Recon.EPI2DCorr, 'yes' ), ...
                                strcmpi( MR.Parameter.Recon.EPICorrectionMethod, 'Linear' ), isocenterMPS, ACMatrix, ...
                                strcmpi( MR.Parameter.Recon.EPICorrPerLocation, 'yes' ));
                            
                            MR.Parameter.SetEpiCorrData(epi_corr_data);
                        else
                            warning( 'The EPI correction parameters are not recalculated because they were manually set. If you want them to be recalculated call: r.Parameter.EpiCorrData.Reset');
                        end
                        
                        MR.Parameter.UpdateImageInfo = 0;
                        MR.DataClass.Convert2Cell;
                        
                        if isempty( MR.Parameter.LabelLookupTable )
                            MR.Parameter.LabelLookupTable = cellfun( @(x)[], MR.Data, 'UniformOutput', 0);
                        else
                            MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                        end
                        
                        if ~isempty( MR.Parameter.EPICorrData )
                            MR.Data = cellfun( @(x,y)EPI.epi_corr( x, MR.Parameter.EPICorrData, ...
                                MR.Parameter.Labels, y, strcmpi( MR.Parameter.Recon.EPICorrectionMethod, 'Linear' ), ...
                                strcmpi( MR.Parameter.Recon.EPICorrPerLocation, 'yes' ) ), ...
                                MR.Data, MR.Parameter.LabelLookupTable , 'UniformOutput', false );
                        end
                        
                        if recalculate_EPI_corr
                            MR.Parameter.EPICorrData.ResetFlags;
                        end
                        
                        MR.Parameter.UpdateImageInfo = 1;
                        MR.Data = Helper.UnconvertCell( MR.Data );
                        
                        MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                        
                    elseif strcmpi( MR.Parameter.DataFormat, 'ExportedRaw' )
                        if any( strcmpi( MR.Parameter.Scan.FastImgMode, {'EPI', 'TFEEPI' } ))
                            
                            phx_profiles = (MR.Parameter.Labels.Index.typ == 3);
                            
                            sign = double( MR.Parameter.Labels.Index.sign( phx_profiles ) );
                            sign( sign == -1 ) = 0;                          
                            MR.Parameter.Labels.Index.ky( phx_profiles ) = sign;
                            
                            % Read and Sort the EPI correction data and perform basic
                            % corrections
                            Parameter2Read = Parameter2ReadPars;
                            Parameter2Read.SetMax;
                            Parameter2Read.typ  = 3;
                            Parameter2Read.mix  = 0;
                            Parameter2Read.echo = 0;                           
            
                            [corr_data, corr_data_ind] = Reader.read_raw( MR.Parameter, Parameter2Read, ...
                            Parameter2Read.typ, Parameter2Read.mix, Parameter2Read.echo, [0,1], 0, 0, [], [] );
                            no_zero_pad = [true, true];
                            corr_data = Sorter.sort_exported_raw( corr_data, MR.Parameter.Labels, corr_data_ind, [0,1], [0,1], 1, 1, no_zero_pad, 1, false);
                                                   
                            corr_data = corr_data(end:-1:1,:,:,:,:);
                            coils = MR.Parameter.Parameter2Read.chan;
                            data = MR.Data;
                            data_ind = MR.Parameter.LabelLookupTable;
                            for c = 1:size(MR.Data,4)
                                for si = [-1,1]
                                    sign_ind = zeros( size(data_ind) );
                                    sign_ind(data_ind~= 0) = double( MR.Parameter.Labels.Index.sign( data_ind(data_ind~= 0) )) == si & ...
                                        double( MR.Parameter.Labels.Index.chan( data_ind(data_ind~= 0) )) == coils(c);
                                    data(:,sign_ind ) = bsxfun( @times, data(:,sign_ind ) ,corr_data(:,(si==1)+1,1,c ));
                                    
                                end
                            end
                            MR.Data = data;
                        end
                    end
                    
                end
                MR.Parameter.ReconFlags.isdepicorr = 1;
            end
        end
        function RingingFilter( MR )
            % RingingFilter: Removes ringing artifacts by applying a hamming filter on the k-space
            % data
            %
            % Syntax:     r.RingingFilter;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.RingingFilter')">Parameter.Recon.RingingFilter</a>: {'yes'} | 'no'
            %             Enables/disables the ringing filter.
            %           - <a href="matlab:helpwin('ReconParsDoc.RingingFilterStrength')">Parameter.Recon.RingingFilterStrength</a>:
            %             A value between 0 and 1 which specifies the strength of the ringing
            %             filter. 0 means no filter is applied and 1 is the maximum filter strength.
            %             The filter strength can be specified in every direction (Measurement,
            %             Phase Encoding, Slice Encoding).
            %
            % Location:   k-space
            %
            % Formats:    Raw | ExportedRaw
            %
            % Description/Algorithm: Sharp edges in k-space can result in ringing artifacts in
            %             image-space. Such edges are introduced by zero-filling or by sampling a
            %             low resultion image. The ringing filter tries the soften the edges by
            %             applying a hamming filter on the k-space data, lowering the values on the
            %             edge. The amount by which the edge-values are lowered is specified by the
            %             filter strength. A filter strength of 0.5 for example will bring down the
            %             edge value by 50%. A strength of brings the value down to 0 and therefore
            %             removing the edges completely. However a large filter strength might cause
            %             image blurring. Therefore it is important to find an optimal compromise
            %             between artifact removal and image blurring. The default filter strength
            %             is set to 0.25.
            %
            % Notes:      The ringing filter should be applied before k-space zero filling to work
            %             effectively.
            
            if strcmpi( MR.Parameter.Recon.RingingFilter, 'yes' )
                
                % Prerequisites: Check if the data has the right form
                if ~MR.Parameter.ReconFlags.isread
                    error( 'Please read the data first' );
                end
                if ~MR.Parameter.ReconFlags.issorted
                    error( 'Please sort the data first' );
                end
                if any( MR.Parameter.ReconFlags.isimspace )
                    error( 'The ringing filter can only be applied in k-space' );
                end
                
                if ~isempty(MR.Data)
                    % MR.Data can be a cell or an array depending on the dataset.
                    % To make our function work with all kinds of data, we convert
                    % MR.Data to a cell in any case and the let the hamming filter
                    % work on every cell element (using cellfun)
                    MR.Parameter.UpdateImageInfo = 0;
                    MR.DataClass.Convert2Cell;
                    
                    % For radial scans the zero filling is already done within
                    % the gridder. Therefore to apply the right filter size we
                    % have to know which part of the k-space has actually been
                    % sampled. This is stored in Kx, Ky, Kz Ranges.
                    if strcmpi( MR.Parameter.Gridder.Preset, 'radial' )
                        sampled_size = cell(size(MR.Data,1), size(MR.Data,2));
                        for ci = 1:size(MR.Data,1)
                            for cj = 1:size(MR.Data,2)
                                if isempty(MR.Data{ci,cj})
                                    sampled_size{ci, cj} = [];
                                    continue;
                                end
                                sampled_size{ci, cj} = [ min( [ size( MR.Data{ci,cj},1),  length( MR.Parameter.Encoding.WorkEncoding.KxRange{ci,cj}(1):MR.Parameter.Encoding.WorkEncoding.KxRange{ci,cj}(2) )]), ...
                                    min( [ size( MR.Data{ci,cj},2) length( MR.Parameter.Encoding.WorkEncoding.KxRange{ci,cj}(1):MR.Parameter.Encoding.WorkEncoding.KxRange{ci,cj}(2) )]) ];
                                if ~isempty( MR.Parameter.Encoding.WorkEncoding.KzRange{ci,cj} )
                                    sampled_size{ci, cj} = [sampled_size{ci,cj}, ...
                                        length( MR.Parameter.Encoding.WorkEncoding.KzRange{ci,cj}(1):MR.Parameter.Encoding.WorkEncoding.KzRange{ci,cj}(2) ) ];
                                else
                                    sampled_size{ci,cj} = [sampled_size{ci,cj}, 1];
                                end
                            end
                        end
                        
                    else
                        sampled_size = cell( size(MR.Data ) );
                    end
                    
                    MR.Data = cellfun( @(x, y)Filter.hamming_filter( x, MR.Parameter.Recon.RingingFilterStrength, y ), MR.Data, sampled_size, 'UniformOutput', 0 );
                    
                    % Convert MR.Data back to original state
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                end
                
            end
        end
        function ConcomitantFieldCorrection( MR )
            % ConcomitantFieldCorrection: Performes a concomitant field correction on the current data.
            %
            % Syntax:     r.ConcomitantFieldCorrection;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.ConcomitantFieldCorrection')">Parameter.Recon.ConcomitantFieldCorrection</a>: {'yes'} | 'no'
            %             Enables/disables the concomitant field correction.
            %           - <a href="matlab:helpwin('MRparameterDoc.Labels')">Parameter.Labels.ConcomFactors</a>: holds the correction factors.
            %           - Geometry parameters: The correction has to be performed in the scanner fixed
            %             xyz-coordinate-system. Therefore all necessary geometry parameters for the
            %             coordinate transformation (Transform function) have to be available.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.isconcomcorrected')">Parameter.ReconFlags.isconcomcorrected</a>
            %
            % Location:   Image-space.
            %
            % Formats:    Raw
            %
            % Description/Algorithm: Maxwell’s equations imply that imaging gradients are accompanied by
            %             higher order spatially varying fields (concomitant fields) that can cause
            %             artifacts in MR imaging. The lowest order concomitant fields depend
            %             quadratically on the imaging gradient amplitude and inversely on the
            %             static field strength.
            %             On the scanner these extra magnetic fields are expressed along the
            %             directions of the physical gradients (scanner fixed xyz-coordinate system).
            %             Therefore the images have to be transformed into that system before the
            %             correction. The correction itself adds a quadratic phase to the images
            %             defined by the correction parameters.
            %
            % Notes:    - The concomitant field correction will not work on data measured without the
            %             ReconFrame patch
            %           - The concomitant field correction is only executed for phase contrast flow scans.
            %             Only then correction parameters are provided by the scanner.
            
            if strcmpi( MR.Parameter.Recon.ConcomitantFieldCorrection, 'yes' )
                if isfield( MR.Parameter.Labels, 'ConcomFactors' ) && ...
                        any( MR.Parameter.Labels.ConcomFactors ~= 0 )
                    
                    if any( ~MR.Parameter.ReconFlags.isimspace )
                        error( 'Error in ConcomitantFieldCorrection: The ConcomitantFieldCorrection has to be performed in image space');
                    end
                    if MR.Parameter.ReconFlags.isconcomcorrected
                        error( 'Error in GeometryCorrection: The data is already concomitant field corrected' );
                    end
                    
                    MR.DataClass.Convert2Cell;
                    
                    A = MR.Transform( 'ijk', 'xyz' );
                    stacks_read = MR.Parameter.Parameter2Read.loca;
                    if isfield( MR.Parameter.Labels, 'StackIndex' )
                        stacks = MR.Parameter.Labels.StackIndex( stacks_read+1 );
                    end
                    
                    MR.Data(1,:) = cellfun( @(x)Flow.concom_corr( x, A, MR.Parameter.Labels.ConcomFactors, stacks, MR.Parameter.Parameter2Read.extr1, single(MR.Parameter.Labels.GeoCorrPars) ), ...
                        MR.Data( 1,:), 'UniformOutput', 0 );
                    
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    MR.Parameter.ReconFlags.isconcomcorrected = 1;
                end
            end
        end
        function GeometryCorrection( MR )
            % GeometryCorrection: Corrects for gradient non-linearities over large FOV's
            %
            % Syntax:     r.GeometryCorrection;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.GeometryCorrection')">Parameter.Recon.GeometryCorrection</a>: {'yes'} | 'no'
            %             Enables/disables the geometry correction.
            %           - <a href="matlab:helpwin('MRparameterDoc.Labels')">Parameter.Labels.GeoCorrPars</a>:
            %             Holds the correction factors.
            %           - Geometry parameters: The correction has to be performed in the scanner fixed
            %             xyz-coordinate-system. Therefore all necessary geometry parameters for the
            %             coordinate transformation (Transform function) have to be available.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.isgeocorrected')">Parameter.ReconFlags.isgeocorrected</a>
            %
            % Location:   Image-space.
            %
            % Formats:    Raw
            %
            % Description/Algorithm: Imaging gradients are not perfectly linear over the whole
            %             possible FOV. The further away we are from the isocentre the stronger the
            %             non-linear effects become, which will result in a slight distortion of the
            %             resulting image. This effect is corrected with the correction parameters
            %             provided by the scanner.
            %
            % Notes:    - The correction parameters are scanner specific and vary with the gradient
            %             system.
            
            if ~isempty(MR.Data) && strcmpi( MR.Parameter.Recon.GeometryCorrection, 'yes' )
                if ~isfield( MR.Parameter.Labels, 'GeoCorrPars' )
                    warning( 'MATLAB:MRecon', 'Cannot execute the geometry correction. Not enough parameter. Please update the scanner patch' );
                    return;
                end
                if any( ~MR.Parameter.ReconFlags.isimspace )
                    error( 'Error in GeometryCorrection: The data must be in image-space' );
                end
                if MR.Parameter.ReconFlags.isgeocorrected
                    error( 'Error in GeometryCorrection: The data is already geometry corrected' );
                end
                
                if ~isempty( MR.Parameter.Scan.ijk ) && ~isempty( MR.Parameter.Scan.xyz )
                    if isfield( MR.Parameter.Labels, 'StackIndex' )
                        try
                            stack_nr = MR.Parameter.Labels.StackIndex( MR.Parameter.Parameter2Read.loca + 1 );
                        catch
                            stack_nr = [];
                        end
                    else
                        stack_nr = [];
                    end
                    
                    MR.Parameter.UpdateImageInfo = 0;
                    MR.DataClass.Convert2Cell;
                    
                    MR.Data = cellfun( @(x)GeoCorr.geometry_correction( x, MR.Parameter, stack_nr),  ...
                        MR.Data, 'UniformOutput', 0 );
                    
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    MR.Parameter.ReconFlags.isgeocorrected = 1;
                end
            end
        end
        function FlowPhaseCorrection( MR )
            if ~isempty(MR.Data) && strcmpi( MR.Parameter.Recon.FlowPhaseCorrection, 'yes' )
                if any( ~MR.Parameter.ReconFlags.isimspace )
                    error( 'Error in FlowPhaseCorrection: The FlowPhaseCorrection has to be performed in image space');
                end
                
                if ~isempty( MR.Parameter.Scan.Venc ) &&  ...
                        any( any( MR.Parameter.Scan.Venc ~= 0 )) && ...
                        strcmpi( MR.Parameter.Recon.CoilCombination, 'pc' )
                    
                    MR.DataClass.Convert2Cell;
                    MR.Data(1,:) = cellfun( @(x)Flow.fit_flow_phase( x ), MR.Data(1,:) , 'UniformOutput', 0);
                    MR.Data = Helper.UnconvertCell( MR.Data );
                end
            end
        end
        
        % ---------------------------------------------------------------%
        % Data Gridding
        % ---------------------------------------------------------------%
        function GridData( MR )
            % GridData: Grids k-space data on a predefined trajectory.
            %
            % Syntax:     r.GridData;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.Gridding')">Parameter.Recon.Gridding</a>: {'yes'} | 'no'
            %             Enables/disables the gridding.
            %           - <a href="matlab:helpwin('GridderParsDoc')">Parameter.Gridder</a>:
            %             The gridder parameter.
            %           - <a href="matlab:helpwin('MRparameterDoc.Labels')">Parameter.Labels.NusEncNrs</a>:
            %             The non-uniform sample position for EPI scans. In EPI acquisitions the
            %             samples in readout direction are not uniformly spread due to the sampling
            %             during gradient ramp up. Therefore EPI scans are usually gridded in
            %             readout direction using the k-space positions specified here.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.isgridded')">Parameter.ReconFlags.isgridded</a>
            %
            % Location:   k-space after SortData.
            %
            % Formats:    Raw | ExportedRaw
            %
            % Description/Algorithm: The fast fourier transform (FFT) requires the data to be on a
            %             regular cartesian grid. Therefore in order to use the FFT on non-cartesian
            %             data such as radial or spiral acquistions, we have to grid the k-space
            %             data from the nominal trajectory to a cartesian grid. There are several
            %             predefined trajectories already built-in into MRecon, including:
            %             2D-radial, 3D-stack of radials, 3D-kooshball, 2D/3D spiral, EPI ramp
            %             sampling. On the other side, user defined trajectories can easily be
            %             defined by setting the <a href="matlab:helpwin('GridderParsDoc')">Gridder Parameters</a>. Please see the example below
            %             on how to to that.
            %             If the <a href="matlab:helpwin('GridderParsDoc')">Gridder Parameters</a> are not set manually MRecon automatically detects
            %             the type of the scan and decides if it should be gridded or not.
            %
            % Notes:    - For non-cartesian scans the zero-filling process in k-space and
            %             oversampling removal is automatically performed during gridding, by
            %             chosing the grid appropriately. See example 3)
            %           - The Gridder parameter can also be filled using the
            %             <a href="matlab:helpwin('MRecon.GridderCalculateTrajectory')">GridderCalculateTrajectory</a> function.
            %
            % Examples:   1) Reconstruct a golden angle radial acquisition by changing the radial
            %                angles in the gridder parameters:
            %
            %                Golden angle scans exhibit a radial angle increment of 111.246°, while
            %                the rest is exactly the same as a "normal" radial trajectory. Therefore
            %                all we have to do to reconstruct a golden angle scan is to change the
            %                radial angles in the gridder parameter:
            %
            %                Create a new reconstruction object with a golden angle scan:
            %
            %                   r = MRecon('golden_angle_scan.raw');
            %
            %                Calculate the total number of phase encoding profiles:
            %
            %                   nr_profs = r.Parameter.Scan.Samples(2)*r.Parameter.Scan.Samples(3)
            %
            %                Set the radial angles to an increment of 111.246°
            %
            %                   r.Parameter.Gridder.RadialAngles = 0:pi/180*111.246:(nr_profs-1)*pi/180*111.246;
            %
            %                Perform the recon
            %
            %                   r.Perform;
            %
            %             2) Define a completely new trajectory and weights.
            %
            %                In this example we mirror the final image in readout direction by
            %                gridding a cartesian scan to a inverted kx trajectory. To do that we
            %                define a completely new trajectory and pass it to the gridder
            %                parameters:
            %
            %                Create a new reconstruction object using a cartesian scan:
            %
            %                   r = MRecon('cartesian.raw');
            %
            %                Obtain the number of samples in every direction:
            %
            %                   no_kx_samples = r.Parameter.Scan.Samples(1);
            %                   no_ky_samples = r.Parameter.Scan.Samples(2);
            %                   no_kz_samples = r.Parameter.Scan.Samples(3);
            %
            %                Create a regular cartesian grid. These are the sampled k-space
            %                coordinates of a cartesian scan. By definition the k-space coordinates
            %                of a sampled trajectory are in the range of -floor(n/2):ceil(n/2)-1,
            %                where n is the number of samples:
            %
            %                   [kx,ky,kz] = ndgrid( -floor(no_kx_samples / 2):ceil( no_kx_samples /2)-1, ...
            %                       -floor(no_ky_samples / 2):ceil( no_ky_samples /2)-1, ...
            %                       -floor(no_kz_samples / 2):ceil( no_kz_samples /2)-1 );
            %
            %                Modify the k-space grid by inverting the readout direction:
            %
            %                   kx_inverted = -kx;
            %
            %                Create the modified k-space. These are the coordinates we grid on. By
            %                definition the coordinates must be of dimension: no_samples x
            %                no_y_profiles x no_z_profiles x 3 (See <a href="matlab:helpwin('GridderParsDoc')">Gridder Parameters</a> for more
            %                information):
            %
            %                   k = cat( 4, kx_inverted, ky, kz );
            %
            %                Assign the new trajectory to the Gridder:
            %
            %                   r.Parameter.Gridder.Kpos = k;
            %
            %                Assign gridder weights (sampling density). In this example the samples
            %                are uniformly distributed on the whoe k-space (all weights are set to 1):
            %
            %                   r.Parameter.Gridder.Weights = ones(no_kx_samples, no_ky_samples, no_kz_samples);
            %
            %                Perform the reconstruction:
            %
            %                   r.Perform;
            %
            %                Show the data. Observe that the image is flipped in readout direction
            %                compared to the ungridded reconstruction (a flip in k-space corresponds
            %                to a flip in image-space):
            %
            %                   r.ShowData;
            %
            %             3) Gridder and matrix sizes:
            %
            %                Oversampling:
            %
            %                During the gridding process we have the possibility to sample on an
            %                arbitrary grid. This means that we can for example sample on a finer or
            %                coarser grid than the acquired one, which corresponds to add or remove
            %                oversampling and therefore changing the FOV in the reconstructed image. The
            %                parameter which specifies the amount of oversampling introduced by the
            %                gridder is <a href="matlab:helpwin('GridderParsDoc.GridOvsFactor')">Parameter.Gridder.GridOvsFactor</a>.
            %                For example if we set:
            %
            %                   r.Parameter.Gridder.GridOvsFactor = 2;
            %
            %                then we grid to half the sampled k-space distance and therefore to
            %                twice the FOV. A value smaller than 1 corresponds to a coarser grid
            %                than acquired and therefore a smaller FOV.
            %
            %                In radial scans the acquired oversampling factor in readout direction
            %                is always 2 (hardcoded in the pulse programming enviroment). However If
            %                we check the oversampling factor in r.Parameter.Encoding.KxOversampling
            %                we ca usually find a value like 1.25. This is the oversampling factor
            %                AFTER gridder, which means that we should set:
            %
            %                   r.Parameter.Gridder.GridOvsFactor = r.Parameter.Encoding.KxOversampling / 2;
            %
            %                This is done automatically done for radial scans, and means that we
            %                grid to a coarser grid leaving an oversampling factor of 1.25.
            %
            %                Zero Filling:
            %                Apart from gridding on a finer grid we also have the possibility to
            %                grid on a larger/smaller k-space and therefore changing the resolution
            %                of the scan. Gridding to a larger k-space means that we perform zero
            %                filling whlie gridding to a smaller one is reducing the resolution of
            %                the scan. The parameter which specifies the size of the k-space to be
            %                gridded on is: <a href="matlab:helpwin('GridderParsDoc.OutputMatrixSize')">Parameter.Gridder.OutputMatrixSize</a>.
            %
            %                The output size is the matrix size after gridding. Lets assume we have
            %                a radial scan with 256 acquired samples. If we set:
            %
            %                   r.Parameter.Gridder.GridOvsFactor = 1;
            %                   r.Parameter.Gridder.OutputMatrixSize = [512,512,1];
            %
            %                The the output matrix has the size 512x512 which means that we have
            %                zero filled the 256x256 sampled k-space to that value. However if we
            %                set:
            %
            %                   r.Parameter.Gridder.GridOvsFactor = 2;
            %                   r.Parameter.Gridder.OutputMatrixSize = [512,512,1];
            %
            %                Then we have gridded it to a grid twice as fine but NOT zero-filled. To
            %                zero-fill and oversample we have to set:
            %
            %                   r.Parameter.Gridder.GridOvsFactor = 2;
            %                   r.Parameter.Gridder.OutputMatrixSize = [1024,1024,1];
            %
            %                Summary:
            %                For radial scans the zero filling and oversampling removal usually
            %                takes place directly in the gridder. Below is an example from a real
            %                radial experiment:
            %
            %                   Number of acquired samples: 464 (oversampling factor = 2)
            %                   r.Parameter.Encoding.KxOversampling = 1.25
            %                   r.Parameter.Encoding.XRes = 240
            %
            %                   --> Matrix size after gridding = 240 * 1.25 = 300
            %                   --> Gridding oversampling = 1.25 / 2 = 0.625
            %
            %               --> Set gridder parameters to:
            %
            %                   r.Parameter.Gridder.GridOvsFactor = 0.625;
            %                   r.Parameter.Gridder.OutputMatrixSize = [300,300,1];
            
            if strcmpi( MR.Parameter.Recon.Gridding, 'yes' )
                if ~MR.Parameter.ReconFlags.issorted && ~strcmpi( MR.Parameter.Gridder.Preset, 'Epi' )
                    error( 'Error in Gridder: Please sort the data first' );
                end
                if MR.Parameter.ReconFlags.isgridded
                    error( 'Error in Gridder: The data is already gridded' );
                end
                if any( MR.Parameter.ReconFlags.isimspace )
                    error( 'Error in Gridder: Please grid the data in k-space' );
                end
                
                if any( strcmpi( MR.Parameter.Gridder.Preset, {'Radial', 'Spiral', 'Epi', 'Cartesian'} )) || ~isempty( MR.Parameter.Gridder.Kpos )
                    if ~isempty(MR.Data)
                        MR.DataClass.Convert2Cell;
                        MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                        [MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable] = Helper.check_cell_sizes( MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable, MR.Data );
                        
                        MR.Parameter.UpdateImageInfo = 0;
                        
                        MR.Parameter.Gridder.InitWorkingPars(size(MR.Data, 2));
                        
                        for cell_index = 1:size(MR.Data, 2)
                            if strcmpi( MR.Parameter.Gridder.Preset, 'Spiral') && ~isempty(find(MR.Parameter.Parameter2Read_Original.typ == Label.MRECON_TYPE_TRAJ_DATA, 1))
                                MR.Parameter.Gridder.Trajectory = Gridder.read_trajectory_from_file(MR.Parameter.Filename, MR.Parameter.DataType, MR.Parameter.Labels);
                            end
                            
                            MR.Parameter.Gridder.AssignWorkingPars(cell_index);
                            
                            if ~strcmpi( MR.Parameter.Gridder.WorkingPars{cell_index}.Preset, 'none' )
                                grid_ovs = [];
                                
                                if isempty( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos )
                                    if isfield( MR.Parameter.Labels, 'NusEncNrs' )
                                        nus_enc_nrs = MR.Parameter.Labels.NusEncNrs;
                                        if size(nus_enc_nrs, 1) < size(nus_enc_nrs, 2)
                                            nus_enc_nrs = nus_enc_nrs';
                                        end
                                    else
                                        nus_enc_nrs = [];
                                    end
                                    [MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, grid_ovs, MR.Parameter.Gridder.WorkingPars{cell_index}.RadialAngles] = ...
                                        GridderTrajectory.calculate_trajectory( MR.Data{1,cell_index}, MR.Parameter.Encoding.WorkEncoding.KxOversampling{cell_index}, ...
                                        MR.Parameter.Encoding.WorkEncoding.KxRange{1,cell_index}, MR.Parameter.Gridder.WorkingPars{cell_index}, ...
                                        nus_enc_nrs, MR.Parameter.Scan.SENSEFactor, MR.Parameter.Encoding.WorkEncoding.Echo{cell_index}, ...
                                        MR.Parameter.Labels.Samples);                                       %
                                else
                                    if ndims( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos ) < 4 || ( size(MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, 4) ~= 3)
                                        error( 'Error in Gridder: Parameter.Gridder.Kpos has to be a 4-dimensional matrix of size no_samples x no_profiles x no_slices x 3' );
                                    end
                                end
                                
                                if isempty( MR.Parameter.Gridder.WorkingPars{cell_index}.GridOvsFactor )
                                    MR.Parameter.Gridder.WorkingPars{cell_index}.GridOvsFactor = grid_ovs;
                                end
                                
                                if isempty( MR.Parameter.Gridder.WorkingPars{cell_index}.Weights )
                                    MR.Parameter.Gridder.WorkingPars{cell_index}.Weights = GridderTrajectory.calculate_weights( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, MR.Parameter.Gridder.WorkingPars{cell_index});
                                end
                                
                            end
                            
                            if ~isempty( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos )
                                
                                if isempty( MR.Parameter.Gridder.WorkingPars{cell_index}.OutputMatrixSize )
                                    MR.Parameter.Gridder.WorkingPars{cell_index}.OutputMatrixSize = [MR.Parameter.Encoding.WorkEncoding.XRes{1,cell_index}*MR.Parameter.Encoding.WorkEncoding.KxOversampling{1,cell_index}, ...
                                        MR.Parameter.Encoding.WorkEncoding.YRes{1,cell_index}*MR.Parameter.Encoding.WorkEncoding.KyOversampling{1,cell_index}, ...
                                        MR.Parameter.Encoding.WorkEncoding.ZRes{1,cell_index}*MR.Parameter.Encoding.WorkEncoding.KzOversampling{1,cell_index}];
                                    
                                    
                                    if strcmpi( MR.Parameter.Gridder.WorkingPars{cell_index}.Preset, 'Epi' ) && isfield( MR.Parameter.Labels, 'Samples' )
                                        MR.Parameter.Gridder.WorkingPars{cell_index}.OutputMatrixSize = [round(MR.Parameter.Labels.Samples(1)*MR.Parameter.Encoding.WorkEncoding.KxOversampling{1,cell_index}),...
                                            round(MR.Parameter.Labels.Samples(2)*MR.Parameter.Encoding.WorkEncoding.KyOversampling{1,cell_index}), ...
                                            round(MR.Parameter.Labels.Samples(3)*MR.Parameter.Encoding.WorkEncoding.KzOversampling{1,cell_index})];
                                    end
                                end
                                
                                if isempty( MR.Parameter.Gridder.WorkingPars{cell_index}.GridOvsFactor )
                                    MR.Parameter.Gridder.WorkingPars{cell_index}.GridOvsFactor = 1;
                                end
                                
                                if isempty(MR.Parameter.Gridder.WorkingPars{cell_index}.Weights)
                                    MR.Parameter.Gridder.WorkingPars{cell_index}.Weights = ones( size(MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos,1),size(MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos,2),size(MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos,3));
                                end
                                
                                if (size( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, 1 ) ~= size( MR.Parameter.Gridder.WorkingPars{cell_index}.Weights, 1 )) || ...
                                        (size( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, 2 ) ~= size( MR.Parameter.Gridder.WorkingPars{cell_index}.Weights, 2 )) || ...
                                        (size( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, 3 ) ~= size( MR.Parameter.Gridder.WorkingPars{cell_index}.Weights, 3 ))
                                    error( 'Error in Gridder: Parameter.Gridder.Weights has to be a matrix of size: %d x %d x %d', size( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, 1 ), size( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, 2 ), size( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, 3) );
                                end
                                if strcmpi( MR.Parameter.Gridder.WorkingPars{cell_index}.Preset, 'radial' ) && ~isempty( MR.Parameter.Gridder.WorkingPars{cell_index}.RadialAngles ) && ( size( MR.Parameter.Gridder.WorkingPars{cell_index}.RadialAngles, 1 ) ~= size( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, 2 ) )
                                    error( 'Error in Gridder: The number of rows in Parameter.Gridder.RadialAngles has to be %d', size( MR.Parameter.Gridder.WorkingPars{cell_index}.Kpos, 2 ) );
                                end
                                
                            end
                            
                            % Check if Gridder parameter have the right sizes
                            if size(MR.Data(1,:), 2) ~= length( MR.Parameter.Gridder.WorkingPars )
                                error( 'Error in Gridder: Parameter.Gridder.WorkingPars has to be a cell array of size %d x %d', size(MR.Data(1,:),1), size(MR.Data(1,:),2) );
                            end
                        end
                        
                        
                        MR.Data(1,:) = cellfun( @(x, y)Gridder.grid_data( x, y ), ...
                            MR.Data(1,:), MR.Parameter.Gridder.WorkingPars, 'UniformOutput', false );
                        
                        % Delete label lookup table
                        if ~strcmpi( MR.Parameter.Gridder.Preset, 'Epi' )
                            MR.Parameter.LabelLookupTable = cellfun( @(x)[], MR.Parameter.LabelLookupTable , 'UniformOutput',0);
                        end
                        
                        MR.Parameter.ReconFlags.isgridded = 1;
                        
                        
                        MR.Parameter.UpdateImageInfo = 1;
                        MR.Data = Helper.UnconvertCell( MR.Data );
                        MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                    end
                end
            end
        end
        function GridderCalculateTrajectory( MR )
            % GridderCalculateTrajectory: Calculated the gridder trajectory for the current scan and
            % fills the <a href="matlab:helpwin('GridderParsDoc')">gridder parameter</a>
            %
            % Syntax:     r.GridderCalculateTrajectory;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('GridderParsDoc')">Parameter.Gridder</a>:
            %             The gridder parameter.
            %           - <a href="matlab:helpwin('MRparameterDoc.Labels')">Parameter.Labels.NusEncNrs</a>:
            %             The non-uniform sample position for EPI scans. In EPI acquisitions the
            %             samples in readout direction are not uniformly spread due to the sampling
            %             during gradient ramp up. Therefore EPI scans are usually gridded in
            %             readout direction using the k-space positions specified here.
            %
            % Location:   k-space after SortData.
            %
            % Formats:    Raw | ExportedRaw
            %
            % Description/Algorithm: When the gridder parameter are not defined by the user, they
            %             are calculated internally and will not be visible or changable. If one
            %             wants to observe or change the gridding parameter, this function can be
            %             called which calculates the trajectory for the gridding process and fills
            %             the gridder parameter accordingly. The function can be called for any
            %             preset scan (radial, spiral, EPI), before calling the GridData function.
            %             This is useful if the user wants to modify a preset trajectory without
            %             calculating it by himself. To do so one would call this function,
            %             modify the gridder parameter followed by the call of <a href="matlab:helpwin('MRecon.GridData')">GridData</a>.
            
            if ~MR.Parameter.ReconFlags.isread
                error( 'Error in CalculateTrajectory: Please read the data first' );
            end
            if ~MR.Parameter.ReconFlags.issorted
                error( 'Error in CalculateTrajectory: Please sort the data first' );
            end
            
            if ~strcmpi( MR.Parameter.Gridder.Preset, 'none' )
                MR.DataClass.Convert2Cell;
                grid_ovs = [];
                
                if strcmpi( MR.Parameter.Gridder.Preset, 'Spiral') && ~isempty(find(MR.Parameter.Parameter2Read_Original.typ == 4, 1))
                    MR.Parameter.Gridder.Trajectory = Gridder.read_trajectory_from_file(MR.Parameter.Filename, MR.Parameter.DataType, MR.Parameter.Labels);
                end
                
                if isfield( MR.Parameter.Labels, 'NusEncNrs' )
                    nus_enc_nrs = MR.Parameter.Labels.NusEncNrs;
                    if size(nus_enc_nrs, 1) < size(nus_enc_nrs, 2)
                        nus_enc_nrs = nus_enc_nrs';
                    end
                else
                    nus_enc_nrs = [];
                end
                
                if isempty( MR.Parameter.Gridder.Kpos )
                    [MR.Parameter.Gridder.Kpos, grid_ovs, MR.Parameter.Gridder.RadialAngles] = cellfun( @(x,y,z,v)GridderTrajectory.calculate_trajectory( x, y, z, MR.Parameter.Gridder, nus_enc_nrs, MR.Parameter.Scan.SENSEFactor, v, MR.Parameter.Labels.Samples ), ...
                        MR.Data(1,:), MR.Parameter.Encoding.WorkEncoding.KxOversampling(1,:), ...
                        MR.Parameter.Encoding.WorkEncoding.KxRange(1,:) , MR.Parameter.Encoding.WorkEncoding.Echo(1), 'UniformOutput', false );
                else
                    MR.Parameter.Gridder.Kpos = Helper.Convert2Cell( MR.Parameter.Gridder.Kpos );
                end
                
                if isempty( MR.Parameter.Gridder.GridOvsFactor )
                    MR.Parameter.Gridder.GridOvsFactor = grid_ovs;
                else
                    MR.Parameter.Gridder.GridOvsFactor = Helper.Convert2Cell( MR.Parameter.Gridder.GridOvsFactor );
                end
                
                if isempty( MR.Parameter.Gridder.Weights )
                    MR.Parameter.Gridder.Weights = cellfun( @(x,y,z,v)GridderTrajectory.calculate_weights( x, MR.Parameter.Gridder), ...
                        MR.Parameter.Gridder.Kpos, 'UniformOutput', false );
                else
                    MR.Parameter.Gridder.Weights = Helper.Convert2Cell( MR.Parameter.Gridder.Weights );
                end
                
                if isempty( MR.Parameter.Gridder.OutputMatrixSize )
                    MR.Parameter.Gridder.OutputMatrixSize = cellfun( @(x,ox, y, oy, z, oz)[x*ox, y*oy, z*oz], ...
                        MR.Parameter.Encoding.WorkEncoding.XRes(1,:), MR.Parameter.Encoding.WorkEncoding.KxOversampling(1,:), ...
                        MR.Parameter.Encoding.WorkEncoding.YRes(1,:), MR.Parameter.Encoding.WorkEncoding.KyOversampling(1,:), ...
                        MR.Parameter.Encoding.WorkEncoding.ZRes(1,:), MR.Parameter.Encoding.WorkEncoding.KzOversampling(1,:), ...
                        'UniformOutput', false);
                    
                    if strcmpi( MR.Parameter.Gridder.Preset, 'Epi' ) && isfield( MR.Parameter.Labels, 'Samples' )
                        MR.Parameter.Gridder.OutputMatrixSize = cellfun( @(ox,oy,oz)[round(MR.Parameter.Labels.Samples(1)*ox), round(MR.Parameter.Labels.Samples(2)*oy), round(MR.Parameter.Labels.Samples(3)*oz)], ...
                            MR.Parameter.Encoding.WorkEncoding.KxOversampling(1,:), ...
                            MR.Parameter.Encoding.WorkEncoding.KyOversampling(1,:), ...
                            MR.Parameter.Encoding.WorkEncoding.KzOversampling(1,:), ...
                            'UniformOutput', false);
                    end
                else
                    MR.Parameter.Gridder.OutputMatrixSize = Helper.Convert2Cell( MR.Parameter.Gridder.OutputMatrixSize );
                end
                
                for cell_index = 1:size(MR.Data, 2)
                    MR.Parameter.Gridder.AssignWorkingPars(cell_index);
                end
                
                MR.Parameter.Gridder.Kpos = Helper.UnconvertCell(MR.Parameter.Gridder.Kpos);
                MR.Parameter.Gridder.RadialAngles = Helper.UnconvertCell(MR.Parameter.Gridder.RadialAngles);
                MR.Parameter.Gridder.Weights = Helper.UnconvertCell(MR.Parameter.Gridder.Weights);
                MR.Parameter.Gridder.GridOvsFactor = Helper.UnconvertCell(MR.Parameter.Gridder.GridOvsFactor);
                MR.Parameter.Gridder.OutputMatrixSize = Helper.UnconvertCell(MR.Parameter.Gridder.OutputMatrixSize);
                
                MR.Data = Helper.UnconvertCell( MR.Data );
                
            end
        end
        function GridderNormalization( MR )
            % GridderNormalization: Corrects for the signal weighting in image-space originating
            % from the convolution with a bessel kernel in the gridding process.
            %
            % Syntax:     r.GridderNormalization;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('GridderParsDoc')">Parameter.Gridder</a>:
            %             The gridder parameter.
            %           - <a href="matlab:helpwin('ReconFlags.isgridded')">Parameter.ReconFlags.isgridded</a>:
            %             Used to check if the data was gridded.
            %           - <a href="matlab:helpwin('EncodingPars.WorkEncoding.KxOversampling')">Parameter.Encoding.WorkEncoding.KxOversampling</a>:
            %             The resulting oversampling after gridding, needed in the normalization.
            %
            % Location:   Image-space.
            %
            % Formats:    Raw | ExportedRaw
            %
            % Description/Algorithm: Gridding is essentially a convolution with a bessel kernel in
            %             k-space. In image-space this corresponds to a multiplication with the
            %             fourier transform of the bessel-kernel leaving a signal weighing on the
            %             reconstructed image. This function corrects for that by dividing the image
            %             by the fourier-transform of the bessel kernel.
            
            if strcmpi( MR.Parameter.Recon.Gridding, 'yes' ) && MR.Parameter.ReconFlags.isgridded && ~strcmpi( MR.Parameter.Gridder.Preset, 'epi' )
                if ~strcmpi( MR.Parameter.Gridder.Preset, 'none' )
                    if ~all( MR.Parameter.ReconFlags.isimspace )
                        error( 'The data can only be normalized in image-space' );
                    end
                    if ~MR.Parameter.ReconFlags.isgridded
                        error( 'Please grid the data first' );
                    end
                    
                    if ~isempty(MR.Data)
                        if strcmpi( MR.Parameter.Gridder.Normalize, 'yes' )
                            MR.Parameter.UpdateImageInfo = 0;
                            MR.DataClass.Convert2Cell;
                            
                            % Check if the Gridder gridded a 3d or 2d trajectory (for a 2d
                            % trajectory we dont have to do the normalizazion in slice encoding
                            % direction)
                            is_3d = cell(1, length(MR.Parameter.Gridder.WorkingPars));
                            for i = 1:length(MR.Parameter.Gridder.WorkingPars)
                                is_3d{i} = size(MR.Parameter.Gridder.WorkingPars{i}.Kpos, 3) > 1;
                            end
                            
                            MR.Data(1,:) = cellfun( @(x, y, ox, oy, oz, xr, yr, zr)Gridder.gridder_normalization( x, MR.Parameter.Gridder.KernelWidth, y, ox, oy, oz, xr, yr, zr), ...
                                MR.Data(1,:), is_3d, MR.Parameter.Encoding.WorkEncoding.KxOversampling(1,:), ...
                                MR.Parameter.Encoding.WorkEncoding.KyOversampling(1,:), MR.Parameter.Encoding.WorkEncoding.KzOversampling(1,:), ...
                                MR.Parameter.Encoding.WorkEncoding.XRange(1,:), MR.Parameter.Encoding.WorkEncoding.YRange(1,:), MR.Parameter.Encoding.WorkEncoding.ZRange(1,:), ...
                                'UniformOutput', false );
                            
                            MR.Parameter.UpdateImageInfo = 1;
                            MR.Data = Helper.UnconvertCell( MR.Data );
                            
                        end
                    end
                end
            end
        end
        
        % ---------------------------------------------------------------%
        % k-space / image-space Transformation
        % ---------------------------------------------------------------%
        function K2I( MR )
            % K2I: Performs a fourier transformation from k-space to image-space
            %
            % Syntax:     r.K2I;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('EncodingPars.FFTDims')">Parameter.Encoding.FFTDims</a>:
            %             Specifies in which dimension the fourier transformation should be
            %             executed. This is a 3 elements vector specifying the dimension in matrix
            %             notation.
            %           - <a href="matlab:helpwin('ReconFlags.isimspace')">Parameter.ReconFlags.isimspace</a>:
            %             Used to check which dimensions are still in k-space
            %           - <a href="matlab:helpwin('EncodingPars.FFTShift')">Parameter.Encoding.FFTShift</a>:
            %             Specifies in which dimension the image should be shifted after the fourier
            %             transform. This is a 3 elements vector specifying the shift in matrix
            %             notation.
            %           - <a href="matlab:helpwin('EncodingPars.XRange')">Parameter.Encoding.XRange</a>, <a href="matlab:helpwin('EncodingPars.YRange')">Parameter.Encoding.YRange</a>, <a href="matlab:helpwin('EncodingPars.ZRange')">Parameter.Encoding.ZRange</a>:
            %             Define the shifts of the image after the fourier transform, for every
            %             dimension.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlags.isimspace')">Parameter.ReconFlags.isimspace</a>
            %
            % Location:   k-space
            %
            % Formats:    Raw | ExportedRaw | Cpx | Rec | Bruker
            %
            % Description/Algorithm: K2I first checks which dimensions should be transformed,
            %             defined in <a href="matlab:helpwin('EncodingPars.FFTDims')">Parameter.Encoding.FFTDims</a>. Afterwards it checks which of
            %             these dimensions are in k-space and then only transforms these ones. The
            %             fourier transformation along one dimension is given by:
            %
            %                   img = sqrt(n).* fftshift( ifft(ifftshift(kspace ,dim),[],dim), dim);
            %
            %             where n are the number of samples along the specific dimension.
            %
            %             After the fourier transform from k- to image-space, the images are shifted
            %             in image-space as defined in Parameter.Encoding.X/Y/ZRange. Thereby the
            %             images are shifted by the amount which leads to a symmetric range around 0.
            %             E.g. if the range is given as [-42, 21] then the image is shifted by 10
            %             pixel to the right. ( [-42, 21] + 10 = [-32, 31] )
            
            if all( MR.Parameter.ReconFlags.isimspace )
                error( 'The data is already in image space');
            end
            
            if ~isempty(MR.Data)
                
                fftdims = find( MR.Parameter.Encoding.WorkEncoding.FFTDims & ...
                    MR.Parameter.ReconFlags.isimspace == 0);
                
                MR.Parameter.UpdateImageInfo = 0;
                MR.DataClass.Convert2Cell;
                
                % spectro begin ----------------------------
                if MR.isSpectro
                    MR.Data = cellfun( @(x)FFT.k2i(x, fftdims, 1, true), MR.Data, 'UniformOutput', false );
                else
                    % spectro end ------------------------------
                    MR.Data = cellfun( @(x)FFT.k2i(x, fftdims, 1), MR.Data, 'UniformOutput', false );
                end
                
                if any( MR.Parameter.Encoding.WorkEncoding.FFTShift ) && any( strcmpi( MR.Parameter.DataFormat,{'ExportedRaw', 'Raw' }))
                    
                    if ~MR.Parameter.Encoding.WorkEncoding.FFTShift(1) || ~MR.Parameter.Encoding.WorkEncoding.FFTDims(1)
                        xrange = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.XRange, 'UniformOutput',0);
                    else
                        xrange = MR.Parameter.Encoding.WorkEncoding.XRange;
                        MR.Parameter.Encoding.WorkEncoding.FFTShift(1) = 0;
                    end
                    if ~MR.Parameter.Encoding.WorkEncoding.FFTShift(2) || ~MR.Parameter.Encoding.WorkEncoding.FFTDims(2)
                        yrange = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.YRange, 'UniformOutput',0);
                    else
                        yrange = MR.Parameter.Encoding.WorkEncoding.YRange;
                        MR.Parameter.Encoding.WorkEncoding.FFTShift(2) = 0;
                    end
                    
                    % If the image-space shift has been disabled then pass an empty range to the
                    % shift function (no shift)
                    if ~MR.Parameter.Encoding.WorkEncoding.FFTShift(3) || ~MR.Parameter.Encoding.WorkEncoding.FFTDims(3)
                        zrange = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.ZRange, 'UniformOutput',0);
                        zres = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.ZRes, 'UniformOutput',0);
                    else
                        zrange = MR.Parameter.Encoding.WorkEncoding.ZRange;
                        zres = MR.Parameter.Encoding.WorkEncoding.ZRes;
                        MR.Parameter.Encoding.WorkEncoding.FFTShift(3) = 0;
                    end
                    
                    MR.Data = cellfun( @(x, xr, yr, zr, zre)FFT.shift_image( x, xr, yr, zr, zre) , ...
                        MR.Data, ...
                        xrange, ...
                        yrange, ...
                        zrange, ...
                        zres, 'UniformOutput', 0 );
                end
                
                MR.Parameter.UpdateImageInfo = 1;
                MR.Data = Helper.UnconvertCell( MR.Data );
                
                MR.Parameter.ReconFlags.isimspace(fftdims) = 1;
            end
            
        end
        function K2IM( MR )
            % K2IM: Performs a fourier transformation from k-space to image-space in measurement
            % direction
            %
            % Syntax:     r.K2IM;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('EncodingPars.FFTDims')">Parameter.Encoding.FFTDims</a>:
            %             Specifies if the FFT in measurement direction should be executed. This is
            %             a 3 elements vector specifying the dimension in matrix notation.
            %           - <a href="matlab:helpwin('ReconFlags.isimspace')">Parameter.ReconFlags.isimspace</a>:
            %             Used to check if the measurement direction is still in k-space
            %           - <a href="matlab:helpwin('EncodingPars.FFTShift')">Parameter.Encoding.FFTShift</a>:
            %             Specifies if the image should be shifted in measurement direction after
            %             the fourier transform. This is a 3 elements vector specifying the
            %             shift in matrix notation.
            %           - <a href="matlab:helpwin('EncodingPars.XRange')">Parameter.Encoding.XRange</a>:
            %             Define the shifts of the image after the fourier transform in measurement
            %             direction.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlags.isimspace')">Parameter.ReconFlags.isimspace</a>
            %
            % Location:   k-space
            %
            % Formats:    Raw | ExportedRaw | Cpx | Rec | Bruker
            %
            % Description/Algorithm: K2IM first checks if the measurement direction should be
            %             transformed, defined in <a href="matlab:helpwin('EncodingPars.FFTDims')">Parameter.Encoding.FFTDims</a>. Afterwards it
            %             checks if the measurement direction is still in k-space and if yes
            %             transforms it. The fourier transformation along one dimension is given by:
            %
            %                   img = sqrt(n).* fftshift( ifft(ifftshift(kspace ,dim),[],dim), dim);
            %
            %             where n are the number of samples along the specific dimension.
            %
            %             After the fourier transform from k- to image-space, the image is shifted
            %             in measurement direction as defined in Parameter.Encoding.XRange.
            %             Thereby the images are shifted by the amount which makes the range symmetric
            %             around 0. E.g. if the range is given as [-42, 21] then the image is
            %             shifted by 10 pixel to the right. ( [-42, 21] + 10 = [-32, 31] )
            if MR.Parameter.ReconFlags.isimspace(1)
                error( 'The data has already been fourier transformed in measurement direction');
            end
            
            if ~isempty(MR.Data)
                MR.Parameter.UpdateImageInfo = 0;
                MR.DataClass.Convert2Cell;
                
                % spectro begin ----------------------------
                if MR.isSpectro
                    MR.Data = cellfun( @(x)FFT.k2i(x, 1, 1, true), MR.Data, 'UniformOutput', false );
                else
                    % spectro end ------------------------------
                    MR.Data = cellfun( @(x)FFT.k2i(x, 1, 1), MR.Data, 'UniformOutput', false );
                end
                
                if MR.Parameter.Encoding.WorkEncoding.FFTShift(1)
                    
                    if ~MR.Parameter.Encoding.WorkEncoding.FFTShift(1) || ~MR.Parameter.Encoding.WorkEncoding.FFTDims(2)
                        xrange = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.XRange, 'UniformOutput',0);
                    else
                        xrange = MR.Parameter.Encoding.WorkEncoding.XRange;
                        MR.Parameter.Encoding.WorkEncoding.FFTShift(1) = 0;
                    end
                    yrange = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.YRange, 'UniformOutput',0);
                    zrange = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.ZRange, 'UniformOutput',0);
                    zres = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.ZRes, 'UniformOutput',0);
                    
                    xrange = xrange( 1:size(MR.Data,1), 1:size(MR.Data,2) );
                    yrange = yrange( 1:size(MR.Data,1), 1:size(MR.Data,2) );
                    zrange = zrange( 1:size(MR.Data,1), 1:size(MR.Data,2) );
                    zres = zres( 1:size(MR.Data,1), 1:size(MR.Data,2) );
                    
                    MR.Data = cellfun( @(x, xr, yr, zr, zre)FFT.shift_image( x, xr, yr, zr, zre) , ...
                        MR.Data, ...
                        xrange, ...
                        yrange, ...
                        zrange, ...
                        zres, 'UniformOutput', 0 );
                    
                    MR.Parameter.Encoding.WorkEncoding.FFTShift = [0, 1, 1].*MR.Parameter.Encoding.WorkEncoding.FFTShift;
                end
                
                MR.Parameter.UpdateImageInfo = 1;
                MR.Data = Helper.UnconvertCell( MR.Data );
                MR.Parameter.ReconFlags.isimspace = [1, MR.Parameter.ReconFlags.isimspace(2:3) ];
            end
            
        end
        function K2IP( MR )
            % K2IP: Performs a fourier transformation from k-space to image-space in phase encoding
            % directions
            %
            % Syntax:     r.K2IP;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('EncodingPars.FFTDims')">Parameter.Encoding.FFTDims</a>:
            %             Specifies if the FFT in phase encoding directions should be executed. This
            %             is a 3 elements vector specifying the dimension in matrix notation.
            %           - <a href="matlab:helpwin('ReconFlags.isimspace')">Parameter.ReconFlags.isimspace</a>:
            %             Used to check if the phase encoding directions are still in k-space
            %           - <a href="matlab:helpwin('EncodingPars.FFTShift')">Parameter.Encoding.FFTShift</a>:
            %             Specifies if the image should be shifted in phase encoding directions after
            %             the fourier transform. This is a 3 elements vector specifying the
            %             shift in matrix notation.
            %           - <a href="matlab:helpwin('EncodingPars.XRange')">Parameter.Encoding.YRange</a>, <a href="matlab:helpwin('EncodingPars.ZRange')">Parameter.Encoding.ZRange</a>:
            %             Define the shifts of the image after the fourier transform in phase
            %             encoding directions.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlags.isimspace')">Parameter.ReconFlags.isimspace</a>
            %
            % Location:   k-space
            %
            % Formats:    Raw | ExportedRaw | Cpx | Rec | Bruker
            %
            % Description/Algorithm: K2IP first checks if the phase encoding directions should be ,
            %             transformed, defined in <a href="matlab:helpwin('EncodingPars.FFTDims')">Parameter.Encoding.FFTDims</a>. Afterwards it
            %             checks if the phase encoding directions are still in k-space and if yes
            %             transforms them. The fourier transformation along one dimension is given by:
            %
            %                   img = sqrt(n).* fftshift( ifft(ifftshift(kspace ,dim),[],dim), dim);
            %
            %             where n are the number of samples along the specific dimension.
            %
            %             After the fourier transform from k- to image-space, the images are shifted
            %             in phase encoding directions as defined in Parameter.Encoding.Y/ZRange.
            %             Thereby the images are shifted by the amount which results in a symmetric
            %             range around 0. E.g. if the range is given as [-42, 21] then the image is
            %             shifted by 10 pixel to the right. ( [-42, 21] + 10 = [-32, 31] )
            
            if all( MR.Parameter.ReconFlags.isimspace(2:3) )
                error( 'The data has alredy been fourier transformed in phase encoding direction');
            end
            
            if ~isempty(MR.Data)
                
                fftdims = find( MR.Parameter.Encoding.WorkEncoding.FFTDims & ...
                    MR.Parameter.ReconFlags.isimspace == 0);
                fftdims( fftdims == 1) = [];
                
                MR.Parameter.UpdateImageInfo = 0;
                MR.DataClass.Convert2Cell;
                for i = 1:length( fftdims)
                    % spectro begin ----------------------------
                    if MR.isSpectro
                        MR.Data = cellfun( @(x)FFT.k2i(x, fftdims(i), 1, true), MR.Data, 'UniformOutput', false );
                    else
                        % spectro end ------------------------------
                        MR.Data = cellfun( @(x)FFT.k2i(x, fftdims(i), 1), MR.Data, 'UniformOutput', false );
                    end
                end
                
                if any( MR.Parameter.Encoding.WorkEncoding.FFTShift(2:3) )
                    
                    xrange = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.XRange, 'UniformOutput',0);
                    
                    if ~MR.Parameter.Encoding.WorkEncoding.FFTShift(2) || ~MR.Parameter.Encoding.WorkEncoding.FFTDims(2)
                        yrange = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.YRange, 'UniformOutput',0);
                    else
                        yrange = MR.Parameter.Encoding.WorkEncoding.YRange;
                        MR.Parameter.Encoding.WorkEncoding.FFTShift(2) = 0;
                    end
                    if ~MR.Parameter.Encoding.WorkEncoding.FFTShift(3) || ~MR.Parameter.Encoding.WorkEncoding.FFTDims(3)
                        zrange = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.ZRange, 'UniformOutput',0);
                        zres = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.ZRes, 'UniformOutput',0);
                    else
                        zrange = MR.Parameter.Encoding.WorkEncoding.ZRange;
                        zres = MR.Parameter.Encoding.WorkEncoding.ZRes;
                        MR.Parameter.Encoding.WorkEncoding.FFTShift(3) = 0;
                    end
                    
                    xrange = xrange( 1:size(MR.Data,1), 1:size(MR.Data,2) );
                    yrange = yrange( 1:size(MR.Data,1), 1:size(MR.Data,2) );
                    zrange = zrange( 1:size(MR.Data,1), 1:size(MR.Data,2) );
                    zres = zres( 1:size(MR.Data,1), 1:size(MR.Data,2) );
                    
                    MR.Data = cellfun( @(x, xr, yr, zr, zre)FFT.shift_image( x, xr, yr, zr, zre) , ...
                        MR.Data, ...
                        xrange, ...
                        yrange, ...
                        zrange, ...
                        zres, 'UniformOutput', 0 );
                    MR.Parameter.Encoding.WorkEncoding.FFTShift = [1, 0, 0].*MR.Parameter.Encoding.WorkEncoding.FFTShift;
                end
                
                
                MR.Parameter.UpdateImageInfo = 1;
                MR.Data = Helper.UnconvertCell( MR.Data );
                
                MR.Parameter.ReconFlags.isimspace = [MR.Parameter.ReconFlags.isimspace(1),1,1 ];
            end
        end
        function I2K( MR )
            % I2K: Performs a fourier transformation from image-space to k-space
            %
            % Syntax:     r.I2K;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('EncodingPars.FFTDims')">Parameter.Encoding.FFTDims</a>:
            %             Specifies in which dimension the fourier transformation should be
            %             executed. This is a 3 elements vector specifying the dimension in matrix
            %             notation.
            %           - <a href="matlab:helpwin('ReconFlags.isimspace')">Parameter.ReconFlags.isimspace</a>:
            %             Used to check which dimensions are already in image space
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlags.isimspace')">Parameter.ReconFlags.isimspace</a>
            %
            % Location:   Image-Space
            %
            % Formats:    Raw | ExportedRaw | Cpx | Rec | Bruker
            %
            % Description/Algorithm: I2K first checks which dimensions should be transformed,
            %             defined in <a href="matlab:helpwin('EncodingPars.FFTDims')">Parameter.Encoding.FFTDims</a>. Afterwards it checks which of
            %             these dimensions are in image-space and then only transforms these ones.
            %             The fourier transformation along one dimension is given by:
            %
            %                   kspace = 1/sqrt(n)*fftshift(fft(ifftshift( img, dim ),[],dim),dim);
            %
            %             where n are the number of samples along the specific dimension.
            
            if all( ~MR.Parameter.ReconFlags.isimspace )
                error( 'Error in I2K: The data is already in k-space');
            end
            
            fftdims = find( MR.Parameter.Encoding.WorkEncoding.FFTDims & ...
                MR.Parameter.ReconFlags.isimspace == 1);
            
            if ~isempty(MR.Data)
                MR.Parameter.UpdateImageInfo = 0;
                MR.DataClass.Convert2Cell;
                
                % spectro begin ----------------------------
                if MR.isSpectro
                    MR.Data = cellfun( @(x)FFT.i2k(x, fftdims, 1, true ), MR.Data, 'UniformOutput', false );
                else
                    % spectro end ------------------------------
                    MR.Data = cellfun( @(x)FFT.i2k(x, fftdims, 1 ), MR.Data, 'UniformOutput', false );
                end
                
                MR.Parameter.UpdateImageInfo = 1;
                MR.Data = Helper.UnconvertCell( MR.Data );
                
                MR.Parameter.ReconFlags.isimspace(fftdims) = 0;
            end
        end
        
        % ---------------------------------------------------------------%
        % Coil Combination
        % ---------------------------------------------------------------%
        function CombineCoils( MR )
            % CombineCoils: Combines the individual coil images into one single image. Several
            % combination methods have been implemented in MRecon (see below).
            %
            % Syntax:     r.CombineCoils;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.CoilCombination')">Parameter.Recon.CoilCombination</a>: {'sos'} | 'pc' | 'svd' | 'snr_weight'
            %             When set on 'sos' a sum-of-squares combination is performed. Note that the image
            %             phase is lost with this method. 'pc' is used in phase contrast flow
            %             reconstructions where the magnitude is a sum-of-squares combination and the
            %             phase is a magnitude weighted combination of the individual phase images. 'svd' and 'snr_weight'
            %             are intended for spectroscopy data. For spectro 'svd' is the default because it works
            %             without prior phasing the individual coil signals. For 'snr_weight' the individual channels
            %             should already be correctly phased either by external phase correction or by the
            %             implemented eddy current correction procedure.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.iscombined')">Parameter.ReconFlags.iscombined</a>
            %
            % Location:   Image-space.
            %
            % Formats:    Raw | ExportedRaw | Cpx | Bruker
            %
            % Notes:    - For a non-phase-contrast scan the image phase is lost after CombineCoils. If
            %             you want to preserve the phase please perform a SENSE/CLEAR reconstruction see
            %             the SENSEUnfold function. This does not apply to spectroscopy data where the phase
            %             is always retained.
            %           - The CombineCoils function destroys the label lookup
            %             table (Parameter.LabelLookupTable) since the profiles cannot be related to the
            %             labels after combining the coils.
            %           - The coils have to be in the 4th dimension of the Data
            %             array to be combined.
            
            if ~isempty(MR.Data)
                if ~strcmpi( MR.Parameter.Recon.CoilCombination, 'no' )
                    if MR.Parameter.ReconFlags.iscombined
                        error( 'Error in CombineCoils: Data is already combined');
                    end
                    
                    MR.Parameter.UpdateImageInfo = 0;
                    MR.DataClass.Convert2Cell;
                    
                    % spectro begin ----------------------------
                    if MR.isSpectro
                        if ( size(MR.Data{1}, MRecon.dim.coil) > 1 )
                            if ~isempty(MR.Parameter.Encoding.YRes) && ~isempty(MR.Parameter.Encoding.ZRes)
                                error('SVD coil combination works at the moment only on single voxel data');
                            end
                        end
                        if ~MR.Parameter.ReconFlags.issorted
                            error('The Data should be sorted prior to coil combination');
                        end
                    else
                        if any( ~MR.Parameter.ReconFlags.isimspace )
                            error( 'Error in CombineCoils: Coil combination has to be performed on image space data');
                        end
                    end
                    % spectro end ------------------------------
                    
                    is_multichannel = cellfun( @(x) size(x, 4) ~= 1, MR.Data);
                    switch lower(MR.Parameter.Recon.CoilCombination)
                        case 'sos'
                            MR.Data(is_multichannel) = cellfun( @(x)sos(x, 4, 1), MR.Data(is_multichannel), 'UniformOutput', false );                          
                        case 'pc'
                            MR.Data(is_multichannel) = cellfun( @(x)sos(x, 4, 2), MR.Data(is_multichannel), 'UniformOutput', false );                           
                            % spectro begin ----------------------------
                        case 'svd'
                            if ~MR.isSpectro
                                error('SVD based coil combination is only available for spectroscopy scans');
                            end
                            MR.MRImpl.SpectroCombineCoils('svd');
                        case 'snr-weight'
                            if ~MR.isSpectro
                                error('SNR weighted coil combination is only available for spectroscopy scans');
                            end
                            if ~MR.Parameter.ReconFlags.isecc
                                error('SNR weighted coil combination can only work properly with phased data. Use eddy current correction before');
                            end
                            MR.MRImpl.SpectroCombineCoils('snr-weight');
                            % spectro end ------------------------------
                            
                    end
                    
                    MR.Parameter.UpdateImageInfo = 1;
                    
                    % Delete label lookup table
                    MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                    MR.Parameter.LabelLookupTable = cellfun( @(x)[], MR.Data , 'UniformOutput',0);
                    MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                    
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    
                    MR.Parameter.ReconFlags.iscombined = 1;
                end
                
            end
        end
        
        % ---------------------------------------------------------------%
        % SENSE / Partial Fourier
        % ---------------------------------------------------------------%
        function PartialFourier( MR )
            % PartialFourier: Performs a partial fourier reconstruction (homodyne method).
            %
            % Syntax:     r.PartialFourier;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.PartialFourier')">Parameter.Recon.PartialFourier</a>: {'yes'} | 'no'
            %             Enables/disables the partial fourier reconstruction.
            %           - <a href="matlab:helpwin('MRparameterDoc.DataFormat')">Parameter.DataFormat</a>:
            %             Used to check if the right data format is passed to the function
            %           - <a href="matlab:helpwin('ReconFlagParsDoc.isimspace')">Parameter.ReconFlags.isimspace</a>:
            %             Used to check if the data is in k- or image-space.
            %           - <a href="matlab:helpwin('EncodingPars.KxRange')">Parameter.Encoding.KxRange</a>, <a href="matlab:helpwin('EncodingPars.KyRange')">Parameter.Encoding.KyRange</a>, <a href="matlab:helpwin('EncodingPars.KzRange')">Parameter.Encoding.KzRange</a>:
            %             Specifies the sampling ranges in k-space. These parameters are used to
            %             check if the scan is acquired using partial fourier and define the fully
            %             sampled k-space center.
            %           - <a href="matlab:helpwin('ScanPars.SENSEFactor')">Parameter.Scan.SENSEFactor</a>:
            %             The SENSE factor used to define the fully sampled k-space center.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.ispartialfourier')">Parameter.ReconFlags.ispartialfourier</a>
            %
            % Location:   k-space | image-space
            %
            % Formats:    Raw | ExportedRaw | Bruker
            %
            % Description/Algorithm: In partial fourier imaging only a fraction of k-space is
            %             acquired for improved scan- or echo-time. If the k-space is partly sampled
            %             in readout direction we speak of partial echo and halfscan for a partly
            %             sampled space in phase encoding direction. In the first case the echo time
            %             is reduced while in the second one the scantime is shortened. Performing a
            %             standard imaging reconstruction on the partly sampled k-space however
            %             leads to image bluring and reduced SNR. Therefore a partial fourier
            %             reconstruction is performed which utilizes the k-space symmetry to
            %             reconstruct images of better quality. In MRecon a so called homodyne
            %             reconstruction is implemented: Noll, D.C.; Nishimura, D.G.; Macovski, A.;
            %             "Homodyne detection in magnetic resonance imaging" Medical Imaging,
            %             IEEE Transactions on, vol.10, no.2, pp.154-163, Jun 1991.
            %
            % Notes:    - In principle the partial fourier reconstruction can either be performed in
            %             k- or image space. However if executed in k-space it is not possible to
            %             perform a SENSE or CLEAR reconstruction afterwards since the phase
            %             information will be lost. Therefore we suggest to execute it in image
            %             space after SENSEUnfold which works in any case.
            
            if strcmpi( MR.Parameter.Recon.PartialFourier, 'yes' )
                if ~MR.Parameter.ReconFlags.issorted
                    error( 'Error in partial Fourier: Please sort the data first' );
                end
                
                if ~isempty(MR.Data)
                    is_kspace  = sum( MR.Parameter.ReconFlags.isimspace ) == 0;
                    is_imspace = sum( MR.Parameter.ReconFlags.isimspace ) == 3;
                    
                    % The data has to be fully in k-space or image
                    % space
                    if ~is_imspace && ~is_kspace
                        error( 'Error in partial fourier: The data has to be fully in k-space or image-space' );
                    end
                    
                    % If a SENSE recon is performed then partial
                    % fourier has to be executed in image space
                    if is_kspace && ...
                            ( ~isempty( MR.Parameter.Scan.SENSEFactor )   && ...
                            ~isempty( MR.Parameter.Recon.Sensitivities) && ...
                            strcmpi( MR.Parameter.Recon.SENSE, 'yes') )
                        error( 'Error in partial fourier: The partial fourier reconstruction has to be performed in image space if SENSE is enabled' );
                    end
                    
                    if is_imspace && ~isempty( MR.Parameter.Scan.SENSEFactor)
                        if MR.Parameter.ReconFlags.isrotated
                            error( 'The partial Fourier reconstruction cannot be performed on rotated images. Please call the RotateImage function after SENSEUnfold' );
                        end
                        m_ovs = cellfun( @(x)1/MR.Parameter.Scan.SENSEFactor(1), MR.Parameter.Encoding.WorkEncoding.KxRange, 'UniformOutput', 0);
                        MR.Parameter.Encoding.WorkEncoding.KxRange = cellfun( @(x,y) Helper.set_k_ranges( x, y ), MR.Parameter.Encoding.WorkEncoding.KxRange, ...
                            m_ovs, 'UniformOutput',0 );
                        p_ovs = cellfun( @(x)1/MR.Parameter.Scan.SENSEFactor(2), MR.Parameter.Encoding.WorkEncoding.KyRange, 'UniformOutput', 0);
                        MR.Parameter.Encoding.WorkEncoding.KyRange = cellfun( @(x,y) Helper.set_k_ranges( x, y ), MR.Parameter.Encoding.WorkEncoding.KyRange, ...
                            p_ovs, 'UniformOutput',0 );
                        s_ovs = cellfun( @(x)1/MR.Parameter.Scan.SENSEFactor(3), MR.Parameter.Encoding.WorkEncoding.KzRange, 'UniformOutput', 0);
                        MR.Parameter.Encoding.WorkEncoding.KzRange = cellfun( @(x,y) Helper.set_k_ranges( x, y ), MR.Parameter.Encoding.WorkEncoding.KzRange, ...
                            s_ovs, 'UniformOutput',0 );
                    end
                    
                    if any( strcmpi( MR.Parameter.DataFormat, {'Raw', 'ExportedRaw', 'Bruker'} ) )
                        if all( MR.Parameter.ReconFlags.ispartialfourier )
                            error( 'Error in partial fourier: The partial fourier correction has already been applied' );
                        end
                        
                        MR.Parameter.UpdateImageInfo = 0;
                        MR.DataClass.Convert2Cell;
                        
                        % Check if we have to perform a partial fourier
                        % recon at all
                        isPF = logical(cellfun( @(x,kx, ky, kz)PartialFourier.is_partial_fourier( x, kx, ky, kz ), ...
                            MR.Data(1,:), MR.Parameter.Encoding.WorkEncoding.KxRange(1,:), ...
                            MR.Parameter.Encoding.WorkEncoding.KyRange(1,:), ...
                            MR.Parameter.Encoding.WorkEncoding.KzRange(1,:)));
                        
                        if any(isPF)
                            apply_filter = ~MR.Parameter.ReconFlags.ispartialfourier(1);
                            
                            % Do the forward transform here to save memory (crappy
                            % Matlab memory handling)
                            
                            % The homodyne filter hasn't been applied yet
                            
                            if is_imspace
                                orig_phase = cellfun( @(x)angle(x), MR.Data(1,isPF) , 'UniformOutput', 0 );
                                if apply_filter
                                    % Do the fft individually for each dimension to save
                                    % memory
                                    MR.I2K;
                                    MR.DataClass.Convert2Cell;                                   
                                end
                            end
                            
                            if apply_filter
                                [MR.Data(1,isPF), MR.MRImpl.PFLowRes] = cellfun( @(d, x, y, z)PartialFourier.partial_fourier_filter( d, x, y, z ), ...
                                    MR.Data(1,isPF), MR.Parameter.Encoding.WorkEncoding.KxRange(1,isPF), ...
                                    MR.Parameter.Encoding.WorkEncoding.KyRange(1,isPF), ...
                                    MR.Parameter.Encoding.WorkEncoding.KzRange(1,isPF), 'UniformOutput', false );
                                
                                if is_kspace
                                    MR.MRImpl.PFLowRes = cellfun( @(x, xr, yr, zr, zre)FFT.shift_image( x, xr, yr, zr, zre) , ...
                                        MR.MRImpl.PFLowRes, ...
                                        MR.Parameter.Encoding.WorkEncoding.XRange(1,isPF), ...
                                        MR.Parameter.Encoding.WorkEncoding.YRange(1,isPF), ...
                                        MR.Parameter.Encoding.WorkEncoding.ZRange(1,isPF), ...
                                        MR.Parameter.Encoding.WorkEncoding.ZRes(1,isPF), 'UniformOutput', 0 );
                                end
                            end
                            
                            if any(isPF)
                                if apply_filter || is_kspace
                                    % Do the backward transform here to save memory (crappy
                                    % Matlab memory handling)
                                    MR.K2I;
                                    MR.DataClass.Convert2Cell;                                   
                                end
                                
                                MR.Data(1,isPF) = cellfun( @(x,y)PartialFourier.partial_fourier_multiply( x, y ), MR.Data(1,isPF), MR.MRImpl.PFLowRes, 'UniformOutput', 0);
                                
                                if is_kspace                                   
                                    MR.I2K;
                                    MR.DataClass.Convert2Cell;
                                else
                                    MR.Data(1,isPF) = cellfun( @(x, y)abs(x).* (cos(y) + 1i.*sin(y)) , MR.Data(1,isPF), orig_phase, 'UniformOutput', 0 );
                                end
                            end
                        end
                        MR.Parameter.UpdateImageInfo = 1;
                        MR.Data = Helper.UnconvertCell( MR.Data );
                    end
                end
                MR.Parameter.ReconFlags.ispartialfourier = [1, 1];
            end
        end
        function PartialFourierFilter( MR )
            % PartialFourierFilter: Applies the homodyne filter in k-space
            %
            % Syntax:     r.PartialFourierFilter;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.PartialFourier')">Parameter.Recon.PartialFourier</a>: {'yes'} | 'no'
            %             Enables/disables the partial fourier reconstruction.
            %           - <a href="matlab:helpwin('MRparameterDoc.DataFormat')">Parameter.DataFormat</a>:
            %             Used to check if the right data format is passed to the function
            %           - <a href="matlab:helpwin('ReconFlagParsDoc.isimspace')">Parameter.ReconFlags.isimspace</a>:
            %             Used to check if the data is in k- or image-space.
            %           - <a href="matlab:helpwin('EncodingPars.KxRange')">Parameter.Encoding.KxRange</a>, <a href="matlab:helpwin('EncodingPars.KyRange')">Parameter.Encoding.KyRange</a>, <a href="matlab:helpwin('EncodingPars.KzRange')">Parameter.Encoding.KzRange</a>:
            %             Specifies the sampling ranges in k-space. These parameters are used to
            %             check if the scan is acquired using partial fourier and define the fully
            %             sampled k-space center.
            %           - <a href="matlab:helpwin('ScanPars.SENSEFactor')">Parameter.Scan.SENSEFactor</a>:
            %             The SENSE factor used to define the fully sampled k-space center.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.ispartialfourier')">Parameter.ReconFlags.ispartialfourier</a>
            %
            % Location:   k-space
            %
            % Formats:    Raw | ExportedRaw | Bruker
            %
            % Description/Algorithm: This function applies the homodyne filter in k-space, used for
            %             the partial fourier reconstruction. If this function is called in k-space
            %             before K2I it saves the need for two extra FFT's when PartialFourier is
            %             executed and thus increases the performance of the recon. Calling this
            %             function is optional. If it has not been called the filter will be applied
            %             during PartialFourier. For more information about the partial fourier
            %             reconstruction refer to the documentation of PartialFourier.
            
            if strcmpi( MR.Parameter.Recon.PartialFourier, 'yes' )
                if ~isempty(MR.Data)
                    
                    % Do not apply the filter if SENSE/CLEAR is enabled (with SENSE the coils are already
                    % combined when the partial fourier recon is applied. Here the coils are not
                    % combined yet.
                    if isempty( MR.Parameter.Recon.Sensitivities )
                        
                        if ~MR.Parameter.ReconFlags.issorted
                            error( 'Error in PartialFourierFilter: Please sort the data first' );
                        end
                        
                        is_kspace  = sum( MR.Parameter.ReconFlags.isimspace ) == 0;
                        if ~is_kspace
                            error( 'Error in PartialFourierFilter: The homodyne filter has to be applied in k-space' );
                        end
                        
                        if any( strcmpi( MR.Parameter.DataFormat, {'Raw', 'ExportedRaw', 'Bruker'} ) )
                            if MR.Parameter.ReconFlags.ispartialfourier(1)
                                error( 'Error in PartialFourierFilter: The homodyne filter has already been applied' );
                            end
                            
                            MR.Parameter.UpdateImageInfo = 0;
                            MR.DataClass.Convert2Cell;
                            
                            % Check if we have to perform a partial fourier
                            % recon at all
                            isPF = logical(cellfun( @(x,kx, ky, kz)PartialFourier.is_partial_fourier( x, kx, ky, kz ), ...
                                MR.Data(1,:), MR.Parameter.Encoding.WorkEncoding.KxRange(1,:), ...
                                MR.Parameter.Encoding.WorkEncoding.KyRange(1,:), ...
                                MR.Parameter.Encoding.WorkEncoding.KzRange(1,:)));
                            
                            if any(isPF)
                                [MR.Data(1,isPF), MR.MRImpl.PFLowRes] = cellfun( @(d, x, y, z)PartialFourier.partial_fourier_filter( d, x, y, z), ...
                                    MR.Data(1,isPF), MR.Parameter.Encoding.WorkEncoding.KxRange(1,isPF), ...
                                    MR.Parameter.Encoding.WorkEncoding.KyRange(1,isPF), ...
                                    MR.Parameter.Encoding.WorkEncoding.KzRange(1,isPF), ...
                                    'UniformOutput', false );
                                
                                MR.MRImpl.PFLowRes = cellfun( @(x, xr, yr, zr, zre)FFT.shift_image( x, xr, yr, zr, zre) , ...
                                    MR.MRImpl.PFLowRes, ...
                                    MR.Parameter.Encoding.WorkEncoding.XRange(1,isPF), ...
                                    MR.Parameter.Encoding.WorkEncoding.YRange(1,isPF), ...
                                    MR.Parameter.Encoding.WorkEncoding.ZRange(1,isPF), ...
                                    MR.Parameter.Encoding.WorkEncoding.ZRes(1,isPF), 'UniformOutput', 0 );
                            end
                            
                            MR.Parameter.UpdateImageInfo = 1;
                            MR.Data = Helper.UnconvertCell( MR.Data );
                            
                        end
                        MR.Parameter.ReconFlags.ispartialfourier(1) = 1;
                    end
                end
            end
        end
        function SENSEUnfold( MR )
            % SENSEUnfold: Performs a SENSE reconstruction using the coil sensitivity information
            % provided by the user.
            %
            % Syntax:     r.SENSEUnfold;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.SENSE')">Parameter.Recon.SENSE</a>: {'yes'} | 'no'
            %             Enables/disables the SENSE unfolding.
            %           - <a href="matlab:helpwin('ReconParsDoc.Sensitivities')">Parameter.Recon.Sensitivities</a>:
            %             The sensitivity object (MRsense) used in the unfolding process. The
            %             sensitivity object holds the coil sensitivities, the noise covariance
            %             matrix as well as regularisation images.
            %           - <a href="matlab:helpwin('ReconParsDoc.SENSERegStrength')">Parameter.Recon.SENSERegStrength</a>:
            %             The regularization strength used during the SENSE unfolding. A higher
            %             regularizazion results in a stronger suppression of low signal values in
            %             the unfolded images.
            %           - <a href="matlab:helpwin('ScanParsDoc.SENSEFactor')">Parameter.Scan.SENSEFactor</a>:
            %             The SENSE reduction factor.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.isunfolded')">Parameter.ReconFlags.isunfolded</a>
            %
            % Location:   image-space
            %
            % Formats:    Raw | ExportedRaw
            %
            % Description/Algorithm: Please see the <a href="matlab:open('MRsense quickstart.pdf')">MRsense manual</a> for more infomation on the
            % implemented algorithms and examples.
            
            if ~isempty( MR.Data ) && ...
                    ~isempty( MR.Parameter.Scan.SENSEFactor ) && ...
                    strcmpi( MR.Parameter.Recon.SENSE, 'yes' ) && ...
                    ~isempty( MR.Parameter.Recon.Sensitivities )
                
                if isempty( MR.Parameter.Recon.Sensitivities.Sensitivity )
                    error( 'Error in SENSEUnfold: No Sensitivities available. Run the MRsense Perform function' );
                end
                if any( ~MR.Parameter.ReconFlags.isimspace )
                    error( 'Error in SENSEUnfold: SENSE Unfolding has to be performed on image space data');
                end
                if MR.Parameter.ReconFlags.isunfolded
                    error( 'Error in SENSEUnfold: Data is already Unfolded');
                end
                if isempty( MR.Parameter.Recon.Sensitivities )
                    error( 'Error in SENSEUnfold: Coil sensitivities not found! Please add sensitivity infomation in Parameter.Recon.Sensitivities...');
                end
                if isempty( MR.Parameter.Scan.SENSEFactor )
                    error( 'Error in SENSEUnfold: SENSE factors not found! Please add them in Parameter.Scan.SENSEFactor...');
                end
                for i = 1:3
                    if MR.Parameter.ReconFlags.isoversampled(i) == 0 && MR.Parameter.Scan.SENSEFactor(i) > 1
                        error( 'Error in SENSEUnfold: The Oversampling cannot be removed in the undersampled direction before the unfolding process');
                    end
                end
                if MR.Parameter.ReconFlags.iszerofilled(2)
                    error( 'Error in SENSEUnfold: The unfodling has to be performed before image space zero filling');
                end
                if MR.Parameter.ReconFlags.iscombined
                    error( 'Error in SENSEUnfold: The unfodling has to be performed before coil combination');
                end
                
                % The radon and iradon function is used in the fiex mex file. Add a call here (which is never executed)
                % so that matlab knows that it will be used and includes it in the exe file during the compilation
                if (0)
                    MR.Data = radon( MR.Data );
                    MR.Data = iradon( MR.Data, 0:179 );
                    fiex_config;
                end
                
                if ~isempty(MR.Data)
                    s_psi = [];
                    
                    MR.Parameter.UpdateImageInfo = 0;
                    MR.DataClass.Convert2Cell;
                    
                    res = cellfun( @(x)[round(size(x,1)*max([1,MR.Parameter.Scan.SENSEFactor(1,1)])),round(size(x,2)*max([1,MR.Parameter.Scan.SENSEFactor(1,2)])),round(size(x,3)*max([1,MR.Parameter.Scan.SENSEFactor(1,3)]))], ...
                        MR.Data, 'UniformOutput',0 );
                    
                    if( MR.Parameter.ReconFlags.iszerofilled(1) )
                        MR.Parameter.Encoding.WorkEncoding.KxOversampling = cellfun( @(x, y)x(1)/max( [1,y] ), ...
                            res, MR.Parameter.Encoding.WorkEncoding.XRes, 'UniformOutput',0 );
                        MR.Parameter.Encoding.WorkEncoding.KyOversampling = cellfun( @(x, y)x(2)/max( [1,y] ), ...
                            res, MR.Parameter.Encoding.WorkEncoding.YRes, 'UniformOutput',0 );
                        MR.Parameter.Encoding.WorkEncoding.KzOversampling = cellfun( @(x, y)x(3)/max( [1,y] ), ...
                            res, MR.Parameter.Encoding.WorkEncoding.ZRes, 'UniformOutput',0 );
                    end
                    
                    if ~isempty(MR.Parameter.Recon.Sensitivities.ExtraOversampling)
                        extra_ovs = MR.Parameter.Recon.Sensitivities.ExtraOversampling;
                        res = cellfun( @(x)[round(x(1)*max([1,extra_ovs(1)])),round(x(2)*max([1,extra_ovs(2)])),round(x(3)*max([1,extra_ovs(3)]))], ...
                            res, 'UniformOutput',0 );
                    end
                    
                    body_ref = [];
                    coil_ref = [];
                    mc = metaclass(MR.Parameter.Recon.Sensitivities);
                    if strcmpi( mc.Name, 'MRsense' ) || strcmpi( mc.Name, 'MRsenseSin' ) || strcmpi( mc.Name, 'MRsenseCpx' )
                        sens_image = MR.Parameter.Recon.Sensitivities.Sensitivity;
                        body_ref   = MR.Parameter.Recon.Sensitivities.ReformatedBodycoilData;
                        coil_ref   = MR.Parameter.Recon.Sensitivities.ReformatedCoilData;
                        s_psi = MR.Parameter.Recon.Sensitivities.Psi;
                    elseif isnumeric(MR.Parameter.Recon.Sensitivities)
                        sens_image = MR.Parameter.Recon.Sensitivities;
                    else
                        error( 'Unknown input type for the reference scan' );
                    end
                    
                    % set the sensitivity maps of channels which were not
                    % measured to 0 (per location)
                    try
                        for i = 1:size(coil_ref, 8 )
                            cur_stack = MR.Parameter.Labels.StackIndex(i)+1;
                            channel_numbers_measured = MR.Parameter.Labels.CoilNrsPerStack{cur_stack};
                            set2zero =  find( ismember( MR.Parameter.Recon.Sensitivities.ChannelNumbers, channel_numbers_measured ) == 0);
                            sens_image( :,:,:,set2zero, :,:,:,i,:,:,:,:,:) = 0;
                            coil_ref( :,:,:,set2zero, :,:,:,i,:,:,:,:,:) = 0;
                        end
                    catch
                    end
                    
                    % match the channel numbers in the sensitivity maps
                    % with the ones in the data (if possible)
                    try
                        P = MR.Parameter.Parameter2Read.Copy;
                        P.Update(MR.Parameter.Labels.Index);
                        [~,chan_ind] = ismember(P.chan, MR.Parameter.Recon.Sensitivities.ChannelNumbers );
                        sens_image = sens_image(:,:,:,chan_ind,:,:,:,:,:,:,:,:);
                        coil_ref = coil_ref(:,:,:,chan_ind,:,:,:,:,:,:,:,:);
                        s_psi = s_psi( chan_ind, chan_ind);
                    catch
                    end
                    
                    % check if the number of locations in the data matches
                    % the one of the sensitivities
                    nr_locas = max( cellfun( @(x)size(x, 8), MR.Data(1,:) ));
                    if nr_locas > size( sens_image, 8)
                        if MR.Parameter.ReconFlags.iscombined
                            error( 'Error in SENSEUnfold: The number of locations in the data exceeds the one in the sensitivities');
                        end
                    end
                    if nr_locas < size( sens_image, 8)
                        if length(MR.Parameter.Parameter2Read.loca) == nr_locas
                            sens_image = sens_image( :,:,:,:,:,:,:,MR.Parameter.Parameter2Read.loca+1,:,:,:,:);
                            coil_ref = coil_ref( :,:,:,:,:,:,:,MR.Parameter.Parameter2Read.loca+1,:,:,:,:);
                            body_ref = body_ref( :,:,:,:,:,:,:,MR.Parameter.Parameter2Read.loca+1,:,:,:,:);
                        else
                            error( 'Error in SENSEUnfold: The number of locations in the data and Sensitivity maps is different. Cannot determine which locations in the sensitivities belong to the data');
                        end
                    end
                                      
                    psi = [];
                    if isempty(psi)
                        if ~isempty( s_psi )
                            psi = single(s_psi);
                        else
                            psi = single( eye( size( sens_image,4) ) );
                        end
                    end
                    
                    % When array compression is switched on compress the sensitivities and the
                    % coil_ref
                    if strcmpi( MR.Parameter.Recon.ArrayCompression, 'yes' )
                        cur_nr_chans = size(sens_image, 4);
                        for ac_in = 1:size( MR.Parameter.Recon.ACMatrix, 3 )
                            row_ind = min( [size(MR.Parameter.Recon.ACMatrix, 1), find( isnan( MR.Parameter.Recon.ACMatrix(:,1,ac_in ) ),1 )-1]);
                            col_ind = min( [size(MR.Parameter.Recon.ACMatrix, 2), find( isnan( MR.Parameter.Recon.ACMatrix(1,:,ac_in ) ),1 )-1]);
                            A = MR.Parameter.Recon.ACMatrix(1:row_ind, 1:col_ind,ac_in);
                            
                            if size(A,2) ~= cur_nr_chans
                                continue
                            end
                            sens_image = permute( sens_image, [4,1,2,3,5,6,7,8,9,10,11,12,13] );
                            coil_ref = permute( coil_ref, [4,1,2,3,5,6,7,8,9,10,11,12,13] );
                            siz = size(sens_image); siz(1) = size(A,1);
                            sens_image = A*sens_image(:,:);
                            coil_ref = A*coil_ref(:,:);
                            sens_image = reshape(sens_image, siz);
                            coil_ref = reshape(coil_ref, siz);
                            sens_image = permute( sens_image, [2,3,4,1,5,6,7,8,9,10,11,12,13] );
                            coil_ref = permute( coil_ref, [2,3,4,1,5,6,7,8,9,10,11,12,13] );
                            psi = A*psi*A';
                        end
                    end
                    
                    MR.Data = cellfun( @(x,y)SENSE.sense_recon(x, sens_image, psi, y, coil_ref, body_ref, MR.Parameter.Recon.SENSERegStrength), MR.Data(1,:), res(1,:), 'UniformOutput', false );
                    
                    % remove the extra oversampling
                    if ~isempty(MR.Parameter.Recon.Sensitivities.ExtraOversampling)
                        extra_ovs = MR.Parameter.Recon.Sensitivities.ExtraOversampling;
                        [MR.Data, MR.Parameter.LabelLookupTable] = cellfun( @(x)ImageProduction.rem_ovs( x, [], extra_ovs(1), extra_ovs(2), extra_ovs(3) ), ...
                            MR.Data,...
                            'UniformOutput', 0 );
                    end
                    
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    MR.Parameter.ReconFlags.isunfolded = 1;
                    
                end
                
            end
        end
        
        % ---------------------------------------------------------------%
        % Image Production
        % ---------------------------------------------------------------%
        function ZeroFill( MR )
            % ZeroFill: Zero pads the data either in k-space or image-space
            %
            % Syntax:     r.ZeroFill;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.kSpaceZeroFill')">Parameter.Recon.kSpaceZeroFill</a>: {'yes'} | 'no'
            %             Enables/disables the k-space zero filling.
            %           - <a href="matlab:helpwin('ReconParsDoc.ImageSpaceZeroFill')">Parameter.Recon.ImageSpaceZeroFill</a>: {'yes'} | 'no'
            %             Enables/disables the image-space zero filling.
            %           - <a href="matlab:helpwin('ReconFlagParsDoc.isimspace')">Parameters.ReconFlags.isimspace</a>:
            %             Used to check if the data is in k- or image-space.
            %           - <a href="matlab:helpwin('ReconFlagParsDoc.isoversampled')">Parameters.ReconFlags.isoversampled</a>:
            %             Used to check if and in which directions the data is oversampled.
            %           - <a href="matlab:helpwin('EncodingParsDoc.XRes')">Parameter.Encoding.XRes</a>, <a href="matlab:helpwin('EncodingParsDoc.YRes')">Parameter.Encoding.YRes</a>, <a href="matlab:helpwin('EncodingParsDoc.ZRes')">Parameter.Encoding.ZRes</a>:
            %             The data is zero filled to these dimension in k-space (times the remaining oversampling factor).
            %           - <a href="matlab:helpwin('EncodingParsDoc.XReconRes')">Parameter.Encoding.XReconRes</a>, <a href="matlab:helpwin('EncodingParsDoc.YReconRes')">Parameter.Encoding.YReconRes</a>, <a href="matlab:helpwin('EncodingParsDoc.ZReconRes')">Parameter.Encoding.ZReconRes</a>:
            %             The data is zero filled to these dimension in image-space (times the remaining oversampling factor).
            %           - <a href="matlab:helpwin('EncodingParsDoc.KxOversampling')">Parameter.Encoding.KxOversampling</a>, <a href="matlab:helpwin('EncodingParsDoc.KyOversampling')">Parameter.Encoding.KyOversampling</a>, <a href="matlab:helpwin('EncodingParsDoc.KzOversampling')">Parameter.Encoding.KzOversampling</a>:
            %             The oversamping factors in x/y/z-directions.
            %           - <a href="matlab:helpwin('ReconFlagParsDoc.issorted')">Parameters.ReconFlags.issorted</a>:
            %             Used to check if the data is already sorted.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.iszerofilled')">Parameter.ReconFlags.iszerofilled</a>
            %
            % Location:   k-space | Image-space.
            %
            % Formats:    Raw | ExportedRaw | Cpx | Bruker
            %
            % Description/Algorithm: The ZeroFill function zero-pads the current data to the
            %             dimensions specified in the Encoding parameters (see above). Thereby we
            %             have a different dimension depending whether we are in k-space or
            %             image-space. For example the dimension in y-direction after zero filling
            %             in k-space is given by:
            %
            %               yres = r.Parameter.Encoding.YRes * cur_ovs
            %
            %             where the current oversampling cur_ovs is given by:
            %
            %               cur_ovs = r.Parameter.Encoding.KyOversampling  (if oversampling is not removed)
            %             or
            %               cur_ovs = 1                                    (if oversampling is removed)
            %
            %             Whether and in which directions the data is oversampled is defined by the
            %             oversampling ReconFlag (<a href="matlab:helpwin('ReconFlagParsDoc.isoversampled')">Parameters.ReconFlags.isoversampled</a>).
            %
            %             Zero filling in k-space is used to correct for anisotropic voxel sizes and
            %             the artificially increase image resolution. Zero filling in k-space is
            %             used to correct for a rectangular field-of-view and to produce quadratic
            %             images.
            
            if ~MR.Parameter.ReconFlags.isread
                error( 'Error in ZeroFill: Please read the data first' );
            end
            if ~MR.Parameter.ReconFlags.issorted
                error( 'Error in ZeroFill: Please sort the data first' );
            end
            
            if ~isempty(MR.Data)
                MR.Parameter.UpdateImageInfo = 0;
                MR.DataClass.Convert2Cell;
                
                if isempty( MR.Parameter.LabelLookupTable )
                    MR.Parameter.LabelLookupTable = cellfun( @(x)[], MR.Data, 'UniformOutput', 0);
                else
                    MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                end
                
                if any( MR.Parameter.ReconFlags.isimspace )
                    if strcmpi( MR.Parameter.Recon.ImageSpaceZeroFill, 'yes' )
                        xres = cellfun( @(x,y) round(x*y), MR.Parameter.Encoding.WorkEncoding.XReconRes, MR.Parameter.Encoding.WorkEncoding.KxOversampling, 'UniformOutput', 0);
                        yres = cellfun( @(x,y) round(x*y), MR.Parameter.Encoding.WorkEncoding.YReconRes, MR.Parameter.Encoding.WorkEncoding.KyOversampling, 'UniformOutput', 0);
                        zres = cellfun( @(x,y) round(x*y), MR.Parameter.Encoding.WorkEncoding.ZReconRes, MR.Parameter.Encoding.WorkEncoding.KzOversampling, 'UniformOutput', 0);
                        
                        xres_temp = cell(size(xres));
                        yres_temp = cell(size(yres));
                        zres_temp = cell(size(zres));
                        for i = 1:size( xres, 2 )
                            res = [ xres(:,i), yres(:,i), zres(:,i) ];
                            mps = MR.Parameter.Encoding.WorkEncoding.MPS;
                            res= [ res(:, mps(1)), res(:, mps(2)), res(:, mps(3)) ];
                            xres_temp(:,i) = res(:,1);
                            yres_temp(:,i) = res(:,2);
                            zres_temp(:,i) = res(:,3);
                        end
                        xres = xres_temp;
                        yres = yres_temp;
                        zres = zres_temp;
                        
                        xres = xres(1:size(MR.Data,1), 1:size(MR.Data,2));
                        yres = yres(1:size(MR.Data,1), 1:size(MR.Data,2));
                        zres = zres(1:size(MR.Data,1), 1:size(MR.Data,2));
                        % check if the matrix in r.Data is larger than what it
                        % should be after after zero filling.
                        xres_check = find(cellfun( @(x,y) max( [0, size(x, 1) > y]), MR.Data, xres), 1);
                        yres_check = find(cellfun( @(x,y) max( [0, size(x, 2) > y]), MR.Data, yres), 1);
                        zres_check = find(cellfun( @(x,y) max( [0, size(x, 3) > y]), MR.Data, zres), 1);
                        wrong_dim = min( [xres_check, yres_check, zres_check]);
                        if ~isempty( wrong_dim )
                            error( 'Error in ZeroFill: The current Data matrix is larger than what it should be after zero filling!\nThe Data matrix size is %d x %d x %d and you want to zero fill it to %d x %d x %d\nCheck the values of XReconRes, YReconRes, ZReconRes in Parameter.Encoding', ...
                                size(MR.Data{wrong_dim}, 1), size(MR.Data{wrong_dim}, 2), size(MR.Data{wrong_dim}, 3), max([xres{wrong_dim},1]), max([yres{wrong_dim},1]), max([zres{wrong_dim},1]) );
                        end
                        
                        factor = [ max( [1, xres{1,1}./size(MR.Data{1,1}, 1)]), max( [1, yres{1,1}./size(MR.Data{1,1}, 2)]), max( [1, zres{1,1}./size(MR.Data{1,1}, 3)]) ];
                        MR.Parameter.UpdateCurFOV( factor );
                        
                        MR.Data = cellfun( @(x, xr, yr, zr, l)zero_fill( x, xr, yr, zr), MR.Data , xres, yres, zres, 'UniformOutput', 0);
                        if ~isempty(MR.Parameter.LabelLookupTable{1})
                            MR.Parameter.LabelLookupTable = cellfun( @( x, yr, zr, l)zero_fill( x, [], yr, zr), MR.Data , yres, zres, 'UniformOutput', 0);
                        end
                        
                        MR.Parameter.ReconFlags.iszerofilled(2) = 1;
                    end
                else
                    if strcmpi( MR.Parameter.Recon.kSpaceZeroFill, 'yes' )
                        xres = cellfun( @(x,y) round(x*y), MR.Parameter.Encoding.WorkEncoding.XRes, MR.Parameter.Encoding.WorkEncoding.KxOversampling, 'UniformOutput', 0);
                        yres = cellfun( @(x,y) round(x*y), MR.Parameter.Encoding.WorkEncoding.YRes, MR.Parameter.Encoding.WorkEncoding.KyOversampling, 'UniformOutput', 0);
                        zres = cellfun( @(x,y) round(x*y), MR.Parameter.Encoding.WorkEncoding.ZRes, MR.Parameter.Encoding.WorkEncoding.KzOversampling, 'UniformOutput', 0);
                        
                        if ~isempty( MR.Parameter.Scan.SENSEFactor ) && ~all( MR.Parameter.Scan.SENSEFactor(1,:) == 1 ) && ~MR.Parameter.ReconFlags.isunfolded
                            xres = cellfun( @(x) round(x/MR.Parameter.Scan.SENSEFactor(1, min([1, length(MR.Parameter.Scan.SENSEFactor)] ) )) , xres, 'UniformOutput', 0);
                            yres = cellfun( @(x) round(x/MR.Parameter.Scan.SENSEFactor(1, min([2, length(MR.Parameter.Scan.SENSEFactor)] ) )) , yres, 'UniformOutput', 0);
                            zres = cellfun( @(x) round(x/MR.Parameter.Scan.SENSEFactor(1, min([3, length(MR.Parameter.Scan.SENSEFactor)] ) )) , zres, 'UniformOutput', 0);
                        end
                        xres = xres(1:size(MR.Data,1), 1:size(MR.Data,2));
                        yres = yres(1:size(MR.Data,1), 1:size(MR.Data,2));
                        zres = zres(1:size(MR.Data,1), 1:size(MR.Data,2));
                        
                        xres_temp = cell(size(xres));
                        yres_temp = cell(size(yres));
                        zres_temp = cell(size(zres));
                        for i = 1:size( xres, 2 )
                            res = [ xres(:,i), yres(:,i), zres(:,i) ];
                            mps = MR.Parameter.Encoding.WorkEncoding.MPS;
                            res= [ res(:, mps(1)), res(:, mps(2)), res(:, mps(3)) ];
                            xres_temp(:,i) = res(:,1);
                            yres_temp(:,i) = res(:,2);
                            zres_temp(:,i) = res(:,3);
                        end
                        xres = xres_temp;
                        yres = yres_temp;
                        zres = zres_temp;
                        
                        % check if the matrix in r.Data is larger than what it
                        % should be after after zero filling.
                        xres_check = find(cellfun( @(x,y) max( [0, size(x, 1) > y]), MR.Data, xres), 1);
                        yres_check = find(cellfun( @(x,y) max( [0, size(x, 2) > y]), MR.Data, yres), 1);
                        zres_check = find(cellfun( @(x,y) max( [0, size(x, 3) > y]), MR.Data, zres), 1);
                        wrong_dim = min( [xres_check, yres_check, zres_check]);
                        if ~isempty( wrong_dim )
                            error( 'Error in ZeroFill: The current Data matrix is larger than what it should be after zero filling!\nThe Data matrix size is %d x %d x %d and you want to zero fill it to %d x %d x %d\nCheck the values of XRes, YRes, ZRes in Parameter.Encoding', ...
                                size(MR.Data{wrong_dim}, 1), size(MR.Data{wrong_dim}, 2), size(MR.Data{wrong_dim}, 3), max([xres{wrong_dim},1]), max([yres{wrong_dim},1]), max([zres{wrong_dim},1]) );
                        end
                        
                        % check if we have to zero fill at all
                        x_size_change = any(cellfun( @(x,y) size(x, 1) ~= max([y,1]), MR.Data, xres));
                        y_size_change = any(cellfun( @(x,y) size(x, 2) ~= max([y,1]), MR.Data, yres));
                        z_size_change = any(cellfun( @(x,y) size(x, 3) ~= max([y,1]), MR.Data, zres));
                        sizes_change = x_size_change | y_size_change | z_size_change;
                        
                        
                        if( sizes_change)                            
                            MR.Data = cellfun( @(x, xr, yr, zr, l)zero_fill( x, xr, yr, zr), MR.Data , xres, yres, zres, 'UniformOutput', 0);
                            MR.Parameter.LabelLookupTable = cellfun( @( x, yr, zr)zero_fill( x, [], yr, zr), MR.Parameter.LabelLookupTable , yres, zres, 'UniformOutput', 0);
                        end
                        MR.Parameter.ReconFlags.iszerofilled(1) = 1;
                    end
                end
                
                MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                MR.Parameter.UpdateImageInfo = 1;
                MR.Data = Helper.UnconvertCell( MR.Data );
            end
            
        end
        function RemoveOversampling( MR, varargin )
            % RemoveOversampling: Removes the oversampling from the image by cropping the data in
            % image-space. For spectroscopy an additional algorithm based on Linear Prediction (LP) is available
            % for the readout direction. It is performed in the time domain and leads to much reduced truncation
            % artifacts at the beginning and the end of the FID [Fuchs A et al., ISMRM 2012].
            %
            % Syntax:     r.RemoveOversampling;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.RemoveMOversampling')">Parameter.Recon.RemoveMOversampling</a>: {'yes'} | 'no'
            %             Enables/disables the oversampling removal in readout (measurement)
            %             direction
            %           - <a href="matlab:helpwin('ReconParsDoc.RemovePOversampling')">Parameter.Recon.RemovePOversampling</a>: {'yes'} | 'no'
            %             Enables/disables the oversampling removal in phase-encoding direction(s)
            %           - <a href="matlab:helpwin('ReconFlagParsDoc.isimspace')">Parameters.ReconFlags.isimspace</a>:
            %             Used to check if the data is in k- or image-space.
            %           - <a href="matlab:helpwin('EncodingParsDoc.KxOversampling')">Parameter.Encoding.KxOversampling</a>, <a href="matlab:helpwin('EncodingParsDoc.KyOversampling')">Parameter.Encoding.KyOversampling</a>, <a href="matlab:helpwin('EncodingParsDoc.KzOversampling')">Parameter.Encoding.KzOversampling</a>:
            %             The oversamping factors in x/y/z-directions. After calling the
            %             RemoveOversampling function the Data matrix is smaller by these factors.
            %           - <a href="matlab:helpwin('ReconFlagParsDoc.issorted')">Parameters.ReconFlags.issorted</a>:
            %             Used to check if the data is already sorted.
            %           - <a href="matlab:helpwin('SpectroPars.Downsample')">Parameter.Spectro.Downsample</a>:
            %             Enable/disable the LP based algorithm for spectroscopy.
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.isoversampled')">Parameter.ReconFlags.isoversampled</a>
            %
            % Location:   k-space | Image-space.
            %
            % Formats:    Raw | ExportedRaw | Cpx | Bruker
            %
            % Description/Algorithm: Generally the receiver samples data denser in readout direction
            %             than specified in the user interface. The denser sampling produces a
            %             larger field-of view (FOV) and prevents fold-over artifacts in readout
            %             direction. Oversampling in readout direction does not increase the scan
            %             time as the receiver can sample at a much higher frequency than desired.
            %             In 3D acquisitions the data is usually oversampled in slice encoding
            %             direction as well to correct for improper slice excitation of the outer
            %             slices. Oversamlling in slice encoding direction does affect the scan time
            %             however.
            %             Since the denser sampling basically results in the acquisiton of a larger
            %             FOV we have to remove the oversampling by cropping the image in
            %             image-space and thereby preserving the signal-to-noise ratio (SNR).
            %             Removing the oversampling in k-space would not be beneficial as we
            %             would lose SNR.
            %             The Linear prediction (LP) based algorithm [Fuchs A et al., ISMRM 2012]
            %             can be used for spectroscopy data to downsample the readout direction. It
            %             operates in the time domain without losing SNR and leads to much less pronounced
            %             truncation artifacts.
            %             In MRecon the RemoveOversampling can be called at three different
            %             locations:
            %
            %             1) Before SortData: The oversampling is only removed in readout direction
            %                                 by fourier transforming the data in that direction and
            %                                 cropping it. This is usually done to reduce the amount
            %                                 of data and save memory.
            %             2) After SortData:  The oversampling is removed in all directions
            %                                 by fourier transforming the data in every oversampled
            %                                 direction and cropping it in image space.
            %             3) In image-space:  The oversampling is removed in all directions by
            %                                 cropping the data.
            %
            %             For spectroscopy data using the LP based algorithm should be in k-space and
            % 			  preferable before SortData. Only the readout direction will be downsampled with
            % 			  that method.
            %
            % Notes:    - The oversampling is not removed for data which has to be gridded when
            %             called before SortData (radial, spiral, EPI).
            %           - The implementation of the LP based algorithm is at the moment much slower and it
            %             requires Matlab's Signal Processing Toolbox to be installed.
            
            if ~MR.Parameter.ReconFlags.isread
                error( 'Error in RemoveOversampling: Please read the data first' );
            end
            
            dims2remove = [ strcmpi( MR.Parameter.Recon.RemoveMOversampling, 'yes' ), strcmpi( MR.Parameter.Recon.RemovePOversampling, 'yes' ), strcmpi( MR.Parameter.Recon.RemovePOversampling, 'yes' )];                                                                                            
            p = inputParser;
            
            
            dimsValidator = @(x) length(x) == 3 && (islogical(x) || (isnumeric(x) && all(x <= 1) && all(x>=0)));
            forceValidator = @(x) islogical(x) || (isnumeric(x) && all(x <= 1) && all(x>=0));
            addParameter(p,'force', false, forceValidator);
            addParameter(p,'dims', dims2remove, dimsValidator);
            parse(p,varargin{:});
            
            
            % If data is in k-space transform it the image space to remove
            % oversampling
            if ~isempty(MR.Data)
                
                if MR.Parameter.ReconFlags.issorted || p.Results.force || ...
                        ( any(strcmpi( MR.Parameter.Scan.AcqMode, {'Cartesian'} )) && ( ~strcmpi( MR.Parameter.Gridder.Preset, {'Epi'} ) || MR.Parameter.ReconFlags.isgridded ) )
                    
                    MR.DataClass.Convert2Cell;
                    if isempty( MR.Parameter.LabelLookupTable )
                        MR.Parameter.LabelLookupTable = cellfun( @(x)[], MR.Data, 'UniformOutput', 0);
                    else
                        MR.Parameter.LabelLookupTable = Helper.Convert2Cell( MR.Parameter.LabelLookupTable );
                    end
                    
                    [MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable] = Helper.check_cell_sizes( MR.Parameter.Encoding.WorkEncoding, MR.Parameter.LabelLookupTable, MR.Data );
                    
                    istransformed = 0;
                    
                    x_ovs = MR.Parameter.Encoding.WorkEncoding.KxOversampling;
                    y_ovs = MR.Parameter.Encoding.WorkEncoding.KyOversampling;
                    z_ovs = MR.Parameter.Encoding.WorkEncoding.KzOversampling;
                    fft_x = 1;
                    fft_y = 1;
                    fft_z = 1;
                    
                    if any( ~MR.Parameter.ReconFlags.isimspace )
                        istransformed = 1;
                        
                        fft_dim_bak = MR.Parameter.Encoding.WorkEncoding.FFTDims;
                        %                 MR.Parameter.Encoding.WorkEncoding.FFTDims = find( MR.Parameter.Encoding.WorkEncoding.Oversampling > 1 );
                        fft_x = any( cellfun( @(x) x > 1 && ~MR.Parameter.ReconFlags.isimspace(1), MR.Parameter.Encoding.WorkEncoding.KxOversampling(1,:) ) );
                        
                        if MR.Parameter.ReconFlags.issorted
                            if ~isempty( MR.Parameter.Encoding.KyRange )
                                fft_y = any( cellfun( @(x) max([x, 1]) > 1 && ~MR.Parameter.ReconFlags.isimspace(2), MR.Parameter.Encoding.WorkEncoding.KyOversampling(1,:) ) );
                            else
                                fft_y = 0;
                            end
                            if ~isempty( MR.Parameter.Encoding.KzRange )
                                fft_z = any( cellfun( @(x) max([x, 1]) > 1 && ~MR.Parameter.ReconFlags.isimspace(3), MR.Parameter.Encoding.WorkEncoding.KzOversampling(1,:) ) );
                            else
                                fft_z = 0;
                            end
                            % if oversampling removal was disabled then
                            % return immediately (saves a fourier trafo)
                            if strcmpi( MR.Parameter.Recon.RemoveMOversampling, 'no' ) && strcmpi( MR.Parameter.Recon.RemovePOversampling, 'no' )
                                return;
                            end
                        else
                            fft_y = 0;
                            fft_z = 0;
                            x_ovs = MR.Parameter.Encoding.WorkEncoding.KxOversampling;
                            y_ovs = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.KyOversampling, 'UniformOutput',0);
                            z_ovs = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.KzOversampling, 'UniformOutput',0);
                            
                            % if x oversampling removal was disabled then
                            % return immediately (saves a fourier trafo)
                            if strcmpi( MR.Parameter.Recon.RemoveMOversampling, 'no' )
                                return;
                            end
                        end
                        
                        if all( [fft_x, fft_y, fft_z] == 0 )
                            istransformed = 0;
                        else
                            MR.Parameter.Encoding.WorkEncoding.FFTDims = [fft_x, fft_y, fft_z];
                        end
                        
                        % spectro begin ----------------------------
                        % Spectro downsampling method operates in the time domain, so try to prevent Fourier
                        % Transformation of the data here
                        if ( MR.isSpectro && strcmpi(MR.Parameter.Spectro.Downsample, 'yes') )
                            
                            if MR.Parameter.ReconFlags.isimspace(1)
                                MR.Parameter.Encoding.FFTDims = [1 0 0];
                                MR.I2K;
                                istransformed = 1;
                            else
                                istransformed = 0;
                                MR.Parameter.Encoding.WorkEncoding.FFTDims = fft_dim_bak;
                            end
                        else
                            % spectro end ------------------------------
                            if  istransformed
                                MR.K2I;
                            end
                        end
                    end
                    
                    MR.DataClass.Convert2Cell;
                    
                    MR.Parameter.UpdateImageInfo = 0;
                    
                    if ~p.Results.dims(1)
                        x_ovs = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.KyOversampling, 'UniformOutput',0);
                    end
                    if ~p.Results.dims(2)
                        y_ovs = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.KyOversampling, 'UniformOutput',0);                       
                    end
                    if ~p.Results.dims(3)                       
                        z_ovs = cellfun( @(x)[], MR.Parameter.Encoding.WorkEncoding.KzOversampling, 'UniformOutput',0);
                    end
                    
                    % adapt the ky_range if the data wasn't sorted yet
                    if ~MR.Parameter.ReconFlags.issorted
                        MR.Parameter.Encoding.WorkEncoding.KxRange = cellfun( @(x,y) Helper.set_k_ranges( x, y ), MR.Parameter.Encoding.WorkEncoding.KxRange, ...
                            x_ovs, 'UniformOutput',0 );
                    end
                    
                    x_ovs_temp = cell(size(x_ovs));
                    y_ovs_temp = cell(size(y_ovs));
                    z_ovs_temp = cell(size(z_ovs));
                    for i = 1:size( x_ovs, 2 )
                        ovs = [ x_ovs(:,i), y_ovs(:,i), z_ovs(:,i) ];
                        mps = MR.Parameter.Encoding.WorkEncoding.MPS;
                        ovs= [ ovs(:, mps(1)), ovs(:, mps(2)), ovs(:, mps(3)) ];
                        x_ovs_temp(:,i) = ovs(:,1);
                        y_ovs_temp(:,i) = ovs(:,2);
                        z_ovs_temp(:,i) = ovs(:,3);
                    end
                    x_ovs = x_ovs_temp;
                    y_ovs = y_ovs_temp;
                    z_ovs = z_ovs_temp;
                    
                    % spectro begin ----------------------------
                    if ( MR.isSpectro && strcmpi(MR.Parameter.Spectro.Downsample, 'yes') )
                        
                        MR.MRImpl.SpectroDownsample( x_ovs );
                        
                    else
                        % spectro end ------------------------------
                        try
                            % Shift the image by 1 pixel depending on the
                            % fat shift. The fat shift flips the axis, so
                            % when the data size is even we have to cut an area which is shifted by one pixel to match the data with different fat shifts  area
                            MR.Data = cellfun( @(x,y,z) Helper.apply_fat_shift_pixel_shift( x, y, z, MR.Parameter.Labels ), MR.Data, ...
                                MR.Parameter.LabelLookupTable, x_ovs, 'UniformOutput',0 );
                        catch
                        end
                        
                        [MR.Data, MR.Parameter.LabelLookupTable] = cellfun( @(x, y, ox, oy, oz)ImageProduction.rem_ovs( x, y, ox, oy, oz ), ...
                            MR.Data, ...
                            MR.Parameter.LabelLookupTable, ...
                            x_ovs, ...
                            y_ovs, ...
                            z_ovs, 'UniformOutput', 0 );
                    end
                    
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    MR.Parameter.LabelLookupTable = Helper.UnconvertCell( MR.Parameter.LabelLookupTable );
                    
                    % Transform it back
                    if istransformed
                        % spectro begin ----------------------------
                        if ( MR.isSpectro && strcmpi(MR.Parameter.Spectro.Downsample, 'yes') )
                            MR.K2I;
                            MR.Parameter.Encoding.FFTDims = fft_dim_bak;
                        else
                            % spectro end ------------------------------
                            MR.I2K;
                            MR.Parameter.Encoding.WorkEncoding.FFTDims = fft_dim_bak;
                        end
                    end
                    
                    if fft_x && p.Results.dims(1)
                        MR.Parameter.ReconFlags.isoversampled(1) = 0;
                    end
                    if fft_y && p.Results.dims(2)
                        MR.Parameter.ReconFlags.isoversampled(2) = 0;
                    end
                    if fft_z && p.Results.dims(3)
                        MR.Parameter.ReconFlags.isoversampled(3) = 0;
                    end
                end
            end
        end
        function ScaleData( MR )
            % ScaleData: Clips the data such that 98% of its histogram is preserved
            %
            % Syntax:     r.ScaleData;
            %
            % Parameters used: None
            %
            % Location:   Image-space
            %
            % Formats:    Raw | ExportedRaw | Cpx | Rec | Bruker
            %
            % Description/Algorithm: A histogram is calculated over all values in the Data array and
            %             the value is found which includes 98% of the histograms integral. The data
            %             is then clipped at that value. This function is usefull to remove unwanted
            %             signal peaks.
            
            if any( ~MR.Parameter.ReconFlags.isimspace )
                error( 'The data must be scaled in image-space' );
            end
            
            if ~isempty(MR.Data)
                MR.Parameter.UpdateImageInfo = 0;
                MR.DataClass.Convert2Cell;
                
                MR.Data = cellfun( @(x)ImageProduction.scale_image( x ), MR.Data, 'UniformOutput', 0 );
                
                MR.Parameter.UpdateImageInfo = 1;
                MR.Data = Helper.UnconvertCell( MR.Data );
            end
            
        end
        function RotateImage( MR )
            % RotateImage: Rotates the images such that they are oriented in the same way as on the
            % scanner console.
            %
            % Syntax:     r.RotateImage;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ReconParsDoc.RotateImage')">Parameter.Recon.RotateImage</a>: {'yes'} | 'no'
            %             Enables/disables the image rotation.
            %           - <a href="matlab:helpwin('ScanParsDoc.ijk')">Parameter.Scan.ijk</a>:
            %             Specifies the coordinate system of the images before the rotation (current
            %             coordinate system).
            %           - <a href="matlab:helpwin('ScanParsDoc.REC')">Parameter.Scan.REC</a>:
            %             Specifies the coordinate system of the images after the rotation
            %             (coordinate system of the Rec images).
            %
            % Recon Flag: <a href="matlab:helpwin('ReconFlagParsDoc.isrotated')">Parameter.ReconFlags.isrotated</a>
            %
            % Location:   Image-space
            %
            % Formats:    Raw
            %
            % Description/Algorithm: After reading the raw data it is oriented in The MPS system,
            %             meaning that the rows in the Data matrix are aligned along the measurement
            %             direction and the columns along the phase encoding direction. The MPS
            %             system however heavily depends on the chosen scan parameters, such as
            %             orientation, fold-over direction, fat shift direction etc. The
            %             reconstructed images, displayed on the scanner console on the other hand
            %             are oriented in a well defined fashion expressed in the patient coordinate
            %             system (right-left (RL), anterior-posterior (AP), feet-head (FH) ). The
            %             RotateImage function transforms the Data from the current coordinate
            %             system (usually MPS) the the final orientation displayed on the scanner
            %             console.
            %
            % Examples:   The Rotate Image basically performes a transformation between the system
            %             specified in r.Parameter.Scan.ijk and r.Parameter.Scan.REC. These
            %             coordinates systems are always expressed in the patient system (RL, AP,
            %             FH). Lets assume these values are given as:
            %
            %                   r.Parameter.Scan.ijk = [RL, FH, AP]
            %                   r.Parameter.Scan.REC = [HF, RL, AP]
            %
            %             To transform the ijk system to the REC system, we have to switch the first
            %             and the secons axis of the Data matrix (permute it):
            %
            %                   [RL, FH, AP]  -->  [FH, RL, AP]
            %
            %             and then flip the first axis (flipdim):
            %
            %                   [FH, RL, AP]  -->  [HF, RL, AP]
            
            if ~isempty(MR.Data) && strcmpi( MR.Parameter.Recon.RotateImage, 'yes' )
                if any( ~MR.Parameter.ReconFlags.isimspace )
                    error( 'Error in RotateImage: The data must be in image-space' );
                end
                
                if ~isempty( MR.Parameter.Scan.ijk ) && ~isempty( MR.Parameter.Scan.REC )
                    
                    if isfield( MR.Parameter.Labels, 'StackIndex' )
                        try
                            stack_nr = MR.Parameter.Labels.StackIndex( MR.Parameter.Parameter2Read.loca + 1 );
                        catch
                            stack_nr = [];
                        end
                    else
                        stack_nr = [];
                    end
                    
                    MR.Parameter.UpdateImageInfo = 0;
                    MR.DataClass.Convert2Cell;
                    
                    MR.Data = cellfun( @(x)ImageProduction.rotate_image_new( x, MR.Parameter.Scan.ijk, MR.Parameter.Scan.REC, stack_nr),  ...
                        MR.Data, 'UniformOutput', 0 );
                    %                     MR.Parameter.Encoding.WorkEncoding.MPS = MR.Parameter.Encoding.WorkEncoding.MPS{1};
                    
                    %
                    for i = 1:size( MR.Parameter.Scan.curFOV, 1 )
                        if i <= size( MR.Parameter.Scan.ijk, 1 )
                            cur_ind = i;
                        else
                            cur_ind = size( MR.Parameter.Scan.ijk, 1 );
                        end
                        P = MRparameter.get_coord_transformation(MR.Parameter.Scan.ijk(cur_ind,:) , MR.Parameter.Scan.REC(cur_ind,:) );
                        MR.Parameter.UpdateCurFOV( [], P, i );
                    end
                    
                    % For flow scans:
                    % Before rotate the flow encoding direction is the same
                    % as in r.Parameter.Scan.MPS and might therefore be
                    % inverted compared to the Philips REC images. After
                    % this function the order of the segments is still
                    % along M-P-S but the direction is always in FH, AP and
                    % RL
                    if any( MR.Parameter.Scan.Venc ~= 0) && MR.Parameter.ReconFlags.issegmentsdivided %flow scan
                        MR.Data = cellfun( @(x)Flow.invert_flow_segments( x, MR.Parameter.Scan.Venc, MR.Parameter.Scan.MPS, stack_nr, strcmpi(MR.Parameter.Scan.PCAcqType, 'Hadamard') ),  ...
                            MR.Data, 'UniformOutput', 0 );
                    end
                    
                    
                    MR.Parameter.Scan.ijk = MR.Parameter.Scan.REC;
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                end
                MR.Parameter.ReconFlags.isrotated = 1;
            end
        end
        function DivideFlowSegments( MR )
            % DivideFlowSegments: Calculates the phase image proportional to the flow, by dividing
            % the two flow segments.
            %
            % Syntax:     r.DivideFlowSegments;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ScanPars.Venc')">Parameter.Scan.Venc</a>:
            %             Used to check if the current scan is a flow acquisition. Has to be unequal
            %             to zero.
            %           - <a href="matlab:helpwin('ReconParsDoc.CoilCombination')">Parameter.Recon.CoilCombination</a>:
            %             Used to check if the phase was preserved during the coil combination. Has
            %             to be set to 'pc'.
            %
            % Location:   image-space
            %
            % Formats:    Raw, ExportedRaw
            %
            % Description/Algorithm: In flow acquisitions two flow segments are acquired, eigher
            %             with negativ and positiv flow encoding (symmetric) or with flow and no
            %             flow encoding (asymmetric). To obtain the quantitative flow values a phase
            %             image proportional to the flow is calculated, by dividing the two
            %             segments. The phase values in the divided image reach from -pi-pi
            %             corresponding to flow values from 0-venc.
            %
            % Notes:    - The two segments are stored in the 10th dimension of the Data array.
            
            if ~isempty(MR.Data) && strcmpi(MR.Parameter.Recon.DivideFlowSegments, 'Yes')
                if any( ~MR.Parameter.ReconFlags.isimspace )
                    error( 'Error in DivideFlowSegments: The DivideFlowSegments has to be performed in image space');
                end
                
                MR.Parameter.UpdateImageInfo = 0;
                
                if ~isempty( MR.Parameter.Scan.Venc ) &&  ...
                        any( any( MR.Parameter.Scan.Venc ~= 0 )) && ...
                        strcmpi( MR.Parameter.Recon.CoilCombination, 'pc' )
                    
                    MR.DataClass.Convert2Cell;
                    MR.Data(1,:) = cellfun( @(x)Flow.divide_segments( x, MR.Parameter.Scan.PCAcqType, MR.Parameter.Recon.TKE ), MR.Data(1,:) , 'UniformOutput', 0);
                    MR.Parameter.UpdateImageInfo = 1;
                    MR.Data = Helper.UnconvertCell( MR.Data );
                    MR.Parameter.ReconFlags.issegmentsdivided = 1;
                end
            end
        end
        function ReconTKE( MR )
            % ReconTKE: Perform Bayes Recon for TKE
            %
            % Syntax:     r.ReconTKE;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ScanPars.Venc')">Parameter.Scan.Venc</a>:
            %             Used to check if the current scan is a flow acquisition. Has to be unequal
            %             to zero.
            %           - <a href="matlab:helpwin('ReconPars.CoilCombination')">Parameter.Recon.CoilCombination</a>:
            %             Used to check if the phase was preserved during the coil combination. Has
            %             to be set to 'pc'.
            %
            % Location:   image-space
            %
            % Formats:    Raw, ExportedRaw
            %
            % Description/Algorithm: Perform Bayes Recon for TKE.
            %
            % Notes:
            if strcmpi(MR.Parameter.Recon.TKE,'Yes')
                if any( ~MR.Parameter.ReconFlags.isimspace )
                    error( 'Error in ReconTKE: The ReconTKE has to be performed in image space');
                end
                if ( MR.Parameter.Recon.FluidDensity <= 0 )
                    error( 'Error in ReconTKE: Please set a fluid density [kg/m^3] in Parameter.Recon.FluidDensity');
                end
                MR.DataClass.Convert2Cell;
                [MR.Data(:,1), MR.Data(:,2)] = cellfun( @(x)Flow.recon_tke( x, MR.Parameter.Recon.FluidDensity, MR.Parameter.Scan.kv, MR.Parameter.Scan.Venc ), MR.Data(1,:) , 'UniformOutput', 0);
                MR.Data = Helper.UnconvertCell( MR.Data ); 
                
                % After recon_tke we have 3 orthogonal flow segments
                max_venc = max(sqrt(sum(MR.Parameter.Scan.Venc.^2, 2)));                                          
                MR.Parameter.Scan.Venc = eye(3).*max_venc;
                MR.Parameter.Scan.kv = eye(3).*pi/max_venc*100;
                
                % the TKE maps were copied to a new cell in MR.Data. Copy
                % Encoding Parameter from cell1 to a new cell2:
                fields = fieldnames(MR.Parameter.Encoding.WorkEncoding);
                for i=1:length(fields)
                    name = char(fields(i));
                    MR.Parameter.Encoding.WorkEncoding.(name)(2) = MR.Parameter.Encoding.WorkEncoding.(name)(1);
                end
                
                if size(MR.Parameter.Scan.Venc, 1) < 4
                    MR.Parameter.Scan.Venc = max(MR.Parameter.Scan.Venc,[], 1);
                end                
            end
        end
        function Average( MR )
            % Average: 	  Averages the data in the 12th dimension of the data array. For spectroscopy the
            %  			  averaging scheme can be customized and can be extended to include also the 5th
            %  			  dimension of the data array.
            %
            % Syntax:     r.Average;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ScanParsDoc.Diffusion')">Parameter.Scan.Diffusion</a>:
            %             Used to check if it is a diffusion scan. Only the magnitude is averaged
            %             for diffusion scans since phase cancelations can occur otherwise.
            %           - <a href="matlab:helpwin('AveragingParsDoc')">Parameter.Spectro.Averaging</a>:
            %             Used to configure customized averaging schemes for spectroscopy.
            %
            % Location:   k-space | image-space
            %
            % Formats:    Raw | ExportedRaw
            %
            % Description/Algorithm: This function averages the data in the 12th dimension of the
            %             data array (where the different averages are stored). This function only has
            %             an effect if immediate averaging is disabled (r.Parameter.Recon.ImmediateAveraging).
            %             If immediate averaging is enabled, then the data is averaged in the
            %             SortData function already. Some data, however should not be averaged in
            %             k-space (e.g. diffusion scans). That is why this function is called in
            %             image space to average the remaining averages which have not been
            %             processed yet.
            %             For spectroscopy more complex averaging schemes can be configured for different
            %             situations. Block wise averaging can be set up and could be useful to increase
            %             SNR in intermediate post processing steps. Furthermore add/subtract schemes for
            %             FIDs and dynamics can be customized to be used for example in spectral editing
            %             experiments.
            %
            % Notes:    - For diffusion scans only the magnitude is averaged.
            
            if strcmpi( MR.Parameter.Recon.Average, 'yes' )
                MR.DataClass.Convert2Cell;
                if ~MR.Parameter.ReconFlags.issorted
                    error( 'Error in Average: Please sort the data first' );
                end
                % spectro begin ----------------------------
                if MR.isSpectro
                    MR.MRImpl.SpectroAverage;
                end
                % spectro end ------------------------------
                average = any( cellfun(@(x)size(x, 12) > 1, MR.Data ) );
                if average && ~MR.isSpectro
                    if strcmpi( MR.Parameter.Scan.Diffusion, 'yes')
                        % diffusion scans should be averages after the epi correction
                        if( ~MR.Parameter.ReconFlags.isdepicorr )
                            warning( 'Warning in Average: Diffusion images should be averaged after the the EPI correction');
                        end
                        
                        MR.Data = cellfun( @(x)mean( abs(x), 12 ), MR.Data, 'UniformOutput', 0 );
                    else
                        MR.Data = cellfun( @(x)mean( x, 12 ), MR.Data, 'UniformOutput', 0 );
                    end
                end
                MR.Data = Helper.UnconvertCell( MR.Data );
                % spectro begin ----------------------------
                MR.Parameter.ReconFlags.isaveraged = 1;
                % spectro end ------------------------------
            end
        end
        
        % ---------------------------------------------------------------%
        % Operations on REC Data
        % ---------------------------------------------------------------%
        function RescaleREC( MR )
            % RescaleREC: Rescales REC data to the original floating point value
            %
            % Syntax:     r.RescaleREC;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('InfoPars.RescaleSlope')">Parameter.ImageInformation.RescaleSlope</a>:
            %             The rescale factor.
            %           - <a href="matlab:helpwin('InfoPars.RescaleIntercept')">Parameter.ImageInformation.RescaleIntercept</a>:
            %             The Rescale offset.
            %
            % Location:   image-space
            %
            % Formats:    Rec
            %
            % Description/Algorithm: The data type in Rec data is 16bit unsigned integer. However it
            %             is possible to restore the original floating point values with the rescale
            %             factor/offset which is stored in the parfile. The floating point value is
            %             given by:
            %
            %             floating_point_value = value_in_recfile * rescale_slope + rescale_intercept
            
            try
                data = MR.Data;
                for i = 1:length( MR.Parameter.ImageInformation(:) )
                    if isempty(MR.Parameter.ImageInformation(i).RescaleSlope) || isempty(MR.Parameter.ImageInformation(i).RescaleIntercept)
                        continue;
                    end
                    if( MR.Parameter.ImageInformation(i).ImageType ~= 3 )
                        data(:,:,i) = data(:,:,i).*MR.Parameter.ImageInformation(i).RescaleSlope + MR.Parameter.ImageInformation(i).RescaleIntercept;
                    else
                        % rescale the phase images
                        data(:,:,i) = data(:,:,i)./4095.*2*pi - pi;
                    end
                end
                MR.Data = data;
            catch exeption                
                warning( 'MATLAB:MRecon', 'could not Rescale the REC images. Error in: \nfunction: %s\nline: %d\nerror: %s', exeption.stack(1).name, exeption.stack(1).line, exeption.message );
            end
        end
        function CreateComplexREC( MR )
            % CreateComplexREC: Creates complex images from Philips REC data.
            %
            % Syntax:     r.CreateComplexREC;
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('MRparameterDoc.DataFormat')">MR.Parameter.DataFormat</a>:
            %             Is used to check if REC data is passed to the function
            %
            % Location:   Image-space.
            %
            % Formats:    Rec
            %
            % Description/Algorithm: Magnitude and phase are stored as seperate images in Philips
            %             REC files. This means that they are also stored as seprate images in the
            %             MRecon Data array when the file is read. This function creates complex
            %             data from the seperate magnitude and phase images and stores it in the Data
            %             array. The complex image is calculated as:
            %
            %                        complex_image = magnitude_image * exp( i*phase_image)
            %
            % Notes:    - Not all REC files contain magnitude or phase images. The stored image
            %             types are selected on the scanner user interface. If the REC data does not
            %             contain both types this function has no effect.
            
            if ~strcmpi( MR.Parameter.DataFormat, 'Rec' )
                error( 'Error in CreateComplexREC: This function can only be executed on REC data');
            end
            
            try
                if size(MR.Data, 7) > 1
                    MR.Data = MR.Data(:,:,:,:,:,:,1,:,:,:,:,:).*exp( 1i.* MR.Data(:,:,:,:,:,:,2,:,:,:,:,:));
                end
            catch exeption
                warning( 'MATLAB:MRecon', 'could not Rescale the REC images. Error in: \nfunction: %s\nline: %d\nerror: %s', exeption.stack(1).name, exeption.stack(1).line, exeption.message );
            end
        end
        
        % ---------------------------------------------------------------%
        % Image Viewing
        % ---------------------------------------------------------------%
        function ShowData( MR )
            % ShowData: Displays the Data array using MRecon's built-in 3D image viewer.
            %
            % Syntax:     r.ShowData;
            %
            % Parameters used: None
            %
            % Location:   k-space | Image-space
            %
            % Formats:    Raw | ExportedRaw | Cpx | Rec | Bruker
            %
            % Description/Algorithm: Opens the built-in 3D image viewer of MRecon and displays the
            %             Data array.
            %
            % Notes:    - The MRecon image viewer is a standalone tool and can be used to visualize
            %             any Matlab matrix. It can be called as:
            %
            %                   image_slide( Matrix )
            %
            %             r.ShowData therefore is aequivalent to:
            %
            %                   image_slide( r.Data )
            
            if isempty( MR.Data )
                error( 'Error in ShowData: The data array is empty')
            end
            if iscell( MR.Data )
                for i = 1:size(MR.Data, 1)
                    for j = 1:size(MR.Data, 2)
                        if ~isempty( MR.Data{i,j} )
                            if isreal( MR.Data )
                                imslide(MR.Data{i,j})
                            else
                                imslide(angle(MR.Data{i,j}))
                                imslide(abs(MR.Data{i,j}))
                            end
                        end
                    end
                end
            else
                if ~isempty( MR.Data)
                    if isreal( MR.Data )
                        imslide(MR.Data)
                    else
                        imslide(angle(MR.Data))
                        imslide(abs(MR.Data))
                    end
                end
            end
        end
        function CreateVideo( MR, Dimension, Filename, Framerate, OutputFormat, StartFrame, EndFrame )
            % CreateVideo: Creates a video from reconstructed images stored in the MRecon object
            %
            % Syntax:     CreateVideo(Dimension, Filename, Framerate, OutputFormat, StartFrame, EndFrame );
            %
            % Inputs:
            %   Dimension       Dimension of the Data array over which the the video is generated.
            %                   Default = 5 (dynamics)
            %   Filename        Name of the output video file. Default = generated from parameters
            %   Framerate       Desired framerate of the output video (frames/s). Default = 15
            %   OutputFormat    the output format('avi' or 'gif'). Default = 'avi'
            %   StartFrame      First frame of the video. Default = 1
            %   EndFrame        Last frame of the video. Default = nr_images
            %
            %
            % Important Notes:  Animated .gifs with framerates higher than approx. 15 fps are not
            %                   played at the right speed with some video player (e.g. Microsoft Power Point)
            %                   Use .avi format for higher framerates.
            
            if( nargin > 7 )
                error( 'Too many input arguments');
            elseif nargin > 6
            elseif nargin > 5
                EndFrame = -1;
            elseif nargin > 4
                StartFrame = 1;
                EndFrame = -1;
            elseif nargin > 3
                OutputFormat = 'avi';
                StartFrame = 1;
                EndFrame = -1;
            elseif nargin > 2
                Framerate = 15;
                OutputFormat = 'avi';
                StartFrame = 1;
                EndFrame = -1;
            elseif nargin > 1
                Filename = '';
                Framerate = 15;
                OutputFormat = 'avi';
                StartFrame = 1;
                EndFrame = -1;
            else
                Dimension = 5;
                Filename = '';
                Framerate = 15;
                OutputFormat = 'avi';
                StartFrame = 1;
                EndFrame = -1;
            end
            
            
            % error checks
            if( Dimension < 1 || Dimension > 13 )
                error( 'The dimension has to be between 1 and 13');
            end
            if( Framerate < 0 )
                error( 'The framerate has to be larger than 0');
            end
            
            selection_string = ['::11111111111';',,,,,,,,,,,,,'];
            selection_string(1,Dimension) = ':';
            selection_string = selection_string(:)';
            selection_string = selection_string(1:end-1);
            eval( ['data = squeeze(abs(MR.Data(', selection_string, ')));']);
            
            % Max is the maximum pixel intensity to which all output frames are scaled (the first few dynamics are omitted as spin saturation effects may occur)
            Max=max(max(max(data(:,:,ceil(size(data,3)/20:size(data,3))))));
            % Output filename is generated automatically from Scan
            if(isempty(Filename))
                Filename=strcat('scan',num2str(MR.Parameter.Scan.AcqNo),'_framerate',num2str(ceil(Framerate)),'fps.',OutputFormat);
            end
            if( EndFrame == -1)
                EndFrame = size(data,3);
            end
            if( EndFrame > size(data, 3) )
                EndFrame = size(data,3);
            end
            
            if( StartFrame <= 0 )
                error( 'The start frame has to be larger than 0');
            end
            if( StartFrame > EndFrame )
                error( 'The start frame has to be smaller than the end frame');
            end
            
            % Export Videos
            switch OutputFormat
                
                case 'avi'
                    vidObj = VideoWriter(Filename); % creates Matlab VideoWriter object
                    vidObj.FrameRate = Framerate;
                    open(vidObj);
                    figure;
                    for frame=StartFrame:EndFrame
                        image(data(:,:,frame),'CDataMapping','scaled'), caxis([0,Max]);set(gca,'DataAspectRatio',[1 1 1]);
                        axis off, colormap gray;
                        f=getframe(gcf);
                        writeVideo(vidObj, f);
                    end
                    close(gcf)
                    close(vidObj);
                    
                case 'gif' % .gif generation loop
                    delaytime = 1/Framerate;
                    for i=StartFrame:EndFrame
                        I = mat2gray(data(:,:,i)); % creates an array of gray values I between [0,1];
                        [X, map] = gray2ind(I,256); % transforms I into an indexed image X, and creates a gray colormap with 256 entries
                        if i == StartFrame % sets the animated gif to run in an endless loop
                            imwrite(X,map,Filename,'gif','LoopCount',Inf,'DelayTime',delaytime);
                        else
                            imwrite(X,map,Filename,'gif','WriteMode','append','DelayTime',delaytime);
                        end
                    end
                otherwise
                    error('Video format not supported');
            end
        end
        
        % ---------------------------------------------------------------%
        % Coordinate System Transformations
        % ---------------------------------------------------------------%
        function varargout = Transform( MR, varargin )
            % Transform: Transforms coordinates from one scanner coordinate system to another
            %
            % Syntax:     [xyz_transformed, A] = r.Transform( xyz, from_system, to_system, stack_nr);
            %               or
            %             A = r.Transform( from_system, to_system, stack_nr);
            %
            %
            % Parameters used:
            %           - <a href="matlab:helpwin('ScanParsDoc.ijk')">Parameter.Scan.ijk</a>, <a href="matlab:helpwin('ScanParsDoc.MPS')">Parameter.Scan.MPS</a>, <a href="matlab:helpwin('ScanParsDoc.xyz')">Parameter.Scan.xyz</a>, <a href="matlab:helpwin('ScanParsDoc.REC')">Parameter.Scan.REC</a>
            %             The different scanner coordinates system expressed in the patient
            %             coordinate system.
            %           - <a href="matlab:helpwin('ScanParsDoc.Angulation')">Parameter.Scan.Angulation</a>, <a href="matlab:helpwin('ScanParsDoc.Offcentre')">Parameter.Scan.Offcentre</a>:
            %             Angulation and offcentre of the scan.
            %           - <a href="matlab:helpwin('ScanParsDoc.curFOV')">Parameter.Scan.curFOV</a>:
            %             The current field of view of the scan. curFOV corresponds to the FOV of
            %             the current data and does change during the reconstruction process.
            %             Removing the oversampling for example reduces the FOV.
            %           - <a href="matlab:helpwin('ScanParsDoc.SliceGap')">Parameter.Scan.SliceGap</a>:
            %             The slice gap of the scan.
            %
            %           - <a href="matlab:helpwin('EncodingParsDoc.XRes')">Parameter.Encoding.XRes</a>, <a href="matlab:helpwin('EncodingParsDoc.YRes')">Parameter.Encoding.YRes</a>, <a href="matlab:helpwin('EncodingParsDoc.ZRes')">Parameter.Encoding.ZRes</a>
            %             <a href="matlab:helpwin('EncodingParsDoc.KyOversampling')">Parameter.Encoding.KyOversampling</a>, <a href="matlab:helpwin('EncodingParsDoc.KzOversampling')">Parameter.Encoding.KyOversampling</a>:
            %             These parameters are only used if the Data matrix is empty. The default
            %             matrix size used in the Transform function, when no data is present, is the
            %             one which we would obtain after the SENSE unfolding (which has to be used
            %             for calculating the coil sensitivities). It is given by:
            %                   [XRes, YRes*KyOversampling, ZRes*KzOversampling]
            %
            % Location:   k-space | image-space
            %
            % Formats:    Raw | Rec
            %
            % Description/Algorithm: Please see the <a href="matlab:open('CoordinateTransformations.pdf')">transformation manual</a> for more infomation on the
            %             different coordinate systems and examples.
            %
            %             To transform a coordinate from one coordinate system to another the
            %             following parameters have to be known:
            %
            %             - Patient Position     (defines coordinate system)
            %             - Patient Orientation  (defines coordinate system)
            %             - Image Orientation    (defines coordinate system)
            %             - Fold Over Direction  (defines coordinate system)
            %             - Fat Shift Direction  (defines coordinate system)
            %
            %             - Angulation           (fixed, used in Transform)
            %             - Offcentre            (fixed, used in Transform)
            %             - Slice Gap            (fixed, used in Transform)
            %             - Resolution           (variable, used in Transform: Resolution = curFOV / size(Data) )
            %             - Image Matrix         (variable, used in Transform for ijk, REC systems)
            %
            %             The first 5 parameters are used to define the different coordinate systems
            %             (ijk, MPS, RAF, xyz, REC), which are stored in Parameter.Scan and used in
            %             the Transform function. Therefore the parameters itself are not needed in
            %             the Transform function anymore.
            %             Basically the parameters which are used in the Transform function can be
            %             divided into 2 groups: fixed and variable parameters. The fixed ones, such
            %             as the image angulation, do not change during the reconstruction process
            %             whereas the variable ones can change. The image resolution for example
            %             changes if we perform a k-space zero filling. In transform the current
            %             image resolution is calculated by dividing the current FOV by the current
            %             matrix size. Therefore if you are performing your own reconstruction and
            %             are not using the built-in MRecon functions (e.g. RemoveOversampling)
            %             always make sure that Parameter.Scan.curFOV corresponds to the current FOV
            %             of your data when calling the Transform function.
            
            if isfield( MR.Parameter.Labels, 'StackIndex' )
                stacks = unique(MR.Parameter.Labels.StackIndex);
            else
                stacks = 0;
            end
            nr_stacks = length(stacks);
            matrix_size = [];
            fov = [];
            angulation = [];
            offcentre = [];
            slice_gap = [];
            
            if isempty( varargin )
                error( 'Error in Transform: please specify the input and output coordinate system' );
            end
            if ischar( varargin{1} )
                matrix_only_mode = 1;
            else
                matrix_only_mode = 0;
            end
            
            if matrix_only_mode
                if length( varargin ) < 2
                    error( 'Error in Transform: please specify the input and output coordinate system' );
                end
                if ~ischar( varargin{2} )
                    error( 'Error in Transform: The output coordinate system must be a string' );
                end
                if length( varargin ) == 2
                    stack = 1:nr_stacks;
                    option_start = 3;
                else
                    if ischar( varargin{3} )
                        option_start = 3;
                        stack = 1:nr_stacks;
                    else
                        stack = varargin{3};
                        option_start = 4;
                    end
                end
                from = varargin{1};
                to = varargin{2};
                x = [1,1,1];
            else
                if length( varargin ) < 3
                    error( 'Error in Transform: please specify the input and output coordinate system' );
                end
                if ~ischar( varargin{3} )
                    error( 'Error in Transform: The output coordinate system must be a string' );
                end
                if length( varargin ) == 3
                    stack = 1;
                    option_start = 4;
                else
                    if ischar( varargin{4} )                        
                        error( 'Error in Transform: The stack number must be numeric' );
                    else
                        stack = varargin{4};
                        option_start = 5;
                    end
                end
                from = varargin{2};
                to = varargin{3};
                x = varargin{1};
            end
            
            if length( varargin ) >= option_start
                for i = option_start:2:length( varargin )
                    if ischar(varargin{i}) && strcmpi( varargin{i}, 'MatrixSize' )
                        if length( varargin ) < i+1
                            error( 'Error in Transform: Please specify an option value for MatrixSize' );
                        else
                            matrix_size = varargin{i+1};
                        end
                    end
                    if ischar(varargin{i}) && strcmpi( varargin{i}, 'FOV' )
                        if length( varargin ) < i+1
                            error( 'Error in Transform: Please specify an option value for FOV' );
                        else
                            fov = varargin{i+1};
                        end
                    end
                    if ischar(varargin{i}) && strcmpi( varargin{i}, 'Angulation' )
                        if length( varargin ) < i+1
                            error( 'Error in Transform: Please specify an option value for Angulation' );
                        else
                            angulation = varargin{i+1};
                        end
                    end
                    if ischar(varargin{i}) && strcmpi( varargin{i}, 'Offcentre' )
                        if length( varargin ) < i+1
                            error( 'Error in Transform: Please specify an option value for Offcentre' );
                        else
                            offcentre = varargin{i+1};
                        end
                    end
                    if ischar(varargin{i}) && strcmpi( varargin{i}, 'SliceGap' )
                        if length( varargin ) < i+1
                            error( 'Error in Transform: Please specify an option value for SliceGap' );
                        else
                            slice_gap = varargin{i+1};
                        end
                    end
                end
            end
            
            [xT, A] = MR.Parameter.Transform( x, from, to, stack, angulation, offcentre, matrix_size, fov, slice_gap );
            if matrix_only_mode
                varargout{1} = A;
            else
                varargout{1} = xT;
                varargout{2} = A;
            end
        end
        
        % ---------------------------------------------------------------%
        % Deep Copy of Class
        % ---------------------------------------------------------------%
        function new = Copy(this)
            % Copy: Deep copy of the MRecon class
            %
            % Syntax:     r_new = r.Copy;
            %
            % Parameters used: None
            %
            % Formats:    All
            %
            % Description/Algorithm: The MRecon class is a handle class which means that assigning it to another
            %             variable is just a pointer operation and does not copy the object.
            %             Therefore the syntax:
            %                                           r_new = r;
            %
            %             does not create a copy of r. Instead r_new does point to the same object
            %             which means that when we change r_new we also change r. To create a deep
            %             copy of our object we have to call the Copy function:
            %
            %                                           r_new = r.Copy;
            
            mc = eval( ['?',class(this)] );
            new = feval(class(this), 'Empty');
            
            % Copy all non-hidden properties.
            p = mc.Properties;
            for i = 1:length(p)
                if ~p{i}.NonCopyable
                    try
                        % Property i has a Copy function
                        new.(p{i}.Name) = this.(p{i}.Name).Copy;
                    catch
                        try
                            % Property i does not have a Copy function
                            new.(p{i}.Name) = this.(p{i}.Name);
                        catch
                        end
                    end
                end
            end
        end
        
        % ---------------------------------------------------------------%
        % Search for parameter name
        % ---------------------------------------------------------------%
        function Search(this, search_text, filter, use_name_text)
            % Search: Search for a parameter or a parameter value within the MRecon class
            %
            % Syntax:     r.Search( search_string );
            %
            % Parameters used: All
            %
            % Formats:    All
            %
            % Description/Algorithm: The search function iterates through all the parameters
            %             and values of the MRecon class and checks if it contains the search string
            %
            % Notes:      For performance reasons the search function does not check values in arrays
            %             with length > 100000
            %
            % Examples:   To search for the field-of-view (FOV) call:
            %
            %                   >> r.Search('FOV')
            %
            %             which results in:
            %
            %                   Parameter.Scan.FOV = 120 63  240
            %                   Parameter.Scan.curFOV = 240  258   63
            %                   Parameter.Labels.FOV = 240 120 63
            %
            %             In the same way you can also search for values. For example
            %
            %                   >> r.Search('120')
            %
            %             results in:
            %
            %                   Parameter.Scan.FOV = 120 63  240
            %                   Parameter.Labels.FOV = 240 120 63
            %
            
            if( nargin == 2)
                filter = '';
                use_name_text = false;
            end
            if( nargin == 3)
                use_name_text = false;
            end
            
            this.Parameter.Search(search_text, filter, use_name_text);
            
            if( isempty(filter) || strcmpi(filter, 'mrecon'))
                disp(' ');
                disp(' ');
                disp( '--------------------------------------' );
                disp( 'Search results in MRecon parameter:' );
                disp( '--------------------------------------' );
                Helper.search( this, search_text, '' );
            end
        end
        
        % ---------------------------------------------------------------%
        % Compare the parameters of 2 objects
        % ---------------------------------------------------------------%
        function equal = Compare(this, other, no_disp)
            % Compare: Compares to MRecon classes and displays the differences
            %
            % Syntax:     r.Compare( r_other );
            %
            % Parameters used: All
            %
            % Formats:    All
            %
            % Description/Algorithm: The compare function iterates through all the parameters values
            %             of the MRecon class and compares them with the values of another reconstruction
            %             object.
            %
            % Examples:   To compare the current reconstruction object (r) with another one (r1) call:
            %
            %                   >> r.Compare(r1);
            
            if nargin < 3
                no_disp = false;
            end
            equal = Helper.compare( this, other, '', no_disp );
        end
        
        function eq = eq(this, other)
            % eq: Equal operator overloading to test for equality between 2
            % MRecon objects (e.g. r1 == r2 )
            %
            % Syntax:     r1 == r2;
            %
            % Parameters used: All
            %
            % Formats:    All
            eq = this.Compare(other);
        end      
        
        % ---------------------------------------------------------------%
        % Various
        % ---------------------------------------------------------------%
        function [MemoryNeeded, MemoryAvailable, MaxDataSize] = GetMemoryInformation( MR )
            P = MR.Parameter.Parameter2Read.Copy;
            
            ProfileMask = (ismember(MR.Parameter.Labels.Index.typ,     P.typ )  & ...
                ismember(MR.Parameter.Labels.Index.mix,   P.mix )  & ...
                ismember(MR.Parameter.Labels.Index.dyn,   P.dyn )  & ...
                ismember(MR.Parameter.Labels.Index.card,  P.card ) & ...
                ismember(MR.Parameter.Labels.Index.loca,  P.loca ) & ...
                ismember(MR.Parameter.Labels.Index.echo,  P.echo)  & ...
                ismember(MR.Parameter.Labels.Index.extr1, P.extr1 )& ...
                ismember(MR.Parameter.Labels.Index.extr2, P.extr2 )& ...
                ismember(MR.Parameter.Labels.Index.ky,    P.ky )   & ...
                ismember(MR.Parameter.Labels.Index.kz,    P.kz ));
            
            % Splitting the calculation of the ProfileMask up in two steps to make it work for
            % exported cpx data as well (exported cpx data does not have an ver and rtop label)
            if isfield( MR.Parameter.Labels.Index, 'aver' )
                ProfileMask = (ProfileMask &...
                    ismember(MR.Parameter.Labels.Index.aver,  P.aver )& ...
                    ismember(MR.Parameter.Labels.Index.rtop, P.rtop) );
            end
            NrImagingProfiles = length( find(ProfileMask) );
            DataSizeAfterRead = NrImagingProfiles*max( MR.Parameter.Encoding.DataSizeByte );
            
            Nr3dChunks = length(P.mix) * length(P.chan) * length(P.dyn) * length(P.card) * length(P.loca) * length(P.echo) * length(P.extr1) * length(P.extr2) * length(P.aver);
            Size3dChunk = max( [1, max(MR.Parameter.Encoding.XRes)] ) * max( [1, max(MR.Parameter.Encoding.YRes)] ) * max( [1, max(MR.Parameter.Encoding.ZRes)] ) * ...
                max( [1, max(MR.Parameter.Encoding.KyOversampling)] ) * max( [1, max(MR.Parameter.Encoding.KzOversampling)] ) ./ ...
                prod(MR.Parameter.Scan.SENSEFactor);
            DataSizeAfterZeroFill = Size3dChunk*Nr3dChunks*2*4;
            
            MaxDataSize = max([DataSizeAfterRead,DataSizeAfterZeroFill]);
            MemoryNeeded = 2.3 * MaxDataSize;
            
            try
                [~, sV] = Helper.memory_mac;
                MemoryAvailable = sV.PhysicalMemory.Available;
            catch
                MemoryAvailable = 1.5*MemoryNeeded;
            end
        end        
               
        % ---------------------------------------------------------------%
        % Spectroscopy
        % ---------------------------------------------------------------%
        function EddyCurrentCorrection( MR )
            % EddyCurrentCorrection: Applies Klose [see Klose U, MRM 1990;14(1)] eddy current correction to spectroscopy data.
            %
            % Syntax: 	  r.EddyCurrentCorrection;
            %
            % Parameters used:
            % 			  None.
            %
            % Location:   k-space | image-space.
            %
            % Formats: 	  Raw (Spectro)
            %
            % Description/Algorithm:
            % 			  Applies Klose [see Klose U, MRM 1990;14(1)] eddy current correction to spectroscopy data.
            % 			  This is only possible if during the spectroscopy exam the spectral correction parameter
            % 			  was turned on in the Examcards. This leads to an additional data set of unsuppressed water
            % 			  that is saved together with the water suppressed raw data. This signal in mix 2 is
            % 			  then used to revert the time dependent phase distortions that arise due to eddy currents.
            % 			  Additionally to the time dependent phase component also any zero or linear order phase is
            % 			  removed and no further phase correction to the data is necessary.
            
            if ~MR.isSpectro
                error('This eddy current correction procedure is just intended for spectroscopy data');
            end
            
            switch lower(MR.Parameter.Recon.EddyCurrentCorrection)
                case 'yes'
                    if ~MR.Parameter.ReconFlags.isread
                        error('Please read the data first' );
                    end
                    if ~MR.Parameter.ReconFlags.issorted
                        error('Please sort the data first');
                    end
                    if MR.Parameter.ReconFlags.isecc
                        error('Data is already eddy current corrected');
                    end
                    if ( MR.Parameter.Encoding.NrMixes < 2 )
                        error('No unsuppressed water scan was found');
                    end
                    
                    fft_dim_bak = MR.Parameter.Encoding.FFTDims;
                    
                    % it should be checked here that the data is present in the timedomain
                    if MR.Parameter.ReconFlags.isimspace(1)
                        MR.Parameter.Encoding.FFTDims = [1 0 0];
                        MR.I2K;
                        istransformed = true;
                    else
                        istransformed = false;
                    end
                    
                    nr_ref_scans = MR.Parameter.Encoding.WorkEncoding.NrFids{2};
                    
                    if ( nr_ref_scans > 1 )
                        warning('Ref scan data is not yet averaged. Using an average/dynamic for correction');
                    end
                    
                    fprintf('Applying eddy current correction ...\n');
                    
                    ref_scan = sum(MR.Data{1,2}, 12);
                    ref_scan_phase = exp(-1i*angle(ref_scan));
                    MR.Data{1,1} = bsxfun( @times, MR.Data{1,1}, ref_scan_phase );
                    MR.Data{1,2} = bsxfun( @times, MR.Data{1,2}, ref_scan_phase );                                      
                    
                    % Transform everything back if necessary
                    if istransformed
                        MR.K2I;
                        MR.Parameter.Encoding.FFTDims = fft_dim_bak;
                    end
                    
                    if ( MR.Parameter.Chunk.CurLoop == MR.Parameter.Chunk.NrLoops )
                        MR.Parameter.ReconFlags.isecc = true;
                    end
                    
                otherwise
                    warning('Eddy correction is deselected, skipping processing step ...');
            end
            
        end
        function WriteSdat( MR, varargin )
            % WriteSdat: Exports the spectroscopy data to a Philips binary SDAT file and a corresponding
            % 			 SPAR text file containing scan parameter informations.
            %
            % Syntax:    r.WriteSdat; or r.WriteSdat( filename );
            %
            % Parameters used:
            % 			 - <a href="matlab:helpwin('ScanParsDoc')">Parameter.Scan</a>:
            % 			   Holds all parameters that will be written to the SPAR file.
            %
            % Location:  k-space.
            %
            % Formats:   Raw (Spectro)
            %
            % Description/Algorithm:
            % 			 Exports the spectroscopy data to a Philips binary SDAT file and a corresponding
            % 			 SPAR text file containing scan parameter informations. If no filename is provided
            % 			 the same name and location of the raw file is used for export. For unsuppressed water
            % 			 scans that may be included in the raw data the extension '_ref' will be added to the
            % 			 filename and a separate pair of SDAT/SPAR files will be saved to the same location.
            % 			 The individual measurements of a spectroscopy experiment need to be averaged to 1 final
            % 			 FID. If dynamics are present or CSI data, the individual FIDs/dynamic or FIDs/CSI voxel
            % 			 are organized in a 2 dimensional data matrix. The binary data is stored as interleaved
            % 			 real and imaginary values in VAX floating point format.
            
            if ~MR.Parameter.ReconFlags.isread
                error('Please read the data first' );
            end
            if ~MR.Parameter.ReconFlags.issorted
                error('Expecting already sorted data');
            end
            % TK20151123 begin
            % original code:
            % if any(MR.Parameter.ReconFlags.isimspace)
            if ~((MR.Parameter.ReconFlags.isimspace(1) == 0) && (MR.Parameter.ReconFlags.isimspace(2) == 1) && (MR.Parameter.ReconFlags.isimspace(3) == 1) )
                % TK20151123 end
                error('Data must be in time domain and real space in order to write sdat');
            end
            if (nargin > 2)
                error('Only one input argument should specify the desired filename');
            end
            
            MR.DataClass.Convert2Cell;
            
            if ( any([size(MR.Data{1},MRecon.dim.coil), size(MR.Data{1}, MRecon.dim.hp), size(MR.Data{1}, MRecon.dim.echo), ...
                    size(MR.Data{1}, MRecon.dim.loc), size(MR.Data{1}, MRecon.dim.mix), size(MR.Data{1}, MRecon.dim.extr1), size(MR.Data{1}, MRecon.dim.extr2)] > 1 ))
                error('Only fully reconstructed SV or CSI data can be saved to SDAT/SPAR files');
            end
            
            [path,filename,~] = fileparts(MR.Parameter.Filename.Data);
            path = [path filesep];
            if( MR.Parameter.IsParameter('RFR_STUDY_DICOM_STUDY_DATE')  )
                scan_date_raw = MR.Parameter.GetValue('RFR_STUDY_DICOM_STUDY_DATE');
                if length(scan_date_raw) < 8
                    scan_date = '01/01/1900';
                else
                    scan_date = [scan_date_raw(7:8) '/' scan_date_raw(5:6) '/' scan_date_raw(1:4)]; 	% format: dd/mm/yyyy
                end
            else
                scan_date_raw = char(regexpi(filename, '_\d{8}_\d{7}', 'match'));
                if length(scan_date_raw) < 9
                    scan_date = '01/01/1900';
                else
                    scan_date = [scan_date_raw(7:8) '/' scan_date_raw(4:5) '/' scan_date_raw(6:9)]; 	% format: dd/mm/yyyy
                end
            end
            
            if (nargin > 1)
                [pp, filename, ~] = fileparts(varargin{1});
                if ~isempty(pp)
                    path = pp;
                else
                    path = [pwd filesep];
                end
            end
            
            par.examination_name = 'unknown';
            par.scan_id = filename(21:end);
            %
            par.scan_date =[scan_date(7:10),'.',scan_date(4:5),'.',scan_date(1:2)];
            
            for mix_cnt = 1:MR.Parameter.Encoding.NrMixes
                par.patient_name 			= 'unknown';
                par.patient_birth_date 		= 'unknown';
                par.patient_position 		= ['"' MR.Parameter.Scan.PatientPosition '"'];
                par.patient_orientation 	= ['"' MR.Parameter.Scan.PatientOrientation '"'];
                par.samples 				= MR.Parameter.Encoding.XRes(mix_cnt);
                par.rows 					= MR.Parameter.Scan.Samples(2)*MR.Parameter.Scan.Samples(3)*MR.Parameter.Encoding.WorkEncoding.NrDyn{1};
                par.synthesizer_frequency 	= MR.Parameter.Labels.ResonanceFreq;
                par.offset_frequency 		= 0;
                par.sample_frequency 		= 32e3 / MR.Parameter.Encoding.KxOversampling(1); 	% at the moment spectro HW can do 32k Hz BW -> might change
                par.echo_nr 				= 1;
                par.mix_number 				= 1;
                par.nucleus 				= '1H';
                par.t0_mu1_direction 		= 0;
                par.echo_time 				= MR.Parameter.Scan.TE;
                par.repetition_time 		= MR.Parameter.Scan.TR(mix_cnt);
                par.averages 				= 1;
                par.volume_selection_enable = '"yes"';
                par.volumes 				= 1;
                par.ap_size 				= MR.Parameter.Scan.FOV(1);
                par.lr_size 				= MR.Parameter.Scan.FOV(2);
                par.cc_size 				= MR.Parameter.Scan.FOV(3);
                par.ap_off_center 			= MR.Parameter.Scan.Offcentre(1);
                par.lr_off_center 			= MR.Parameter.Scan.Offcentre(2);
                par.cc_off_center 			= MR.Parameter.Scan.Offcentre(3);
                par.ap_angulation 			= MR.Parameter.Scan.Angulation(1);
                par.lr_angulation 			= MR.Parameter.Scan.Angulation(2);
                par.cc_angulation 			= MR.Parameter.Scan.Angulation(3);
                par.volume_selection_method = 1;
                par.t1_measurement_enable 	= '"no"'; % no information available
                par.t2_measurement_enable 	= '"no"';
                par.time_series_enable 		= '"no"';
                par.phase_encoding_enable 	= '"yes"';
                par.nr_phase_encoding_profiles = max(MR.Parameter.Scan.Samples(2:3));
                par.ps_ap_off_center 		= 0;
                par.ps_lr_off_center 		= 0;
                par.ps_cc_off_center 		= 0;
                par.ps_ap_angulation 		= 0;
                par.ps_lr_angulation 		= 0;
                par.ps_cc_angulation 		= 0;
                % additional pars that occur in Philips SDAT and may just be necessary
                par.si_ap_off_center 		= par.ap_off_center;
                par.si_lr_off_center 		= par.lr_off_center;
                par.si_cc_off_center 		= par.cc_off_center;
                par.si_ap_off_angulation 	= par.ap_angulation;
                par.si_lr_off_angulation 	= par.lr_angulation;
                par.si_cc_off_angulation 	= par.cc_angulation;
                par.t0_kx_direction 		= 50;
                par.t0_ky_direction 		= 50;
                par.nr_of_phase_encoding_profiles_ky = MR.Parameter.Scan.Samples(3) ;
                par.phase_encoding_direction = '"trans"';  %%ACHTUNG%%%
                par.phase_encoding_fov 		= MR.Parameter.Scan.FOV(1);
                par.slice_thickness 		= par.cc_size;
                par.image_plane_slice_thickness = 0;
                par.slice_distance 			= 0;
                par.nr_of_slices_for_multislice = 1;
                par.Spec_imageinplanetransf = '"plusA-plusB"';
                par.spec_data_type 			= 'cf';
                par.spec_sample_extension 	= '[V]';
                par.spec_num_col 			= par.samples;
                par.spec_col_lower_val 		= 0;
                par.spec_col_upper_val 		= 0;
                par.spec_col_extension 		= '[sec]';
                par.spec_num_row 			= par.rows;
                par.spec_row_lower_val 		= 1;
                par.spec_row_upper_val 		= par.rows;
                par.spec_row_extension 		= '[index]';
                par.num_dimensions 			= 3;
                par.dim1_ext 				= '[sec]';
                par.dim1_pnts 				= par.samples;
                par.dim1_low_val 			=  - par.sample_frequency / 2;
                par.dim1_step 				= 1/par.sample_frequency;
                par.dim1_direction 			= 'mu1';
                par.dim1_t0_point 			= 0;
                par.dim2_ext 				= '[num]';
                par.dim2_pnts 				= par.nr_phase_encoding_profiles;
                par.dim2_low_val 			= 1;
                par.dim2_step 				= 1;
                par.dim2_direction 			= 'x';
                par.dim2_t0_point 			= 50;
                par.dim3_ext 				= '[num]';
                par.dim3_pnts 				= par.nr_of_phase_encoding_profiles_ky;
                par.dim3_low_val 			= 1;
                par.dim3_step 				= 1;
                par.dim3_direction 			= 'y';
                par.dim3_t0_point 			= 50;
                par.echo_acquisition 		= 'FID';
                par.TSI_factor 				= 0;
                par.spectrum_echo_time 		= par.echo_time;
                par.spectrum_inversion_time = 0;
                par.image_chemical_shift 	= 0;
                par.resp_motion_comp_technique = 'NONE';
                par.de_coupling 			= 'NO';
                
                switch mix_cnt
                    case 1
                        fullfilename = sprintf('%s%s', path, filename);
                    case 2
                        fullfilename = sprintf('%s%s_ref',path, filename);
                    otherwise
                        error('Can not handle more than 2 mixes for writing sdat files');
                end
                % according to Niklaus
                fid = permute(squeeze(MR.Data{1,mix_cnt}), [1 3 2]);
                % TK20151123 begin
                % original code:
                % fid = reshape(fid, [par.rows par.samples]);
                fid = reshape(fid, [par.samples par.rows]);
                fid = permute(fid, [2 1]);
                % TK20151123 end
                
                Spectroscopy.write_sdat(fullfilename, fid, par, 'd');
            end
            
            MR.Data = Helper.UnconvertCell( MR.Data );
        end
        function spectro = isSpectro(MR)
            try
                spectro = MR.Parameter.Labels.Spectro;
            catch
                spectro = 0;
            end
        end
              
        % ---------------------------------------------------------------%
        % Set / Get Functions
        % ---------------------------------------------------------------%
        function set.Data(obj,val)
            obj.MRImpl.Data = val;           
        end
        function value = get.Data(obj)
            value = obj.MRImpl.Data;
        end  
        function value = get.Parameter(obj)
            value = obj.MRImpl.Parameter;
        end  
        function props = get.DataProperties(obj)
            props = obj.MRImpl.DataProperties;            
        end
        function value = get.DataClass(obj)
            value = obj.MRImpl.DataClass;
        end                          
    end
    
    methods (Static)
        function Doc()
            % Doc: Opens the MRecon documentation
            %
            % Syntax:     MRecon.Doc;
            web('https://docs.gyrotools.com/mrecon/latest');
        end
        function version = Version()
            % Version: Displays the MRecon version
            MRecon.AddPath;
            version = ['MRecon v', num2str(MRparameter.mVersionMajor), '.', num2str(MRparameter.mVersionMinor), '.', num2str(MRparameter.mVersionBuild)];
            if ~isempty(MRparameter.mVersionSnapshot)
                version = [version, MRparameter.mVersionSnapshot];
            end        
        end
        function Activate(aToken)
            % Activate: Activates MRecon with an activation token
            MRecon.AddPath;
            MReconImpl.Activate(aToken);               
        end
        function ActivateFile(aLicenseFile)
            % ActivateFile: Activates MRecon with an activation file (manual activation)
            MRecon.AddPath;
            MReconImpl.ActivateFile(aLicenseFile);                          
        end
        function CreateActivationFile(output_dir)
            if nargin == 0
                output_dir = './';
            end
            MRparameter.CreateActivationFile(output_dir);            
        end
        function RenewLicense()
            % RenewLicense: Manually renews the license. (Normally there is no need to
            % call this function. The license is automatically renewed when
            % using MRecon)
            MRecon.AddPath;
            MReconImpl.RenewLicense;            
        end
        function CheckForUpdates()
            MRecon.AddPath;
            MReconImpl.CheckForUpdates;
        end
        function CheckoutLicense(duration_days)
            % CheckoutLicense: Checks-out the license for a given number of
            % days. Please note that your license will be locked for that
            % period of time. It cannot be modified or transferred to
            % another computer. 
            MRecon.AddPath;
            MReconImpl.CheckoutLicense(duration_days);
        end
        function LicenseInfo()
            MRecon.AddPath;
            MReconImpl.LicenseInfo;           
        end
        function UserPortal()
            web('https://license.gyrotools.com/portal');
        end
        function Update()
            MRecon.AddPath;
            MReconImpl.Update;           
        end
        function GetApiKey()
            if verLessThan('Matlab', '8.5')
                error('Your Matlab version is not supported anymore. You need at least Matlab 2015a to use MRecon with the license server. Please contact GyroTools');
            end
            MRecon.AddPath;
            api_key = MRparameter.GetApiKey();
            if ~isempty(api_key)
                disp('---------------------------------------------------');
                disp('API Key:');
                disp(api_key);
                disp('---------------------------------------------------');
                warning('This is a personal key associated with your account. For security reasons do not share or expose your API key to anyone.');
            else
                disp('-------------------------------------------------------------------------');
                disp('You license is not eligible for api key access. Please contact GyroTools');
                disp('-------------------------------------------------------------------------');
            end
        end
        function Anonymize(rawfile, pars2anonymize, user_def_values)
            % Anonymize: Anonymizes the rawfile given as input. 
            % A list of parameters which will be anonymized can be gotten with
            % the "GetPars2Anonymize" function. A user defined list can be 
            % passed as input to anonymize a different set of parameters. 
            % Furthermore user defined values can be specified for each 
            % parameter. The values needs to be passed as struct, where the
            % field name corresponds to the parameter name and the field
            % value to the value after anonymization. 
            % WARNING: For implemenation reasons the user defined values will 
            % be cropped/padded to the length of the original value.
            %
            % Examples:   
            %   anonymize a rawfile:
            %
            %       >> MRecon.Anonymize('MyRawfile.raw');
            %
            %   anonymize a rawfile but keep the patients birth date:
            %
            %       >> pars2anonymize = MRecon.GetPars2Anonymize;
            %       >> pars2anonymize{18:19} = []
            %       >> MRecon.Anonymize('MyRawfile.raw', pars2anonymize);
            %
            %   anonymize a rawfile and specify a user defined patient name
            %   and institution name
            %
            %       >> user_def_values.RFR_PATIENT_DICOM_PATIENT_NAME = 'Martin'; 
            %       >> user_def_values.RFR_SERIES_DICOM_INSTITUTION_NAME = 'GyroTools'
            %       >> MRecon.Anonymize('MyRawfile.raw', [], user_def_values);
            MRecon.AddPath;
            P = ParameterReader(rawfile);
            if nargin == 1
                P.Anonymize;
            elseif nargin == 2
                P.Anonymize(pars2anonymize);
            else
                P.Anonymize(pars2anonymize, user_def_values);
            end
        end
        function pars = GetPars2Anonymize()
            % GetPars2Anonymize: Returns a list of parameter names with
            % will be anoynmized when the "Anonymize" function is called
            MRecon.AddPath;
            pars = ParameterReader.get_pars_to_anonymize;
        end
        function AddPath()
            if ~isdeployed
                base_path = which('MRecon.m');
                if isempty( base_path )
                    base_path = which('MRecon.p');
                end
                if ~isempty(base_path)
                    base_path = fileparts(base_path);
                    addpath(genpath(base_path));
                end
            end
        end
        
        function data = GeometryCorrectionFromSinfile(data, sinfile, varargin)
            % GeometryCorrectionFromSinfile: Performs a geometry correction on the data given as input with the parameters from the sin file
            %
            % Syntax:     corrected_data = MRecon.GeometryCorrectionFromSinfile(data, sinfile, options)
            %            
            % Description/Algorithm: A gradient non-linearity correction is performed on
            % the data given as input with the parameters from the sin file.
            %
            % Please note that the images should already be in the correct 
            % radiological display orientation (e.g. for transversal slice
            % orientation A is at the top of the image and L is to the right.)
            %
            % Options: 'patient_position': specifies the patient position.
            % If not given as input then head-first-supine is assumed. The
            % possible patient positions are:
            %       
            %       HFS: Head-first supine
            %       HFP: Head-first prone
            %       HFR: Head-first right
            %       HFL: Head-first left
            %       FFS: Feet-first supine
            %       FFP: Feet-first prone
            %       FFR: Feet-first right
            %       FFL: Feet-first left
            %
            % Examples:            
            %       corrected_data = MRecon.GeometryCorrectionFromSinfile(data, 'my_sinfile.sin', 'patient_position', 'HFR');
            MRecon.AddPath;
            if nargin == 2
                data = MReconImpl.GeometryCorrectionFromSinfile(data, sinfile);    
            elseif nargin == 3
                data = MReconImpl.GeometryCorrectionFromSinfile(data, sinfile, varargin);  
            else
                error('Please specify data and a sinfile');
            end                                                
        end
        function data = GeometryCorrectionUserDef(data, angulation_midslice, offcentre_midslice, patient_position, orientation, resolutionMPS, sCoeffY, cCoeffX, cCoeffZ, r0)
            % GeometryCorrectionUserDef: Performs a geometry correction on the input data with user defined geometry parameters
            %
            % Syntax:     corrected_data = GeometryCorrectionUserDef(data, angulation_midslice, offcentre_midslice, patient_position, orientation, resolutionMPS, sCoeffY, cCoeffX, cCoeffZ, r0)
            %           
            % Inputs: 
            %       data:                Data to be geometry corrected. The dimension order must be the same as in the MRecon Data array (see MRdata)
            %       angulation_midslice: Angulation for the center slice
            %       offcentre_midslice:  Offcentre for the center slice
            %       patient_position:    Patient position [HFS, HFP, HFR, HFL, FFS, FFP, FFR, FFL]
            %       orientation:         Slice orientation [TRA, SAG, COR]
            %       resolutionMPS:       Resolution of the data in mm along the M-P-S axis
            %       sCoeffY:             Geometry correction coefficients for the Y gradient (matrix with size 12x14)
            %       cCoeffX:             Geometry correction coefficients for the X gradient (matrix with size 12x14)
            %       cCoeffZ:             Geometry correction coefficients for the Z gradient (vector with 16 elements)
            %       r0:                  Geometry correction reference radius
            %                       
            % Description/Algorithm: A gradient non-linearity correction is performed on
            % the data given as input with user defined geometry parameters
            %
            % Please note that the data should already be in the correct 
            % radiological display orientation (e.g. for transversal slice
            % orientation A is at the top of the image and L is to the right.)            
            %
            % Examples:   
            %       sinfile = 'MySinfile.sin';
            %       [sCoeffY, cCoeffX, cCoeffZ, r0] = MRecon.GetGeoCorrCoeffsFromSinfile(sinfile);
            %       data = zeros(128,128,128);
            %       data(:, [1:10:128], :) = 1;
            %       corrected_data = MRecon.GeometryCorrectionUserDef(data, [0,0,0], [0, 13.5, -17.24], 'HFS', 'TRA', [3,3,2], sCoeffY, cCoeffX, cCoeffZ, r0);
            MRecon.AddPath;
            data = MReconImpl.GeometryCorrectionUserDef(data, angulation_midslice, offcentre_midslice, patient_position, orientation, resolutionMPS, sCoeffY, cCoeffX, cCoeffZ, r0);                                                   
        end
        function [sCoeffY, cCoeffX, cCoeffZ, r0] = GetGeoCorrCoeffsFromSinfile(sinfile)
            % GetGeoCorrCoeffsFromSinfile: Reads the coefficients needed for a geometry correciton from the sinfile
            %
            % Syntax:     [sCoeffY, cCoeffX, cCoeffZ, r0] = GetGeoCorrCoeffsFromSinfile(sinfile)
            %                        
            % Examples:            
            %       sinfile = 'MySinfile.sin';
            %       [sCoeffY, cCoeffX, cCoeffZ, r0] = MRecon.GetGeoCorrCoeffsFromSinfile(sinfile);
            MRecon.AddPath;
            [sCoeffY, cCoeffX, cCoeffZ, r0] = MReconImpl.GetGeoCorrCoeffsFromSinfile(sinfile);                        
        end
    end       
end
