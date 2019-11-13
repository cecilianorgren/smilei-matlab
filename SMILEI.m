 classdef SMILEI
  %TSeries Load SMILEI simulation data
  %   Does not contain the data, but loads it in an easily accesible manner
  %
  %   SM = SMILEI(h5FilePath,nameListFilePath)
  %   Bx = SM.Bx; % Bx is a (nt x nx x ny) matrix
  %   B = SM.B; % structure with 3 (nt x nx x ny) matrices  
  
  properties (Access = protected) % can only be set from within the class (?)
    % Data can be arbitrary size, so the class contains a pointer to the 
    % data file and each time loads the data with. This is actualy really
    % bad, because it removes all the advantages of having class-specific
    % functions that easily acecsses the data and performs operarations. 
    file_
    namelist_
    info_
    fields_
    iteration_
    twpe_
    twci_
    gridsize_ % this should not be here, just put the grid instead, and use a function to get the size, but I haven't read the namelist yet
    grid_
    wpewce_ = [];
    mime_ = [];
    coordinateSystem_ = '';
  end
  
  properties (Dependent = true)
    file
    namelist
    info
    fields
    iteration
    twpe
    twci
    gridsize
    grid
    wpewce
    mime
    coordinateSystem
  end
  
  properties (Constant = true)    
  end
  
  properties (Constant = true, Hidden = true)
    MAX_TENSOR_ORDER = 2;
    BASIS = {'xyz','xzy'}; % in order to use, need to define transformations between these
    BASIS_NAMES = {'smilei','michael'};
  end
  
  %properties (SetAccess = protected)
  %  representation % can be related to basis
  %end
  
  properties
    userData = []; % anything can be added here
  end
  
  methods
    function obj = SMILEI(h5filePath,nameListFilePath)
      if exist('nameListFilePath','var') && not(isempty(nameListFilePath))
        [wpewce,mime,dxyz] = read_filelist(nameListFilePath);
        obj.wpewce_ = wpewce;
        obj.mime_ = mime;
      end
      
      obj.info = h5info(h5filePath);
      obj.file = h5filePath; 
      obj.iteration = get_iterations(obj);
      %obj.twci = get_iterations(obj);
      obj.twpe = get_twpe(obj);
      obj.twci = get_twpe(obj)*0;
      obj.fields_ = get_fields(obj);
      obj.gridsize_ = get_gridsize(obj);
      %obj.t_ = get_time(h5filePath);
      
      
    
    %  obj.fullDim_ = cell(ndims(data)-1,1); % time should not be included -> -1
    %  obj.representation = cell(ndims(data)-1,1); iDim = 1; % time should not be included -> -1
    end
    
    function [varargout] = subsref(obj,idx)
      %SUBSREF handle indexing
      nargout
      switch idx(1).type
        % Use the built-in subsref for dot notation
        case '.'
          [varargout{1:nargout}] = builtin('subsref',obj,idx);
        case '()'          
          nargout
          obj.iteration_ = builtin('subsref',obj.iteration,idx(1));
          obj.twpe_ = builtin('subsref',obj.twpe,idx(1));
          obj.twci_ = builtin('subsref',obj.twci,idx(1)); 
          if numel(idx) > 1
            obj = builtin('subsref',obj,idx(2:end));
          end          
          
          if isa(obj,'SMILEI') % return smilei object
            [varargout{1:nargout}] = obj;
          elseif isa(obj,'numeric') && numel(obj) == nargout % return subreffed data, for example twpe         
            varargout{1} = obj;
            varargout{2} = obj;
            %nargout == 1;
            %for iarg = 1:nargout            
            %  varargout{iarg} = obj(iarg);
            %end
          end          
          
          
        case '{}'
          error('SMILEI:subsref',...
            'Not a supported subscripted reference.')
      end
    end
    
    function [wpewce,memi,dxyz] = read_namelist(filepath)
      
    end
    function value = get.coordinateSystem(obj)
      value = obj.coordinateSystem_;
    end
        
    function value = length(obj)
      value = numel(obj.iteration);
    end
   
    % Get subset of data, time and/or space
    function obj = ilim(obj,value)
      % Get 
      obj = subsref(obj,value);
      
    end
    
    % Get simulation meta data and parameters
    function out = get_twpe(obj)
      fileInfo = obj.info_;
      nOutput = numel(fileInfo.Groups.Groups);

      for iOutput = 1:nOutput
        time(iOutput) = fileInfo.Groups.Groups(iOutput).Attributes(1).Value;
      end
      out = time;
    end
    function out = get_iterations(obj)
      fileInfo = obj.info_;
      nOutput = numel(fileInfo.Groups.Groups);
      for iOutput = 1:nOutput
        str = fileInfo.Groups.Groups(iOutput).Name;
        split_str = strsplit(str,'/');
        iterations(iOutput) = str2num(split_str{3});
      end

      out = iterations;  
    end
    function out = get_fields(obj)
      fileInfo = obj.info;
      % fields structure is the same for all times
      out = {fileInfo.Groups(1).Groups(1).Datasets.Name};
    end
    function out = get_gridsize(obj)
      fileInfo = obj.info_;
      out = fileInfo.Groups(1).Groups(1).Datasets(1).Dataspace.Size;
    end
    
    % Get fields
    function out = Bx(obj)
      out = get_field(obj,'Bx');
    end
    function out = By(obj)
      out = get_field(obj,'By');
    end
    function out = Bz(obj)
      out = get_field(obj,'Bz');
    end
    function out = Ex(obj)
      out = get_field(obj,'Ex');
    end
    function out = Ey(obj)
      out = get_field(obj,'Ey');
    end
    function out = Ez(obj)
      out = get_field(obj,'Ez');
    end
    function out = Jx_eon(obj)
      out = get_field(obj,'Jx_eon');
    end
    function out = Jy_eon(obj)
      out = get_field(obj,'Jy_eon');
    end
    function out = Jz_eon(obj)
      out = get_field(obj,'Jz_eon');
    end
    function out = Jx_ion(obj)
      out = get_field(obj,'Jx_ion');
    end
    function out = Jy_ion(obj)
      out = get_field(obj,'Jy_ion');
    end
    function out = Jz_ion(obj)
      out = get_field(obj,'Jz_ion');
    end
    function out = Jx(obj)
      % find fields that contains Jx
      index = find(contains(obj.fields,'Jx'));
      Jx = zeros([obj.length,obj.gridsize]);
      for iField = 1:numel(index)
        Jx = Jx + obj.get_field(obj.fields{index(iField)});
      end
      out = Jx;
    end
    function out = Jy(obj)
      % find fields that contains Jx
      index = find(contains(obj.fields,'Jy'));
      J = zeros([obj.length,obj.gridsize]);
      for iField = 1:numel(index)
        J = J + obj.get_field(obj.fields{index(iField)});
      end
      out = J;
    end
    function out = Jz(obj)
      % find fields that contains Jx
      index = find(contains(obj.fields,'Jz'));
      J = zeros([obj.length,obj.gridsize]);
      for iField = 1:numel(index)
        J = J + obj.get_field(obj.fields{index(iField)});
      end
      out = J;
    end
    function out = Rho_eon(obj)
      out = get_field(obj,'Rho_eon');
    end
    function out = Rho_ion(obj)
      out = get_field(obj,'Rho_eon');
    end
    function out = Rho(obj)
      % find fields that contains Rho
      index = find(contains(obj.fields,'Rho'));
      Rho = zeros([obj.length,obj.gridsize]);
      for iField = 1:numel(index)
        Rho = Rho + obj.get_field(obj.fields{index(iField)});
      end
      out = Rho;
    end
    
    % Get and set properties
    function obj = set.coordinateSystem(obj,value)
      if obj.tensorOrder_ < 1 
        error('irf:TSeries:setcoordinateSystem:badInputs',...
          'coordinateSystem can only be set for a tensor')
      end
      if ~ischar(value)
        error('irf:TSeries:setcoordinateSystem:badInputs',...
          'expecting string input')
      end
      obj.coordinateSystem_ = value;
    end
    
    function obj = set.info(obj,value)
      obj.info_ = value;
    end
    function obj = set.file(obj,value)
      obj.file_ = value;
    end
    function obj = set.namelist(obj,value)
      obj.namelist_ = value;
    end
    function obj = set.fields(obj,value)
      obj.fields_ = value;
    end
    function obj = set.iteration(obj,value)
      obj.iteration_ = value;
    end
    function obj = set.twpe(obj,value)
      obj.twpe_ = value;
    end
    function obj = set.twci(obj,value)
      obj.twci_ = value;
    end
    function obj = set.gridsize(obj,value)
      obj.gridsize_ = value;
    end
    function obj = set.grid(obj,value)
      obj.grid_ = value;
    end
    
    function value = get.info(obj)
      value = obj.info_;
    end
    function value = get.file(obj)
      value = obj.file_;
    end
    function value = get.namelist(obj)
      value = obj.namelist_;
    end
    function value = get.fields(obj)
      value = obj.fields_;
    end   
    function value = get.iteration(obj)
      value = obj.iteration_;
    end 
    function value = get.twpe(obj)
      value = obj.twpe_;
    end
    function value = get.twci(obj)
      value = obj.twci_;
    end
    function value = get.gridsize(obj)
      value = obj.gridsize_;
    end
    function value = get.grid(obj)
      value = obj.grid_;
    end            
    
  end
  
  methods (Access=protected)
    function out = get_field(obj,field)
      % valiadate field       
      if not(ismember(field,obj.fields))
        error(sprintf('Field ''%s'' not recognized.',field))
      end
      
      % get iterations
      iterations = obj.iteration;
      nIter = obj.length;
      % initialize matrix
      data = nan([nIter,obj.gridsize]);
      for iIter = 1:nIter
        iter = iterations(iIter);
        str_iter = sprintf('%010.0f',iter);
        data_tmp = h5read(obj.file,['/data/' str_iter '/' field]); % de;
        data(iIter,:,:) = data_tmp;
      end
      out = data;
    end
    
    function Ts = changeBasis(obj, flag)
      % Tranform from one coordinate system to another and return new
      % TimeSeries.
      % flag: = 'xyz>rlp' - Cartesian XYZ to spherical latitude
      %         'rlp>xyz' - Spherical latitude to cartesian XYZ
      %         'xyz>rpz' - Cartesian XYZ to cylindrical
      %         'rpz>xyz' - Cylidrical to cartesian XYZ
      %         'xyz>rtp' - Cartesian XYZ to spherical colatitude
      %         'rtp>xyz' - Spherical colatitude to cartesian XYZ
      %         'rtp>rlp' - Spherical colatitude to spherical latitude
      %         'rlp>rtp' - Spherical latitude to colatitude
      switch lower(flag)
        case 'xyz>rlp'
          [phi, lambda, r] = cart2sph(obj.x.data, obj.y.data, obj.z.data);
          Ts = TSeries(obj.time, [r, lambda, phi], 'vec_rlp');
        case 'rlp>xyz'
          [x, y, z] = sph2cart(obj.phi.data, obj.lambda.data, obj.r.data);
          Ts = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'xyz>rpz'
          [phi, r, z] = cart2pol(obj.x.data, obj.y.data, obj.z.data);
          Ts = TSeries(obj.time, [r, phi, z], 'vec_rpz');
        case 'rpz>xyz'
          [x, y, z] = pol2cart(obj.phi.data, obj.r.data, obj.z.data);
          Ts = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'xyz>rtp'
          [phi, lambda, r] = cart2sph(obj.x.data, obj.y.data, obj.z.data);
          theta = pi/2 - lambda;
          Ts = TSeries(obj.time, [r, theta, phi], 'vec_rtp');
        case 'rtp>xyz'
          lambda = pi/2 - obj.theta.data;
          [x, y, z] = sph2cart(obj.phi.data, lambda, obj.r.data);
          Ts = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'rtp>rlp'
          lambda = pi/2 - obj.theta.data;
          Ts = TSeries(obj.time, [obj.r.data,lambda,obj.phi.data],'vec_rlp');
        case 'rlp>rtp'
          theta = pi/2 - obj.lambda.data;
          Ts = TSeries(obj.time, [obj.r.data,theta,obj.phi.data],'vec_rtp');
        case 'xy>rp'
          [phi, r] = cart2pol(obj.x.data, obj.y.data);
          Ts = TSeries(obj.time, [r, phi], 'vec_rp');
        case 'rp>xy'
          [x, y] = pol2cart(obj.phi.data, obj.r.data);
          Ts = TSeries(obj.time, [x, y], 'vec_xy');
        otherwise
          errStr='Invalid transformation'; error(errStr);
      end
    end
  end
  
end

