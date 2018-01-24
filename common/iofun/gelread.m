function image = gelread(filename)
%GELREAD Read a GEL or 16-bit TIFF file from disk.
%	IM = GELREAD('filename')
%	The TIFF file must have only a single strip which may be
%       stored as either 8 bit or 16 bit and must be uncompressed.
%	Only the image data itself is read from the file.  All other
%	information in the file is ignored.  Image in MATLAB will be
%	double-precision values from 1.0 to 256.0.


bits = 1;
comp = 1;
res_unit = 2;
c_map = [];

if isempty(findstr(filename,'.'))
	filename=[filename,'.tif'];
end;

[file, message] = fopen(filename,'r','l');
if file == -1
  a=['file ',filename,' not found.'];
  error(a);
end

% read header
% read byte order: II = little endian, MM = big endian
byte_order = setstr(fread(file,2,'uchar'));
if ((~strcmp(byte_order','II')) & (~strcmp(byte_order','MM')))
  error('This is not a GEL or TIFF file.');
end

% read in a number which identifies file as TIFF format
if (strcmp(byte_order','II'))
  tiff_id = fread(file,1,'ushort','l');
else
  tiff_id = fread(file,1,'ushort','b');
end

% read the byte offset for the first image file directory (IFD)
if (strcmp(byte_order','II'))
  ifd_pos = fread(file,1,'ulong','l');
else
  ifd_pos = fread(file,1,'ulong','b');
end

% move in the file to the first IFD
fseek(file,ifd_pos,-1);

%read IFD
%read in the number of IFD entries
if (strcmp(byte_order','II'))
  num_entries = fread(file,1,'ushort','l');
else 
  num_entries = fread(file,1,'ushort','b');
end

%read and process each IFD entry
for i = 1:num_entries
  % save the current position in the file
  file_pos = ftell(file);

  % read entry tag
  if (strcmp(byte_order','II'))
    tag = fread(file,1,'ushort','l');
  else
    tag = fread(file,1,'ushort','b');
  end

  % image width - number of columns
  if (tag == 256)
    if (strcmp(byte_order','II'))
      type = fread(file,1,'ushort','l');
      fseek(file,4,0);
      if (type == 3)
        width = fread(file,1,'ushort','l');
      else
        width = fread(file,1,'ulong','l');
      end
    else
      type = fread(file,1,'ushort','b');
      fseek(file,4,0);
      if (type == 3)
        width = fread(file,1,'ushort','b');
      else
        width = fread(file,1,'ulong','b');
      end
    end

  % image height - number of rows
  elseif (tag == 257)
    if (strcmp(byte_order','II'))
      type = fread(file,1,'ushort','l');
      fseek(file,4,0);
      if (type == 3)
        height = fread(file,1,'ushort','l');
      else
        height = fread(file,1,'ulong','l');
      end
    else
      type = fread(file,1,'ushort','b');
      fseek(file,4,0);
      if (type == 3)
        height = fread(file,1,'ushort','b');
      else
        height = fread(file,1,'ulong','b');
      end
    end

  % bits per sample
  elseif (tag == 258)
    if (strcmp(byte_order','II'))
      type = fread(file,1,'ushort','l');
      count = fread(file,1,'ulong','l');
      if (2*count > 4)
        % next field contains an offset
        val_offset = fread(file,1,'ulong','l');
        fseek(file,val_offset,-1);
        bits = fread(file,count,'ushort','l');
      else
        % next field contains values
        bits = fread(file,count,'ushort','l');
      end
    else
      type = fread(file,1,'ushort','b');
      count = fread(file,1,'ulong','b');
      if (2*count > 4)
        % next field contains an offset
        val_offset = fread(file,1,'ulong','b');
        fseek(file,val_offset,-1);
        bits = fread(file,count,'ushort','b');
      else
        % next field contains values
        bits = fread(file,count,'ushort','b');
      end
    end
    num_samples = count;
    if num_samples==3, % RGB image
        r = 24;
    elseif bits(1)==1, % binary image
        r = 1;
    elseif bits(1)==4, % 16 color indexed image
        r = 4;
    elseif bits(1)==8
        r = 8;
    elseif bits(1)==16
        r = 16;
    else
        error('Unsupported number of image depth');
    end

  % photometric interpretation
  elseif (tag == 262)
    if (strcmp(byte_order','II')) 
      fseek(file,6,0);
      photo_type = fread(file,1,'ushort','l');
    else
      fseek(file,6,0);
      photo_type = fread(file,1,'ushort','b');
    end

  % strip offsets
  elseif (tag == 273)
    if (strcmp(byte_order','II')) 
      type = fread(file,1,'ushort','l');
      count = fread(file,1,'ulong','l');
      num_strips = count;
      if ((type-2)*2*count > 4) 
        % next field contains an offset
        val_offset = fread(file,1,'ulong','l'); 
        fseek(file,val_offset,-1);
        if (type == 3)
          strip_offsets = fread(file,count,'ushort','l');
        else 
          strip_offsets = fread(file,count,'ulong','l');
        end
      else
        % next field contains values 
        if (type == 3)
          strip_offsets = fread(file,count,'ushort','l');
        else 
          strip_offsets = fread(file,count,'ulong','l');
        end
      end
    else
      type = fread(file,1,'ushort','b');
      count = fread(file,1,'ulong','b');
      num_strips = count;
      if ((type-2)*2*count > 4)
        % next field contains an offset
        val_offset = fread(file,1,'ulong','b');
        fseek(file,val_offset,-1); 
        if (type == 3)
          strip_offsets = fread(file,count,'ushort','b');
        else 
          strip_offsets = fread(file,count,'ulong','b');
        end
      else
        % next field contains values 
        if (type == 3)
          strip_offsets = fread(file,count,'ushort','b');
        else 
          strip_offsets = fread(file,count,'ulong','b');
        end
      end
    end

  % rows per strip
  elseif (tag == 278)
    if (strcmp(byte_order','II'))
      type = fread(file,1,'ushort','l');
      fseek(file,4,0);
      if (type == 3)
        rows_per_strip = fread(file,1,'ushort','l');
      else
        rows_per_strip = fread(file,1,'ulong','l');
      end
    else
      type = fread(file,1,'ushort','b');
      fseek(file,4,0);
      if (type == 3)
        rows_per_strip = fread(file,1,'ushort','b');
      else
        rows_per_strip = fread(file,1,'ulong','b');
      end
    end
     
  % strip byte counts - number of bytes in each strip after any compression
  elseif (tag == 279)
    if (strcmp(byte_order','II'))
      type = fread(file,1,'ushort','l');
      count = fread(file,1,'ulong','l');
      if ((type-2)*2*count > 4)
        % next field contains an offset
        val_offset = fread(file,1,'ulong','l');
        fseek(file,val_offset,-1);
        if (type == 3)
          strip_bytes = fread(file,count,'ushort','l');
        else
          strip_bytes = fread(file,count,'ulong','l');
        end
      else
        % next field contains values
        if (type == 3)
          strip_bytes = fread(file,count,'ushort','l');
        else
          strip_bytes = fread(file,count,'ulong','l');
        end
      end
    else
      type = fread(file,1,'ushort','b');
      count = fread(file,1,'ulong','b');
      if ((type-2)*2*count > 4)
        % next field contains an offset
        val_offset = fread(file,1,'ulong','b');
        fseek(file,val_offset,-1);
        if (type == 3)
          strip_bytes = fread(file,count,'ushort','b');
        else
          strip_bytes = fread(file,count,'ulong','b');
        end
      else
        % next field contains values
        if (type == 3)
          strip_bytes = fread(file,count,'ushort','b');
        else
          strip_bytes = fread(file,count,'ulong','b');
        end
      end
    end

  % X resolution
  elseif (tag == 282)
    fseek(file,6,0);
    if (strcmp(byte_order','II'))
      val_offset = fread(file,1,'ulong','l');
      fseek(file,val_offset,-1);
      numerator = fread(file,1,'ulong','l');
      denominator = fread(file,1,'ulong','l');
    else
      val_offset = fread(file,1,'ulong','b');
      fseek(file,val_offset,-1);
      numerator = fread(file,1,'ulong','b');
      denominator = fread(file,1,'ulong','b');
    end
      x_res = numerator / denominator;
   
  % Y resolution          
  elseif (tag == 283)
    fseek(file,6,0);
    if (strcmp(byte_order','II'))
      val_offset = fread(file,1,'ulong','l');
      fseek(file,val_offset,-1);
      numerator = fread(file,1,'ulong','l');
      denominator = fread(file,1,'ulong','l');
    else
      val_offset = fread(file,1,'ulong','b');
      fseek(file,val_offset,-1);
      numerator = fread(file,1,'ulong','b');
      denominator = fread(file,1,'ulong','b');
    end
      y_res = numerator / denominator;

  % resolution unit
  elseif (tag == 296)
    fseek(file,6,0);
    if (strcmp(byte_order','II'))
      res_unit = fread(file,1,'ushort','l');
    else
      res_unit = fread(file,1,'ushort','b');
    end

  end
  % move to next FID entry in the file
  fseek(file,file_pos+12,-1);
end

% compute the width of each row in bytes
width_bytes = ceil(width/(8/bits(1)))*num_samples;

% read the strips, decompress them, and place them in a matrix
if (comp == 1)
  count = 0;
  % Preallocate image matrix.
  if isempty(strip_bytes)
    strip_bytes = (width * height) / (num_strips * (8/bits(1)));
  end
  numRows = width_bytes;
  numCols = sum(strip_bytes/width_bytes);
  image = zeros(numRows/2, numCols);
  colIdx = 1;
  for i = 1:num_strips
    fseek(file,strip_offsets(i),-1);
    if (strcmp(byte_order','II'))
      if bits(1) == 16
        strip = fread(file,strip_bytes(i)/2,'ushort','l');
      else   
        strip = fread(file,strip_bytes(i),'schar','l');
      end
    else
      if bits(1) == 16
        strip = fread(file,strip_bytes(i)/2,'ushort','b');
      else   
        strip = fread(file,strip_bytes(i),'schar','b');
      end
    end
    if length(strip)~=(strip_bytes(i)/2), 
      error('End of file reached unexpectedly.');
    end
    stripCols = strip_bytes(i)/width_bytes;
    image(:, colIdx:(colIdx+stripCols-1)) = ...
	reshape(strip, width_bytes/2, stripCols);
    colIdx = colIdx + stripCols;
  end
else
  error('Compressed file not readable');
end
image=image';
fclose(file);
