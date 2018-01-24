function vRGBA = convertRGB2RGBA( R,G,B,A )
%CONVERTRGB2RGBA Converts a 3 RGB plus 1 alpha color descriptor to a single
%RGBA value
% EHarry Nov 2011

vRGBA = [R, G, B, A];
vRGBA = round(vRGBA * 255); % need integer values scaled to range 0-255
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]);

end

