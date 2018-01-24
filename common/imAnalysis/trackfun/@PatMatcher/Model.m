function pm=Model(pm)
%Models the deformation of the template

if(isempty(pm.model.params))
   pm.model.params=zeros(pm.dims,1);
   pm.model.params(pm.dims+1)=1;
end;

pm.A=zeros(pm.dims);
maskA = ones(pm.dims);
maskT = ones(pm.dims,1);

switch(pm.cmd.params)
	case 'all'
   case 'shift'
      maskA = zeros(pm.dims);
   case 'shiftscale'
      maskA=eye(pm.dims);      
end;

for dim = 1:pm.dims,
   pm.model.coord(dim,:)=(pm.template.coord(dim,:)-pm.patCenter(dim))*pm.model.params(pm.dims+1)...
      +pm.model.params(dim);
end;
pm=ModelGrad(pm,0,[],0);
%for l = 1:prod(pm.template.size)
%   pm.model.grad(:,:,l)=eye(pm.dims);
%end;
