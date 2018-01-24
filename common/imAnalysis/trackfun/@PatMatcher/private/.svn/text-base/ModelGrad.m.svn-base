function [pm,pCount]=ModelGrad(pm,pCount,crd,dim)
%Computes the gradient of the deformation of the template
if(dim==pm.dims)
   pCount=pCount+1;
   cellcrd=num2cell(crd);
   ind=sub2ind(pm.template.size,cellcrd{:});
%   ind=0;
%   dSze=prod(pm.template.size);
%   for dim= 1:pm.dims
%      dSze=dSze/pm.template.size(pm.dims+1-dim);
%      ind=ind+(crd(dim)-1)*dSze;
%   end;
   pm.model.grad(:,:,pCount)=[eye(pm.dims) (pm.template.coord(:,ind)-pm.patCenter')];
   return;
end;

dim = dim +1;
for l = 2:(pm.template.size(dim)-1)
   crd(dim)=l;
   [pm pCount]=ModelGrad(pm,pCount,crd,dim);
end;
