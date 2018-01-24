function pm=Radiometric(pm)
%Models the intenisty changes between patch and template
meanPat=mean(pm.patch(:));
stdPat=std(pm.patch(:));
pm.patch(:)=pm.template.std/stdPat*(pm.patch(:)-meanPat)+pm.template.mean;
%pm.patch=pm.patch.*(pm.patch>=0 & pm.patch<=255);
