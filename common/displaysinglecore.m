function flame = displaysinglecore(core,normm,norms,dens,max,unifo)


flame.cores{1} = core;
flame.normm{1} = normm;
flame.norms{1} = norms;
flame.scale{1} = dens;
flame.maxim{1} = max;
flame.unifo{1} = unifo;

displaysingleflame(flame);


end