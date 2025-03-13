function ind=find_closest_pos(mot_par,iv_cyc,ivf,fcrptd)
    [~,i]=min(squeeze(max(squeeze(abs(mot_par(1:3,iv_cyc,~fcrptd)- ...
                                        mot_par(1:3,iv_cyc,ivf))),[],1)));
    ind=find(~fcrptd);
    ind=ind(i);
end