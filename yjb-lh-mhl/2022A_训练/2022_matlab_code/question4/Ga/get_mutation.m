function mu = get_mutation(staus_num,mu)
    %变异
   

    if mu==0.9 && staus_num>=20
        mu=0.02;
    else
         mu= min(mu + staus_num * 0.001, 0.9);
    end
end
