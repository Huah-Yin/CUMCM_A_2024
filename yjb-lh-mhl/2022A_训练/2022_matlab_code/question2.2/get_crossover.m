function MyGama = get_crossover(staus,MyGama)

    if  MyGama==0.9 && staus>=20
        MyGama=0.02;
    else
        MyGama = min(MyGama + staus * 0.001, 0.9);
    end
end
