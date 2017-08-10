function newdata = avg_fret_data(data,ind)

    indtemp = 1:3;
    indtemp(ind) = [];

    data = data(indtemp,:);

    newdata = squeeze(mean(data));

end