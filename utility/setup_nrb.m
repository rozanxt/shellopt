function nrb = setup_nrb(dsn, nrb)
    nrb.coefs(1:3, :, :) = reshape(dsn, 3, size(nrb.coefs, 2), size(nrb.coefs, 3));
end
