function output = own_downsample(Sig, Sps, Npp)
    ds_Sig = downsample(Sig, Sps);
    output = ds_Sig;
end