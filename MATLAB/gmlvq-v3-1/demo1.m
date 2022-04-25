load samplesData/uci-segmentation-sampled.mat
gmlvq = GMLVQ.GMLVQ(fvec, lbl, GMLVQ.Parameters(), 50);
result = gmlvq.runSingle();
result.plot();
