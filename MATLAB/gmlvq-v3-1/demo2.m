load samplesData/twoclass-simple.mat
gmlvq = GMLVQ.GMLVQ(fvec, lbl, GMLVQ.Parameters('rocClass', 1), 30);
result = gmlvq.runSingle();
plot(result); % Both calling styles work
