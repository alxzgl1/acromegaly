% This function runs a single training run on the given data
% @param this GMLVQ.GMLVQ
% @out res GMVLQ.Result
function res = runSingle(this)

GMLVQ.Helpers.setRNG(this.params.rngseed);

% Construct the Run class
run = GMLVQ.Run(this, this.data, GMLVQ.DataPair.empty);
res = run.execute();

end