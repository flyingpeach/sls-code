function covMtx = get_cov_mtx(satParams, eyeSize)

N_SAMPLES = 1e5;
eta       = satParams.eta_;

normDist      = makedist('Normal','mu', 0,'sigma', satParams.noiseStd_);
normDistTrunc = truncate(normDist, -eta(2), eta(2));
samples       = random(normDistTrunc, [N_SAMPLES, 1]);

p1 = sat(samples, eta(1));
p2 = sat(samples, eta(2));

myEye  = eye(eyeSize);

a1     = p1 .* p1;
a2     = p1 .* (p2 - p1);
a3     = (p2 - p1) .* (p2 - p1);

covMtx = [mean(a1)*myEye  mean(a2)*myEye; 
          mean(a2)*myEye  mean(a3)*myEye]; 