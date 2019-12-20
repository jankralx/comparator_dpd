function coefs = DPD_Default_Coefs(coef1, par)

switch par.dpd.model
    case {'MP', 'MP-Timed'}
        U = MP_Matrix(zeros(100,1), par.dpd.P, par.dpd.M);
    case {'DDR2', 'DDR2-Timed'}
        U = DDR2_Matrix(zeros(100,1), par.dpd.P, par.dpd.M);
    otherwise
        error('%s is unsupported DPD model', par.dpd.model);
end

coefs = zeros(size(U,2),1);
coefs(1) = coef1;