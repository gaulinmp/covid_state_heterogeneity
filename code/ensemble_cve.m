


function output = ensemble_cve(wgt, ensembledata)

%% Cross-validation errors for ensemble

Tfor =ensembledata.Tfor_max - ensembledata.Tfor_min+1;

ensemblepred_O = sum(kron(wgt',ones(1,Tfor)).*ensembledata.O_mat);
CVE_O = sum( (ensemblepred_O-ensembledata.Odat).^2 )/(Tfor);

ensemblepred_F = sum(kron(wgt',ones(1,Tfor)).*ensembledata.F_mat);
CVE_F = sum( (ensemblepred_F-ensembledata.Fdat).^2 )/(Tfor);

%ensemblepred_a = sum(kron(wgt',ones(1,11)).*ensembledata.a_mat);
%CVE_a = sum(ensemblepred_a-ensembledata.adat)/(ensembledata.Tfor_max-ensembledata.Tfor_min+1);


output = CVE_O + CVE_F;

