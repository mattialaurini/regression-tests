function catalog = nist_problem_catalog()
%NIST_PROBLEM_CATALOG Metadata for the separable NIST StRD pilot suite.
%
% Official NIST data, starting points and certified values are preserved.
% The bounded domains are study-defined because NIST StRD does not prescribe
% global parameter boxes. The default profile is named "study".

catalog = [ ...
    localBoxBOD(); ...
    localMGH17(); ...
    localLanczos2(); ...
    localKirby2(); ...
    localHahn1(); ...
    localGauss1(); ...
    localNelson(); ...
    localMGH10()];
end

function p = localBase(id,nObs,nPred,modelClass,difficulty,dataUrl,pageUrl)
p = struct();
p.id = string(id);
p.name = string(id);
p.nObservations = nObs;
p.nPredictors = nPred;
p.modelClass = string(modelClass);
p.difficulty = string(difficulty);
p.dataUrl = string(dataUrl);
p.pageUrl = string(pageUrl);
p.dataFile = string(id)+".dat";
p.source = "NIST Statistical Reference Datasets (StRD), SRD 140";
p.domainNote = ["The parameter domains in this package are study-defined " + ...
    "bounded domains; NIST supplies data, starts and certified solutions, " + ...
    "but not canonical global-optimization boxes."];
end

function b = bounds(linearLower,linearUpper,zLower,zUpper)
b = struct("linearLower",double(linearLower(:)), ...
    "linearUpper",double(linearUpper(:)), ...
    "zLower",double(zLower(:)),"zUpper",double(zUpper(:)));
end

function p = localBoxBOD()
p = localBase("BoxBOD",6,1,"Exponential","High", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/BoxBOD.dat", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/boxbod.shtml");
p.modelId = "boxbod";
p.linearIndices = 1;
p.nonlinearIndices = 2;
p.start1 = [1;1];
p.start2 = [100;0.75];
p.certifiedParameters = [2.1380940889e2;5.4723748542e-1];
p.certifiedRSS = 1.1680088766e3;
p.studyBounds = bounds(0,400,0.05,2.0);
p.formula = "b1*(1-exp(-b2*x))";
end

function p = localMGH17()
p = localBase("MGH17",33,1,"Exponential","Average", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/MGH17.dat", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/mgh17.shtml");
p.modelId = "mgh17";
p.linearIndices = [1 2 3];
p.nonlinearIndices = [4 5];
p.start1 = [50;150;-100;1;2];
p.start2 = [0.5;1.5;-1;0.01;0.02];
p.certifiedParameters = [3.7541005211e-1;1.9358469127;-1.4646871366; ...
    1.2867534640e-2;2.2122699662e-2];
p.certifiedRSS = 5.4648946975e-5;
% Ordered, disjoint decay-rate boxes remove the permutation symmetry.
p.studyBounds = bounds([0;0;-3],[1;3;0],[0.005;0.018],[0.018;0.035]);
p.formula = "b1+b2*exp(-b4*x)+b3*exp(-b5*x)";
end

function p = localLanczos2()
p = localBase("Lanczos2",24,1,"Exponential","Average", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Lanczos2.dat", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/lanczos2.shtml");
p.modelId = "lanczos2";
p.linearIndices = [1 3 5];
p.nonlinearIndices = [2 4 6];
p.start1 = [1.2;0.3;5.6;5.5;6.5;7.6];
p.start2 = [0.5;0.7;3.6;4.2;4.0;6.3];
p.certifiedParameters = [9.6251029939e-2;1.0057332849;8.6424689056e-1; ...
    3.0078283915;1.5529016879;5.0028798100];
p.certifiedRSS = 2.2299428125e-11;
% Ordered, disjoint rate boxes remove exponential-component permutations.
p.studyBounds = bounds([0;0;0],[1;5;6],[0.5;2.0;4.0],[1.5;4.5;7.0]);
p.formula = "b1*exp(-b2*x)+b3*exp(-b4*x)+b5*exp(-b6*x)";
end

function p = localKirby2()
p = localBase("Kirby2",151,1,"Rational","Average", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Kirby2.dat", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/kirby2.shtml");
p.modelId = "kirby2";
p.linearIndices = [1 2 3];
p.nonlinearIndices = [4 5];
p.start1 = [2;-0.1;0.003;-0.001;1e-5];
p.start2 = [1.5;-0.15;0.0025;-0.0015;2e-5];
p.certifiedParameters = [1.6745063063;-1.3927397867e-1;2.5961181191e-3; ...
    -1.7241811870e-3;2.1664802578e-5];
p.certifiedRSS = 3.9050739624;
p.studyBounds = bounds([1;-0.2;0.001],[2.5;-0.05;0.004], ...
    [-0.0025;5e-6],[-0.0005;3.5e-5]);
p.formula = "(b1+b2*x+b3*x^2)/(1+b4*x+b5*x^2)";
end

function p = localHahn1()
p = localBase("Hahn1",236,1,"Rational","Average", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Hahn1.dat", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/hahn1.shtml");
p.modelId = "hahn1";
p.linearIndices = [1 2 3 4];
p.nonlinearIndices = [5 6 7];
p.start1 = [10;-1;0.05;-1e-5;-0.05;1e-3;-1e-6];
p.start2 = [1;-0.1;0.005;-1e-6;-0.005;1e-4;-1e-7];
p.certifiedParameters = [1.0776351733;-1.2269296921e-1;4.0863750610e-3; ...
    -1.4262662514e-6;-5.7609940901e-3;2.4053735503e-4;-1.2314450199e-7];
p.certifiedRSS = 1.5324382854;
% A certified-solution neighborhood chosen to keep the denominator positive
% over the complete data range.
p.studyBounds = bounds([0.5;-0.2;0.002;-3e-6], ...
    [1.5;-0.05;0.007;-5e-7], ...
    [-0.0062;2.1e-4;-1.4e-7],[-0.0053;2.8e-4;-1.0e-7]);
p.formula = "(b1+b2*x+b3*x^2+b4*x^3)/(1+b5*x+b6*x^2+b7*x^3)";
end

function p = localGauss1()
p = localBase("Gauss1",250,1,"Exponential/Gaussian","Low", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Gauss1.dat", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/gauss1.shtml");
p.modelId = "gauss1";
p.linearIndices = [1 3 6];
p.nonlinearIndices = [2 4 5 7 8];
p.start1 = [97;0.009;100;65;20;70;178;16.5];
p.start2 = [94;0.0105;99;63;25;71;180;20];
p.certifiedParameters = [98.778210871;1.0497276517e-2;100.48990633; ...
    67.481111276;23.129773360;71.994503004;178.99805021;18.389389025];
p.certifiedRSS = 1.3158222432e3;
p.studyBounds = bounds([80;80;50],[120;120;90], ...
    [0.008;60;15;170;12],[0.013;75;30;190;25]);
p.formula = "b1*exp(-b2*x)+b3*exp(-((x-b4)/b5)^2)+b6*exp(-((x-b7)/b8)^2)";
end

function p = localNelson()
p = localBase("Nelson",128,2,"Exponential","Average", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Nelson.dat", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/nelson.shtml");
p.modelId = "nelson";
p.linearIndices = [1 2];
p.nonlinearIndices = 3;
p.start1 = [2;1e-4;-0.01];
p.start2 = [2.5;5e-9;-0.05];
p.certifiedParameters = [2.5906836021;5.6177717026e-9;-5.7701013174e-2];
p.certifiedRSS = 3.7976833176;
p.studyBounds = bounds([2;0],[3;2e-8],-0.08,-0.03);
p.formula = "log(y)=b1-b2*x1*exp(-b3*x2)";
end

function p = localMGH10()
p = localBase("MGH10",16,1,"Exponential","High", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/MGH10.dat", ...
    "https://www.itl.nist.gov/div898/strd/nls/data/mgh10.shtml");
p.modelId = "mgh10";
p.linearIndices = 1;
p.nonlinearIndices = [2 3];
p.start1 = [2;4e5;2.5e4];
p.start2 = [0.02;4e3;250];
p.certifiedParameters = [5.6096364710e-3;6.1813463463e3;3.4522363462e2];
p.certifiedRSS = 8.7945855171e1;
p.studyBounds = bounds(0.001,0.05,[3e3;200],[9e3;500]);
p.formula = "b1*exp(b2/(x+b3))";
end
