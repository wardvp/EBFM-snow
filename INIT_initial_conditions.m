%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIAL CONDITIONS:
%%% Read initial conditions at starting time from restart file or set 
%%% manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT,IN,OUTFILE] = INIT_initial_conditions(C,grid,io)

OUT = struct;                                                               % structure containing model output variables
IN = struct;                                                                % structure containing model input variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET INITIAL CONDITIONS AT START OF SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start from restart file
if (io.readbootfile)
    disp('Initialize from restart file...');
    cd(io.rebootdir);
    load(io.bootfilein);
    OUT = boot;
    cd(io.homedir);

% Start from preset initial conditions
else
    disp('Initialize from manually set conditions...');
    OUT.Tsurf(1:grid.gpsum,1)                     = 273.15;                 % surface temperature (K)
    OUT.subT(1:grid.gpsum,1:grid.nl)              = 265.0;                  % vertical temperatures (K)
    OUT.subW(1:grid.gpsum,1:grid.nl)              = 0.0;                    % vertical irreducible water content (kg)
    OUT.subS(1:grid.gpsum,1:grid.nl)              = 0.0;                    % vertical slush water content (kg)
    OUT.subD(1:grid.gpsum,1:grid.nl)              = C.soil_density;         % vertical densities (kg m-3)
    OUT.subSOIL(1:grid.gpsum,1:grid.nl)           = 1;                      % vertical snow and soil distribution (0 = snow, 1 = soil)
    OUT.subTmean                                  = OUT.subT;               % vertical annual mean layer temperature (K)
    OUT.timelastsnow(1:grid.gpsum,1)              = 1.0;                    % timestep of last snow fall (days)
    OUT.ys(1:grid.gpsum,1)                        = 500.0;                  % annual snow fall (mm w.e.)
    OUT.subZ(1:grid.gpsum,1:grid.nl)              = grid.max_subZ;          % vertical layer thickness (m)
    OUT.alb_snow(1:grid.gpsum,1)                  = C.alb_fresh;            % snow albedo
    OUT.snowmass(1:grid.gpsum,1)                  = 0.0;                    % cumulative snow mass balance (m w.e.)
end

% Declare non-initialized variables
OUT.cmb_cumulative(1:grid.gpsum,1)                  = 0.0;                  % cumulative climatic mass balance (m w.e.)
OUT.cmb(1:grid.gpsum,1)                             = 0.0;                  % climatic mass balance (m w.e.)
OUT.subK(1:grid.gpsum,1:grid.nl)                    = 0.0;                  % vertical conductivity (m2 s-1)
OUT.subCeff(1:grid.gpsum,1:grid.nl)                 = 0.0;                  % vertical effective heat capacity (J m-3 K)
OUT.subWvol(1:grid.gpsum,1:grid.nl)                 = 0.0;                  % vertical volumetric water content (fraction)
OUT.surfH(1:grid.gpsum,1)                           = 0.0;                  % surface height (m)
OUT.Dfreshsnow(1:grid.gpsum,1)                      = 0.0;                  % fresh snow density (kg m-3)
OUT.tstar(1:grid.gpsum,1)                           = 0.0;                  % albedo decay timescale (days)
OUT.runoff_irr_deep_mean(1:grid.gpsum,1)            = 0.0;                  % runoff of irreducible water below base of the model (m w.e.)

IN.T(1:grid.gpsum,1)                                = 0.0;                  % air temperature (K)
IN.P(1:grid.gpsum,1)                                = 0.0;                  % precipitation (m w.e.)
IN.snow(1:grid.gpsum,1)                             = 0.0;                  % snow precip (m w.e.)
IN.rain(1:grid.gpsum)                               = 0.0;                  % rain precip (m w.e.)
IN.yearsnow(1:grid.gpsum,1:grid.nl)                 = 0.0;                  % annual snow precip (m w.e.)
IN.logyearsnow(1:grid.gpsum,1:grid.nl)              = 0.0;                  % annual snow precip (log m w.e.)
IN.C(1:grid.gpsum,1)                                = 0.0;                  % cloud cover (fraction)
IN.WS(1:grid.gpsum,1)                               = 0.0;                  % wind speed (m s-1)
IN.RH(1:grid.gpsum,1)                               = 0.0;                  % relative humidity (fraction)
IN.q(1:grid.gpsum,1)                                = 0.0;                  % specific humidity (g kg-1)
IN.VP(1:grid.gpsum,1)                               = 0.0;                  % vapor pressure (Pa)
IN.Dair(1:grid.gpsum,1)                             = 0.0;                  % air density (kg m-3)
IN.Pres(1:grid.gpsum,1)                             = 0.0;                  % air pressure (Pa)

OUTFILE = struct;                                                           % output to be saved to files

end