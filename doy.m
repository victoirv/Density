function DOY = doy(Y,M,D,H,MIN,S)
%DOY Specifies the day of year given year, month, day.
%
%   doy = DOY(Y,M,D) or doy = DOY([Y,M,D]) (faster) returns
%   the corresponding integer day of year. 
%
%   doy = DOY(Y,M,D,H,MIN,S) or doy = DOY([Y,M,D,H,MIN,S]) returns
%   the fractional day of year.  
%
%   Other usage assumes unspecified inputs are zero:
%
%   DOY(Y,M,D,H)       or DOY([Y,M,D,H])
%   DOY(Y,M,D,H,MIN)   or DOY([Y,M,D,H,MIN])
%   DOY(Y,M,D,H,MIN,S) or DOY([Y,M,D,H,MIN,S])
%
%   See also YMD, DAY_NUMBER, NUMBER_DAYS, DAYS_PER_YEAR, TIME_UNION,
%   IS_LEAP_YEAR.

% R.S. Weigel, ymd(120,2004).  

if (nargin >=3 )
  if (nargin == 3)
    Y = [Y,M,D]; clear M D
  end
  if (nargin == 4)
    Y = [Y,M,D,H]; clear M D H
  end
  if (nargin == 5)
    Y = [Y,M,D,H,MIN]; clear M D H MIN
  end
  if (nargin == 6)
    Y = [Y,M,D,H,MIN,S]; clear M D H MIN S
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If vector input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin == 1)
  if ( size(Y,2) < 3 ) | ( size(Y,2) > 6 )
    error('Input must have 3-6 columns.');
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sanity checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (max(Y(:,2)) > 12) | (max(Y(:,2)) < 1)
  error('Month must be in range 1-12');
end
if (max(Y(:,3)) < 1)
  error('Day must be greater than 1');
end
%temp = find ( (days_per_month(M,Y(:,1))-D)  < 0 );
temp = [];
if (~isempty(temp)) 
  estr = sprintf('At least one set of Y,M,D values is not consistent\n');
  estr = [estr,'with the number of days in M.  First offender:'];  
  [Y(temp,1),Y(temp,2),Y(temp,3)]
  error('Aborting due to above problem.');
end

if (size(Y,2) > 3)
  if (Y(:,4) > 23) | (Y(:,4) < 0)
    error('Hour must be in range 0-23');
  end
end
if (size(Y,2) > 4)
  if (Y(:,5) > 59) | (Y(:,5) < 0)
    error('Minute must be in range 0-59');
  end
end
if (size(Y,2) > 5)
  if (Y(:,6) > 59) | (Y(:,6) < 0)
    error('Second must be in range 0-59');
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I            = 1 + is_leap_year(Y(:,1));
day_sum(:,1) = [0;31;59;90;120;151;181;212;243;273;304;334];
day_sum(:,2) = [0;31;60;91;121;152;182;213;244;274;305;335];
DOY          = day_sum(sub2ind(size(day_sum),Y(:,2),I)) + Y(:,3);

if (size(Y,2) == 4)
  DOY = DOY + Y(:,4)/24;
end
if (size(Y,2) == 5)
  DOY = DOY + (Y(:,4)*60 + Y(:,5))/1440;
end
if (size(Y,2) == 6)
  DOY = DOY + (Y(:,4)*3600 + Y(:,5)*60 + Y(:,6))/86400;
end
