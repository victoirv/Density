function paperfigures(figuretype)
%Figuretype: 0=dissertation, 1=paper, 2=presentation
if(isempty(figuretype))
    figuretype=0;
end


Density(1,6,[],[],figuretype) %alldata, ccplot, stormavs-dst, F107MD27d
Density(24,6,[],[],figuretype) %RhoBinned case24, mass-gt20
Density(10,6,[],[],figuretype) %dst-50-tak, and significance table
Density(5,6,[],[],figuretype) %dd12
Density(13,6,[],[],figuretype) %dst-day

%Dissertation plots
Density(28,6,[datenum('Mar-10-1989') datenum('Mar-18-1989')],[],figuretype)
Density(28,6,[datenum('Oct-01-1985') datenum('Oct-11-1985')],[],figuretype)
Density(28,6,[datenum('Oct-22-1988') datenum('Nov-03-1988')],[],figuretype)
Density(28,6,[datenum('Jun-05-1989') datenum('Jun-26-1989')],[],figuretype)