#!/usr/bin/perl
#Parses tables to combine into one table for the paper


open(FH1,"<tables/CCtable-GOES2.txt");
my @FH1s=<FH1>;
open(FH2,"<tables/CCtable-GOES5.txt");
my @FH2s=<FH2>;
open(FH3,"<tables/CCtable-GOES6.txt");
my @FH3s=<FH3>;
open(FH4,"<tables/CCtable-GOES7.txt");
my @FH4s=<FH4>;
close(FH1); 
close(FH2);
close(FH3);
close(FH4);

my %H1=();
my %H2=();
my %H3=();
my %H4=();

for $line(@FH1s){
 if($line=~m/\d+/){
 $line=~m/([\w\+]+)\s.*([\+-]\d\.\d\d)[^\d]*/;
 $H1{"$1"}=$2;
 }
}

for $line(@FH2s){
 if($line=~m/\d+/){
 $line=~m/([\w\+]+)\s.*([\+-]\d\.\d\d)[^\d]*/;
 $H2{"$1"}=$2;
 }
}
for $line(@FH3s){
 if($line=~m/\d+/){
 $line=~m/([\w\+]+)\s.*([\+-]\d\.\d\d)[^\d]*/;
 $H3{"$1"}=$2;
 }
}
for $line(@FH4s){
 if($line=~m/\d+/){
 $line=~m/([\w\+]+)\s.*([\+-]\d\.\d\d)[^\d]*/;
 $H4{"$1"}=$2;
 }
}


my @Want=('doy','MLT','Bz','Vsw','Dst','Rhosw','F107','Bz+Vsw','Dst+F107','All');

#Print
open(FH,">tables/perltable.tex");
print FH <<EOF;
\\begin{table}[h]
\\small
\\begin{tabular}{|l|llll|}
\\hline
EOF
print FH " & GOES 2 & GOES 5 & GOES 6 & GOES 7\\\\ \\hline\n";
for $wanted(@Want){
print FH "$wanted & $H1{$wanted} & $H2{$wanted} & $H3{$wanted} & $H4{$wanted} \\\\\n";
}
print FH <<EOF;
\\hline
\\end{tabular}
\\caption{Table of linear model correlations showing the median of 100 random samples. Each sample trained on half of the data (via randomly selected rows of the least squares matrix) and tested on the other half} 
\\label{perltable}
\\end{table}
EOF

close(FH);

