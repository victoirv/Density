#!/usr/bin/perl
#Parses tables to combine into one table for the paper
use warnings;
use strict;

my @Types=('CC','NN');

if($#ARGV+1 > 0) {
 @Types=($ARGV[0]);
 if($ARGV[0]!~m/(NN|CC)/){
  print "Type must be NN or CC\n";
  exit
  }
}

foreach my $Type(@Types){

my $TypeString="linear";
if($Type=~m/NN/)
{$TypeString="nonlinear";}

open(FH1,"<${Type}table-GOES2.txt");
my @FH1s=<FH1>;
open(FH2,"<${Type}table-GOES5.txt");
my @FH2s=<FH2>;
open(FH3,"<${Type}table-GOES6.txt");
my @FH3s=<FH3>;
open(FH4,"<${Type}table-GOES7.txt");
my @FH4s=<FH4>;
close(FH1); 
close(FH2);
close(FH3);
close(FH4);

my %H1=();
my %H2=();
my %H3=();
my %H4=();

my $regex = '^(\S+)\s+- .\d\.\d\d (.\d\.\d\d) \d.\d\d (\d.\d\d).*'; #To catch CC test column and test_sd
if($Type=~m/NN/) {
$regex = '^(\S+)\s+.*(?:\s*.\d\.\d\d){5} (.\d.\d\d) (?:\s*\d.\d\d){5} (\d.\d\d)[^\d]*$'; #To catch NN validation column
}

for my $line(@FH1s){
 if($line=~m/\d+/){
 $line=~/$regex/;
 $H1{"$1"}=$2;
 $H1{"$1-sd"}=$3;
 }
}

for my $line(@FH2s){
 if($line=~m/\d+/){
 $line=~/$regex/;
 $H2{"$1"}=$2;
 $H2{"$1-sd"}=$3;
 }
}
for my $line(@FH3s){
 if($line=~m/\d+/){
 $line=~/$regex/;
 $H3{"$1"}=$2;
 $H3{"$1-sd"}=$3;
 }
}
for my $line(@FH4s){
 if($line=~m/\d+/){
 $line=~/$regex/;
 $H4{"$1"}=$2;
 $H4{"$1-sd"}=$3;
 }
}


my @Want=('DoY','MLT','B_z','V_{sw}','D_{st}','\rho_{sw}','F_{10.7}','B_z+V_{sw}','D_{st}+F_{10.7}','All');

#Print
open(FH,">${Type}perltable.tex");
print FH <<EOF;
\\begin{table}[h]
\\small
\\begin{tabular}{|L|LLLL|}
\\hline
EOF
print FH 'Test & \text{GOES 2} & \text{GOES 5} & \text{GOES 6} & \text{GOES 7}\\\\ \hline'."\n";
for my $wanted(@Want){
print FH "$wanted & $H1{$wanted}\\pm$H1{${wanted}.'-sd'} & $H2{$wanted}\\pm$H2{${wanted}.'-sd'} & $H3{$wanted}\\pm$H3{${wanted}.'-sd'} & $H4{$wanted}\\pm$H4{${wanted}.'-sd'} \\\\\n";
}
print FH <<EOF;
\\hline
\\end{tabular}
\\caption{Table of $TypeString model test correlations showing the median of 100 random samples. Each sample trained on half of the data (via randomly selected rows of the least squares matrix) and tested on the other half} 
\\label{${Type}perltable}
\\end{table}
EOF

close(FH);

}
