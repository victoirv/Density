#!/usr/bin/perl
#Parses tables to combine into one table for the paper
#use warnings;
#use strict;

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


 my @AoF;
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

 my $findex=0;
 @AoF[$findex]=[@FH1s];$findex++;
 @AoF[$findex]=[@FH2s];$findex++;
 @AoF[$findex]=[@FH3s];$findex++;
 @AoF[$findex]=[@FH4s];$findex++;

 my %H1=();

 my $headerrx = 'Vars\s+(.*)';
 my @headers = ();
 my $regex = '^(\S+)\s+- (.\d.*)$'; #Catch everything, split later
 my @Sats=("GOES 2","GOES 5","GOES 6","GOES 7");
 my $satindex=0;

 for my $FHs(0 .. $#Sats){
  for my $line(@{$AoF[$FHs]}){
   if($line=~m/$headerrx/)
   {
    @headers=split(/\s+/,$1);
   }
   if($line=~m/\d\d+/)
   {
    $line=~/$regex/;
    my @vals=split(" ",$2);
    for(my $i=0; $i<scalar(@headers); $i++)
    {
     my $AddString="$Sats[$satindex]-$1-@headers[$i]";
     $H1{$AddString}=@vals[$i];
    }
   }
  }
  $satindex=$satindex+1;
 }

 my @Want=('DoY','MLT','B_z','V_{sw}','D_{st}','\rho_{sw}','F_{10.7}','B_z+V_{sw}','D_{st}+F_{10.7}','All');
 my @WantHs = ();

 if($Type=~/NN/){
  @WantHs=('NNv','+-NNv');
 }
 else
 {
  @WantHs=('CCt','+-CCt');
 }

#Print
 open(FH,">${Type}perltable.tex");
 print FH <<EOF;
 \\begin{table}[h]
 \\small
 \\begin{tabular}{|L|LLLL|}
 \\hline
EOF

 print FH ' & \text{GOES 2} & \text{GOES 5} & \text{GOES 6} & \text{GOES 7}\\\\ \hline'."\n";
 for my $wanted(@Want){
  print FH "$wanted ";
  for my $sat(@Sats){
   print FH "& ";
   for my $wanth(@WantHs){
    my $HashString="$sat-$wanted-$wanth";
    if($wanth=~/\+/)
    {
     print FH "\\pm$H1{$HashString} "
    }
    else{
     print FH "$H1{$HashString}"
    }
   }
  }
 print FH "\\\\\n";
 }

 print FH <<EOF;
 \\hline
 \\end{tabular}
 \\caption{Table of $TypeString model test correlations showing the median of 100 random samples. Each sample trained on half of the data (via randomly selected rows of the least squares matrix) and tested on the other half} 
 \\label{${Type}perltable}
 \\end{table}
EOF

  close(FH);


 #Print difference table
 if($Type=~/NN/){
  @WantHs=(['NNt','NNtr'],['+-NNt','+-NNtr']);
 }
 else
 {
  @WantHs=(['CCt','CCtr'],['+-CCt','+-CCtr']);
 }

 open(FH,">${Type}difftable.tex");
 print FH <<EOF;
 \\begin{table}[h]
 \\small
 \\begin{tabular}{|L|LLLL|}
 \\hline
EOF

 print FH 'T-Tr & \text{GOES 2} & \text{GOES 5} & \text{GOES 6} & \text{GOES 7}\\\\ \hline'."\n";
 for my $wanted(@Want){
  print FH "$wanted ";
  for my $sat(@Sats){
   print FH "& ";
   for my $wanth(@WantHs){
    
    if(@$wanth[0]=~/\+/)
    {
    my $DiffVal=sprintf("%2.2f",sqrt($H1{"$sat-$wanted-@$wanth[0]"}**2+$H1{"$sat-$wanted-@$wanth[1]"}**2));
     print FH "\\pm$DiffVal "
    }
    else{
    my $DiffVal=sprintf("%2.2f",$H1{"$sat-$wanted-@$wanth[0]"}-$H1{"$sat-$wanted-@$wanth[1]"});
     print FH "$DiffVal"
    }
   }
  }
 print FH "\\\\\n";
 }

 print FH <<EOF;
 \\hline
 \\end{tabular}
 \\caption{Table of differences in $TypeString testing-training models, where each correlation is the median correlation of 100 random samples. Each sample trained on half of the data (via randomly selected rows of the least squares matrix) and tested on the other half} 
 \\label{${Type}difftable}
 \\end{table}
EOF

  close(FH);



}
