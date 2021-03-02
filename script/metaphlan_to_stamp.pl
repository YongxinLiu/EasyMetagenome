#!/usr/bin/perl
#Created by: Morgan Langille

#Note: metaphlan file must contain merged files
#Note2: 'all' taxonomy must be output in metaphlan

#metaphlan_to_stamp.pl merged_metaphlan_output.txt > stamp_profile.spf

use warnings;
use strict;

my $header=<>;

chomp($header);

my @cols=split(/\t/,$header);

#Just keep the sample names and get rid of the "ID"
shift(@cols);

my @taxa_ranks=("Kingdom","Phylum","Class","Order","Family","Genus","Species");

#start with assumption that this is not metaphlan2
my $metaphlan2_version=0;

#flag to say if we have printed the header yet.
my $printed_header=0;

LINE:
while(<>){
    chomp;
    if(/^#/){
	#if second line begins with a '#' then we know the file is metaphlan2
	$metaphlan2_version=1;
	#metaphlan2 goes down to the strain level so need to add an extra rank here
	push(@taxa_ranks,'Strain');
	next;
    }
    unless($printed_header){
	print join("\t",@taxa_ranks,@cols),"\n";
	
	#set flag so we don't keep printing the header
	$printed_header=1;
    }

    my @data=split;
    my @taxonomy= split(/\|/,shift(@data));
    my $unclassified_flag=0;
    for my $x (0..$#taxa_ranks){
	if(defined($taxonomy[$x])){
	    if($taxonomy[$x] =~ /unclassified/){
		$taxonomy[$x]='unclassified';
		$unclassified_flag=1;
	    }
	}elsif($unclassified_flag ==1){
	    $taxonomy[$x]='unclassified';
	}else{
	    #not defined and not an unclassified parent so do not include in output
	    next LINE;
	}
    }
    print join("\t",@taxonomy,@data),"\n";
}

