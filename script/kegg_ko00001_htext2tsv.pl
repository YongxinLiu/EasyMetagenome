#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=1 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}
#my (@tmp1,@tmp2); #database in array
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;
#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";

my %count;
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
#my $a,$b,$c,$d;
while (<INPUT>) {
	chomp;
	if (/^A/) {
		$a=$_;
		$a =~ s/ /\t/;
#		print $a,"\n";
	}	
	if (/^B  /) {
		$b=$_;
		$b =~ s/B  //;
		$b =~ s/ /\t/;
#		print $b,"\n";
	}	
	if (/^C    /) {
		$c=$_;
		$c =~ s/C    //;
		$c =~ s/ /\t/;
#		print $c,"\n";
	}	
	if (/^D      /) {
		$d=$_;
		$d =~ s/D      //;
		$d =~ s/  /\t/;
#		print $d,"\n";
		print OUTPUT "$a\t$b\t$c\t$d\n";
	}
}
close INPUT;
close OUTPUT;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

###############################################################################
#Scripts usage and about.
###############################################################################
sub usage {
    die(
        qq!
Usage:    kegg_ko00001_htext2tsv.pl -i inpute_file -o output_file -d database -h header num
Function: Form KEGG annotation htext format to table format
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, liuyongxin_bio\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/6/2
Notes:    
\n!
    )
}