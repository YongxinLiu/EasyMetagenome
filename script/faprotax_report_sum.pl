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
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my %database;
my %database2;
$function="";
$OTU="";
while (<INPUT>) {
	chomp;
	if (/^#/) {
		@tmp=split/\s+/;
		$function=$tmp[1];
		#print $function,"\n";
	}else{
		@tmp=split/\s+/;
		$OTU=$tmp[1];
		#print $OTU,"\n";
		$database{$OTU}{$function}=1;	
		$database2{$function}{$OTU}=1;	
	}
	
}
close INPUT;

# 保存OTU-function列表
open OUTPUT,">$opts{o}.otu_func";
foreach  my $key1 (sort keys %database) {
	foreach my $key2 (sort keys %{$database{$key1}}) { # sort {$a<=>$b;}
		print OUTPUT "$key1\t$key2\n" if defined($database{$key1}{$key2}); # \t$key2\t$database{$key1}{$key2}
	}
}
close OUTPUT;
print "OTU功能注释列表: $opts{o}.otu_func\n";

# 保存功能-OTU列表
open OUTPUT,">$opts{o}.func_otu";
foreach  my $key1 (sort keys %database2) {
	foreach my $key2 (sort keys %{$database2{$key1}}) { 
		print OUTPUT "$key1\t$key2\n" if defined($database2{$key1}{$key2});
	}
}
close OUTPUT;
print "功能包含OTU列表: $opts{o}.func_otu\n";

# 保存功能-OTU列表
open OUTPUT,">$opts{o}.mat";
# 写功能表头
print OUTPUT "OTUID";
foreach  my $key1 (sort keys %database2) {
	print OUTPUT "\t$key1";
}
print OUTPUT "\n";
# 写OTU每行
foreach  my $key1 (sort keys %database) {
	print OUTPUT "$key1";
	foreach my $key2 (sort keys %database2) { 
		if (!defined($database{$key1}{$key2})) {
			$database{$key1}{$key2}=0;
			#print "\t$database{$key1}{$key2}";
			print OUTPUT "\t$database{$key1}{$key2}";
		}else{
			print OUTPUT "\t$database{$key1}{$key2}";
		}
	}
	print OUTPUT "\n";
}
close OUTPUT;
print "OTU功能有无矩阵: $opts{o}.mat\n";



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
Usage:    template.pl -i inpute_file -h header num
Function: 计算faprotax报告的信息，获得OTU-功能，功能-OTU，OTU-功能矩阵
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2018/5/24
Notes:    
\n!
    )
}