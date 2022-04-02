#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/.Modules";
use Parameter::BinList;

my ($HelpFlag,$BinList,$BeginTime);
my $ThisScriptName = basename $0;
my ($OriBam,$FltBam,$IdList,$Samtools);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script was used to extract info from bam with reads\' ID.

 -i      ( Required ) Bam file;
 -o      ( Required ) Filterd bam file;
 -r      ( Required ) A file records the reads\' ID, one ID a line;

 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$OriBam,
	'o=s' => \$FltBam,
	'r=s' => \$IdList,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$OriBam || !$FltBam || !$IdList)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	IfFileExist($OriBam);
	IfFileExist($IdList);
	my $Dir = dirname $FltBam;
	IfDirExist($Dir);
	
	$BinList = BinListGet() if(!$BinList);
	$Samtools = BinSearch("Samtools",$BinList);
}

if(1)
{
	my %tmpHash = ();
	my $TotalNum = 0;
	open(ID,"cat $IdList |") or die $! unless($IdList =~ /\.gz$/);
	open(ID,"zcat $IdList |") or die $! if($IdList =~ /\.gz$/);
	while(my $Line = <ID>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		if(!$tmpHash{$Cols[0]})
		{
			$tmpHash{$Cols[0]} = 1;
			$TotalNum ++;
		}
	}
	close ID;
	print "[ Info ] Total reads id: $TotalNum\n";
	
	open(ORI,"$Samtools view -h $OriBam |") or die $!;
	open(FLT,"| $Samtools view -bh > $FltBam") or die $!;
	while(my $Line = <ORI>)
	{
		if($Line =~ /^@/)
		{
			print FLT $Line;
			next;
		}
		
		my @Cols = split /\t/, $Line;
		print FLT $Line if($tmpHash{$Cols[0]});
	}
	close ORI;
	close FLT;
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
