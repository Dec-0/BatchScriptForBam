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
my ($OriBam,$FltBam,$IdList,$SeqCheckFlag,$RCFlag,$Samtools);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script was used to extract info from bam with reads\' ID.
  和上个版本相比，新增确认碱基序列功能。

 -i      ( Required ) Bam file;
 -o      ( Required ) Filterd bam file;
 -r      ( Required ) A file records the reads\' ID, one ID a line;
 -b      ( Required ) If checking the base seq which will need the sequence;
                      需要提供碱基序列。
 -c      ( Required ) If the seq need to be reverse conplement;

 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$OriBam,
	'o=s' => \$FltBam,
	'r=s' => \$IdList,
	'b!' => \$SeqCheckFlag,
	'c!' => \$RCFlag,
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
	my $CutString = "1";
	$CutString = "1-" if($SeqCheckFlag);
	open(ID,"cat $IdList | cut -f $CutString |") or die $! unless($IdList =~ /\.gz$/);
	open(ID,"zcat $IdList | cut -f $CutString |") or die $! if($IdList =~ /\.gz$/);
	while(my $Line = <ID>)
	{
		chomp $Line;
		
		if(!$tmpHash{$Line})
		{
			$tmpHash{$Line} = 1;
			$TotalNum ++;
		}
		
		# 是否需要反向互补;
		if($RCFlag)
		{
			my @Cols = split /\t/, $Line;
			$Cols[1] = uc($Cols[1]);
			$Cols[1] = reverse($Cols[1]);
			$Cols[1] =~ tr/ATCG/TAGC/;
			$Line = join("\t",@Cols);
			$tmpHash{$Line} = 1 if(!$tmpHash{$Line});
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
		if($SeqCheckFlag)
		{
			next unless($Cols[9]);
			my $tKey = join("\t",$Cols[0],$Cols[9]);
			print FLT $Line if($tmpHash{$tKey});
		}
		else
		{
			print FLT $Line if($tmpHash{$Cols[0]});
		}
	}
	close ORI;
	close FLT;
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
