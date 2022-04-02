#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/.Modules";
use Parameter::BinList;
use BamRelated::BaseAndQual;
use BamRelated::Align;
use VcfRelated::VcfParse;
use SeqRelated::Seq;
use Stats::Base;

my ($HelpFlag,$BinList,$BeginTime);
my $ThisScriptName = basename $0;
my ($List4Var,$Bam,$Dir4Log,$File4ConfirmedFreq,$DebugFlag);
my ($Reference,$Samtools,$Bedtools);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script was used to confirm the freq of variants detected by other tools.
  本脚本用于逐个确定变异检测结果中的每个点，在bam中的各种支持数，并最终统计频率。
  注意：read1和read2假如有重叠的话会检查是否一致，假如不一致不会计入AltDepth但会计入FullDepth。

 -v      ( Required ) A list which records the variants\' info (one snv/indel a line, only need \'chr,from,to,ref,alt\');
                      Ref和Alt可以一样，这样就是确认Ref的深度占总深度的比例；
 -b      ( Required ) Bam file;
 -o      ( Required ) File for the confirmed freq info;

 -d      ( Optional ) Debug flag for the printing of related reads info;
 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'v=s' => \$List4Var,
	'b=s' => \$Bam,
	'o=s' => \$File4ConfirmedFreq,
	'd!' => \$DebugFlag,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$List4Var || !$Bam)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	
	IfFileExist($List4Var);
	IfFileExist($Bam);
	die "[ Error ] Not bam file ($Bam).\n" unless($Bam =~ /\.bam$/);
	$Dir4Log = dirname $File4ConfirmedFreq;
	IfDirExist($Dir4Log);
	
	$BinList = BinListGet() if(!$BinList);
	$Reference = BinSearch("Reference",$BinList);
	$Samtools = BinSearch("Samtools",$BinList);
	$Bedtools = BinSearch("Bedtools",$BinList);
}

if(1)
{
	open(VAR,"cat $List4Var | grep -v ^# | cut -f 1-5 |") or die $! unless($List4Var =~ /\.gz$/);
	open(VAR,"zcat $List4Var | grep -v ^# | cut -f 1-5 |") or die $! if($List4Var =~ /\.gz$/);
	open(LOG,"> $File4ConfirmedFreq") or die $! unless($File4ConfirmedFreq =~ /\.gz$/);
	open(LOG,"| gzip > $File4ConfirmedFreq") or die $! if($File4ConfirmedFreq =~ /\.gz$/);
	while(my $VarLine = <VAR>)
	{
		chomp $VarLine;
		my ($Chr,$From,$To,$Ref,$Alt) = split /\t/, $VarLine;
		if($Ref ne $Alt)
		{
			# 以防只找ref;
			($Chr,$From,$To,$Ref,$Alt) = VarSimplify($Chr,$From,$To,$Ref,$Alt);
		}
		$Ref = "" if($Ref eq "-");
		$Alt = "" if($Alt eq "-");
		if($Ref)
		{
			$To = $From + length($Ref) - 1;
			die "[ Error ] Var ($Chr,$From,$To,$Ref) not correct ($VarLine).\n" unless(RefConfirm($Reference,$Chr,$From,$To,$Ref,$Bedtools));
		}
		else
		{
			# 插入;
			$To = $From;
			$Ref = RefGet($Reference,$Chr,$From,$To,$Bedtools);
			$Alt = $Ref . $Alt;
		}
		my $ExactSeq = $Alt;
		$ExactSeq = "-" unless($ExactSeq);
		
		my %ReadSupportInfo = ();
		my %ReadFlagInfo = ();
		open(ORI,"$Samtools view -F 0xD04 $Bam $Chr\:$From\-$To |") or die $!;
		while(my $BamLine = <ORI>)
		{
			chomp $BamLine;
			# 确认是否覆盖;
			# 在计算insert size和variant是否支持时考虑sclip;
			my $SclipFlag = 1;
			my ($RangeChr,$RangeFrom,$RangeTo) = ReadCoverRange($BamLine,$SclipFlag);
			next unless($RangeFrom <= $From && $RangeTo >= $To);
			
			# exact matching check;
			my $RealSeq = AltGet($Chr,$From,$To,$BamLine,$SclipFlag);
			$RealSeq = "-" unless($RealSeq);
			
			# 记录；
			my @Cols = split /\t/, $BamLine;
			$ReadSupportInfo{$Cols[0]} = "" unless($ReadSupportInfo{$Cols[0]});
			if($ReadSupportInfo{$Cols[0]})
			{
				$ReadSupportInfo{$Cols[0]} .= "," . $RealSeq;
				$ReadFlagInfo{$Cols[0]} .= "," . $Cols[1];
			}
			else
			{
				$ReadSupportInfo{$Cols[0]} = $RealSeq;
				$ReadFlagInfo{$Cols[0]} = $Cols[1];
			}
			
			if($DebugFlag)
			{
				print STDERR join("\t","[ Reads ]",$Chr,$From,$To,$Ref,$Alt,$RealSeq,$BamLine),"\n";
			}
		}
		close ORI;
		
		# 统计;
		my ($FullDp,$AltDp,$Freq) = (0,0,"-");
		my %FRNum4Ref = ("F" => 0,"R" => 0,"F1R2" => 0,"F2R1" => 0);
		my %FRNum4Alt = ("F" => 0,"R" => 0,"F1R2" => 0,"F2R1" => 0);
		foreach my $Read (keys %ReadSupportInfo)
		{
			if($ReadSupportInfo{$Read} =~ /,/)
			{
				my @Base = split /,/, $ReadSupportInfo{$Read};
				my @Flag = split /,/, $ReadFlagInfo{$Read};
				die "[ Error ] More than 2 alt items for one read name ($Read).\n" if($#Base > 1);
				if($Base[0] eq $Base[1])
				{
					$AltDp ++ if($Base[0] eq $ExactSeq);
					$FullDp += 1 if($Base[0] ne "N");
					
					my ($FString1,$FRString1) = FRParseFromFlag($Flag[0]);
					$FRNum4Ref{$FRString1} ++ if($Base[0] eq $Ref);
					$FRNum4Alt{$FRString1} ++ if($Base[0] eq $ExactSeq);
				}
				else
				{
					print STDERR "[ Info ] Read1 and Read2 of $Read support two different alt: $Base[0],$Base[1]\n";
				}
			}
			else
			{
				$FullDp ++ if($ReadSupportInfo{$Read} ne "N");
				$AltDp ++ if($ReadSupportInfo{$Read} eq $ExactSeq);
				
				my ($FString1,$FRString1) = FRParseFromFlag($ReadFlagInfo{$Read});
				$FRNum4Ref{$FString1} ++ if($ReadSupportInfo{$Read} eq $Ref);
				$FRNum4Alt{$FString1} ++ if($ReadSupportInfo{$Read} eq $ExactSeq);
				$FRNum4Ref{$FRString1} ++ if($ReadSupportInfo{$Read} eq $Ref);
				$FRNum4Alt{$FRString1} ++ if($ReadSupportInfo{$Read} eq $ExactSeq);
			}
		}
		$Freq = sprintf("%.4f",$AltDp / $FullDp) if($FullDp > 0);
		my @Cols = split /\t/, $VarLine;
		my ($StrandBiasOfF,$StrandBiasOfFR) = ();
		if($Freq <= 0.5)
		{
			$StrandBiasOfF = StrandBiasOfGATKSB($FRNum4Ref{"F"},$FRNum4Alt{"F"},$FRNum4Ref{"R"},$FRNum4Alt{"R"});
			$StrandBiasOfFR = StrandBiasOfGATKSB($FRNum4Ref{"F1R2"},$FRNum4Alt{"F1R2"},$FRNum4Ref{"F2R1"},$FRNum4Alt{"F2R1"});
		}
		else
		{
			$StrandBiasOfF = StrandBiasOfGATKSB($FRNum4Alt{"F"},$FRNum4Ref{"F"},$FRNum4Alt{"R"},$FRNum4Ref{"R"});
			$StrandBiasOfFR = StrandBiasOfGATKSB($FRNum4Alt{"F1R2"},$FRNum4Ref{"F1R2"},$FRNum4Alt{"F2R1"},$FRNum4Ref{"F2R1"});
		}
		$StrandBiasOfF = 1 if($StrandBiasOfF eq "-");
		$StrandBiasOfF = sprintf("%.4f",$StrandBiasOfF);
		$StrandBiasOfFR = 1 if($StrandBiasOfFR eq "-");
		$StrandBiasOfFR = sprintf("%.4f",$StrandBiasOfFR);
		print LOG join("\t",@Cols[0 .. 4],$Freq,$FullDp,$AltDp,$FRNum4Ref{"F"},$FRNum4Ref{"R"},$FRNum4Alt{"F"},$FRNum4Alt{"R"},$StrandBiasOfF,$FRNum4Ref{"F1R2"},$FRNum4Ref{"F2R1"},$FRNum4Alt{"F1R2"},$FRNum4Alt{"F2R1"},$StrandBiasOfFR),"\n";
	}
	close VAR;
	close LOG;
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
