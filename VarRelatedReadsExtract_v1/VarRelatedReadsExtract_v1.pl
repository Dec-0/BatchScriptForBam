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
use VcfRelated::VcfParse;
use SeqRelated::Seq;
use Sort::ChrPos;

my ($HelpFlag,$BinList,$BeginTime);
my $ThisScriptName = basename $0;
my ($VarList,$Bam,$LogBam,$LogDir);
my ($Reference,$Samtools,$Bedtools);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script was used to extract reads supporting specific variants like snv or indel from bam.

 -v      ( Required ) A list which records the variants\' info (one snv/indel a line, include \'chr,from,to,ref,alt\');
 -b      ( Required ) Bam file;
 -o      ( Required ) The bam file for the filtered reads;

 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'v=s' => \$VarList,
	'b=s' => \$Bam,
	'o=s' => \$LogBam,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$VarList || !$Bam || !$LogBam)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	
	IfFileExist($VarList);
	IfFileExist($Bam);
	die "[ Error ] Not bam file ($Bam).\n" unless($Bam =~ /\.bam$/);
	die "[ Error ] Not bam file ($LogBam).\n" unless($LogBam =~ /\.bam$/);
	$LogDir = dirname $LogBam;
	IfDirExist($LogDir);
	
	$BinList = BinListGet() if(!$BinList);
	$Reference = BinSearch("Reference",$BinList);
	$Samtools = BinSearch("Samtools",$BinList);
	$Bedtools = BinSearch("Bedtools",$BinList);
}

if(1)
{
	# 创建一个临时bed，并过滤出相关的reads信息;
	my $FltBam = $LogBam;
	$FltBam =~ s/bam$/flt.bam/;
	my $FltBed = $LogBam;
	$FltBed =~ s/bam$/flt.bed/;
	# 输入的bed可能有重叠，比如chrX 17151344 17151346和chrX 17151345 17151346，会使bedtools报错“Error: Sorted input specified, but the file - has the following out of order record”;
	#`cat $VarList | grep -v ^# | cut -f 1-3 | $Bedtools sort -i - | $Bedtools merge -i - > $FltBed` unless($VarList =~ /\.gz$/);
	#`zcat $VarList | grep -v ^# | cut -f 1-3 | $Bedtools sort -i - | $Bedtools merge -i - > $FltBed` if($VarList =~ /\.gz$/);
	#`cat $VarList | grep -v ^# | cut -f 1-3 | sort -n -k1.4 -k2 -k3 | $Bedtools merge -i - > $FltBed` unless($VarList =~ /\.gz$/);
	#`zcat $VarList | grep -v ^# | cut -f 1-3 | sort -n -k1.4 -k2 -k3 | $Bedtools merge -i - > $FltBed` if($VarList =~ /\.gz$/);
	if(1)
	{
		my (@Chr,@From,@To,@Other) = ();
		open(VARL,"zcat $VarList | grep -v ^# | cut -f 1-3 |") or die $! if($VarList =~ /\.gz$/);
		open(VARL,"cat $VarList | grep -v ^# | cut -f 1-3 |") or die $! unless($VarList =~ /\.gz$/);
		while(my $Line = <VARL>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			push @Chr, $Cols[0];
			push @From, $Cols[1];
			push @To, $Cols[2];
			push @Other, "";
		}
		close VARL;
		my @tRef = ChrPosAndOther(\@Chr,\@From,\@To,\@Other);
		@Chr = @{$tRef[0]};
		@From = @{$tRef[1]};
		@To = @{$tRef[2]};
		
		my $tmpFile = $FltBed . ".tmp";
		open(TMP,"> $tmpFile") or die $!;
		for my $i (0 .. $#Chr)
		{
			print TMP join("\t",$Chr[$i],$From[$i],$To[$i]),"\n";
		}
		close TMP;
		`cat $tmpFile | $Bedtools merge -i - > $FltBed`;
		`rm $tmpFile`;
	}
	# 除去unmap、secondary、dup、supplementary;
	`$Samtools view -bh -F 0xD04 -L $FltBed $Bam > $FltBam`;
	`$Samtools index $FltBam`;
	
	# 需要先清理掉很多临时bam
	if(1)
	{
		my $tInfo = `ls $LogBam\.tmp.*.bam 2>/dev/null`;
		my @tBam = split /\n/, $tInfo;
		for my $i (0 .. $#tBam)
		{
			next unless(-s $tBam[$i]);
			next unless($tBam[$i] =~ /.tmp.\d+.bam$/);
			`rm $tBam[$i]`;
		}
	}
	open(VAR,"cat $VarList | grep -v ^# |") or die $! unless($VarList =~ /\.gz$/);
	open(VAR,"zcat $VarList | grep -v ^# |") or die $! if($VarList =~ /\.gz$/);
	open(LOG,"| $Samtools view -bh | $Samtools sort - -o $LogBam") or die $!;
	my $VarId = 0;
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
			die "[ Error ] Var not correct ($VarLine).\n" unless(RefConfirm($Reference,$Chr,$From,$To,$Ref,$Bedtools));
			
			$To = $From + length($Ref) - 1;
		}
		else
		{
			# 插入;
			$To = $From;
			$Ref = RefGet($Reference,$Chr,$From,$To,$Bedtools);
			$Alt = $Ref . $Alt;
		}
		my $ExactSeq = $Alt;
		
		unless($VarId)
		{
			# 记录bam头;
			my $Head = `$Samtools view -H $FltBam`;
			print LOG $Head;
			$VarId ++;
		}
		
		open(ORI,"$Samtools view $FltBam $Chr\:$From\-$To |") or die $!;
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
			next unless($RealSeq eq $ExactSeq);
			
			print LOG $BamLine,"\n";
		}
		close ORI;
	}
	close VAR;
	close LOG;
	
	`rm $FltBed $FltBam $FltBam\.bai`;
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
