# Package Name
package BamRelated::Align;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(BasicDepthInfoOnChip ReadsNumCountOnChip ReadsNumCountOnRef ReadsPairNumCount PeakInsertSize FRParseFromFlag);

# 确定区域内的平均深度及Fold-80（-F 0xD04），返回平均深度、Fold-80以及深度统计数组;
sub BasicDepthInfoOnChip
{
	my ($Bam,$Bed,$Samtools) = @_;
	
	die "[ Error ] Bam ($Bam) not exist.\n" unless(-s $Bam);
	die "[ Error ] Bed ($Bed) not exist.\n" unless(-s $Bed);
	my ($TotalBase,$BedSize,$AvgDepth,$CutDepth,$Fold80) = (0,0);
	my @DepthInfo = ();
	open(BAM,"$Samtools depth -aa -d 0 -q 0 -Q 0 -b $Bed $Bam | cut -f 3 |") or die $!;
	while(my $Line = <BAM>)
	{
		chomp $Line;
		$DepthInfo[$Line] = 0 unless($DepthInfo[$Line]);
		$DepthInfo[$Line] ++;
		$TotalBase += $Line;
		$BedSize ++;
	}
	close BAM;
	$AvgDepth = $TotalBase / $BedSize;
	my $AccumCutOff = $BedSize * 0.2;
	my $AccumSize = 0;
	for my $i (0 .. $#DepthInfo)
	{
		$DepthInfo[$i] = 0 unless($DepthInfo[$i]);
		$AccumSize += $DepthInfo[$i];
		if($AccumSize >= $AccumCutOff)
		{
			$CutDepth = $i;
			last;
		}
	}
	$Fold80 = $AvgDepth / $CutDepth if($CutDepth);
	$Fold80 = 0 unless($CutDepth);
	
	return $AvgDepth,$Fold80,\@DepthInfo;
}

# 确定区域内去重后的reads数量，包括read1和read2（-F 0xD04）;
sub ReadsNumCountOnChip
{
	my ($Bam,$Bed,$Samtools) = @_;
	
	die "[ Error ] Bam ($Bam) not exist.\n" unless(-s $Bam);
	die "[ Error ] Bed ($Bed) not exist.\n" unless(-s $Bed);
	my $ReadNum = 0;
	$ReadNum = `$Samtools view -c -F 0xD04 -L $Bed $Bam`;
	chomp $ReadNum;
	
	return $ReadNum;
}
sub ReadsNumCountOnRef
{
	my ($Bam,$Samtools) = @_;
	
	die "[ Error ] Bam ($Bam) not exist.\n" unless(-s $Bam);
	my $ReadNum = 0;
	$ReadNum = `$Samtools view -c -F 0xD04 $Bam`;
	chomp $ReadNum;
	
	return $ReadNum;
}

sub ReadsPairNumCount
{
	my ($Bam,$Samtools) = @_;
	
	die "[ Error ] Bam ($Bam) not exist.\n" unless(-s $Bam);
	my $ReadsPairNum = 0;
	$ReadsPairNum = `$Samtools view -c -F 0xD00 -f 0x40 $Bam`;
	chomp $ReadsPairNum;
	
	return $ReadsPairNum;
}

sub PeakInsertSize
{
	my ($Bam,$Samtools) = @_;
	
	die "[ Error ] Bam ($Bam) not exist.\n" unless(-s $Bam);
	my @SizeInfo = 0;
	open(PI,"$Samtools view -F 0xD0C -f 0x40 $Bam | cut -f 9 |") or die $!;
	while(my $Size = <PI>)
	{
		chomp $Size;
		next if($Size == 0);
		
		$Size = abs($Size);
		next if($Size > 100000);
		$SizeInfo[$Size] = 0 unless($SizeInfo[$Size]);
		$SizeInfo[$Size] ++;
	}
	close PI;
	
	my $PeakSize = 0;
	for my $i (1 .. $#SizeInfo)
	{
		$PeakSize = $i if($SizeInfo[$i] > $SizeInfo[$PeakSize]);
	}
	
	return $PeakSize,\@SizeInfo;
}

sub FRParseFromFlag
{
	my $Flag = $_[0];
	
	# Flags:
  # 0x1     PAIRED        .. paired-end (or multiple-segment) sequencing technology
  # 0x2     PROPER_PAIR   .. each segment properly aligned according to the aligner
  # 0x4     UNMAP         .. segment unmapped
  # 0x8     MUNMAP        .. next segment in the template unmapped
  # 0x10    REVERSE       .. SEQ is reverse complemented
  # 0x20    MREVERSE      .. SEQ of the next segment in the template is reversed
  # 0x40    READ1         .. the first segment in the template
  # 0x80    READ2         .. the last segment in the template
  # 0x100   SECONDARY     .. secondary alignment
  # 0x200   QCFAIL        .. not passing quality controls
  # 0x400   DUP           .. PCR or optical duplicate
  # 0x800   SUPPLEMENTARY .. supplementary alignment

	my $ReadId = 0;
	if(($Flag / 64) % 2 == 1)
	{
		$ReadId = 1;
	}
	elsif(($Flag / 128) % 2 == 1)
	{
		$ReadId = 2;
	}
	die "[ Error ] Not read1 or read2 in $Flag\n" unless($ReadId);
	my $PReadId = 1;
	$PReadId = 2 if($ReadId == 1);
	
	my $Ori = 0;
	if(($Flag / 16) % 2 == 1)
	{
		$Ori = 1;
	}
	my $POri = 0;
	if(($Flag / 32) % 2 == 1)
	{
		$POri = 1;
	}
	
	my ($F,$PF) = ("F","F");
	$F = "R" if($Ori == 1);
	$PF = "R" if($POri == 1);
	
	my $FR = "";
	if($F eq "F")
	{
		$FR = $F . $ReadId . $PF . $PReadId;
	}
	elsif($PF eq "F")
	{
		$FR = $PF . $PReadId . $F . $ReadId;
	}
	
	return ($F,$FR);
}

1;
