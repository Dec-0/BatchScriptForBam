# Package Name
package SeqRelated::Seq;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(RefGet RefGetWithOutWarning ExonCoordInfo CdsOfTranscript SclipExonLocate Nucleo2Amino NucleoBaseRevise GCPercentOfReference WhichGeneInArea FqFilterById);
use VcfRelated::VcfParse;
use Parameter::BinList;
use File::Basename;

use FindBin qw($Bin);

sub RefGet
{
	my ($RefGen,$Chr,$From,$To,$BedtoolsBin) = @_;
	
	die "[ Error ] End position can not be smaller than start position when getting seq ($Chr,$From,$To) from ref($RefGen).\n" if($From > $To);
	die "[ Error ] Reference not exist ($RefGen).\n" unless(-e $RefGen);
	$BedtoolsBin = "bedtools" unless($BedtoolsBin);
	
	# 这里默认输入的是1-based的坐标，而bedtools针对的是0-based的坐标，因此假如想取得坐标[a,b]的碱基，则需要转换成(a - 1, b);
	$From --;
	my $Ref = `echo -e '$Chr\\t$From\\t$To' | $BedtoolsBin getfasta -fi $RefGen -bed - | tail -n 1`;
	chomp $Ref;
	$Ref = uc($Ref);
	
	return $Ref;
}

sub RefGetWithOutWarning
{
	my ($RefGen,$Chr,$From,$To) = @_;
	
	my $Ref = "";
	return $Ref if($From > $To);
	die "[ Error ] Reference not exist ($RefGen).\n" unless(-e $RefGen);
	
	$From --;
	$Ref = `echo -e '$Chr\\t$From\\t$To' | bedtools getfasta -fi $RefGen -bed - | tail -n 1`;
	chomp $Ref;
	$Ref = uc($Ref);
	
	return $Ref;
}

sub ExonCoordInfo
{
	my ($refGene,$Gene) = @_;
	my %GeneInfo = ();
	
	die "[ Error ] refGene not exist ($refGene).\n" unless($refGene && -s $refGene);
	die "[ Error ] No gene name ($Gene).\n" unless($Gene);
	
	my $Return = ();
	$Return = `zcat $refGene | awk '{if(\$13 == "$Gene"){print \$0}}'` if($refGene =~ /\.gz$/);
	$Return = `cat $refGene | awk '{if(\$13 == "$Gene"){print \$0}}'` unless($refGene =~ /\.gz$/);
	my @Items = split /\n/, $Return;
	my ($MinCoord,$MaxCoord) = (0,0);
	for my $i (0 .. $#Items)
	{
		my @Cols = split /\t/, $Items[$i];
		die "[ Error ] Not enough columes ($Items[$i]).\n" unless($#Cols >= 15);
		
		$GeneInfo{$Cols[1]}{"chrom"} = $Cols[2];
		$GeneInfo{$Cols[1]}{"strand"} = $Cols[3];
		$GeneInfo{$Cols[1]}{"txStart"} = $Cols[4];
		$GeneInfo{$Cols[1]}{"txEnd"} = $Cols[5];
		$GeneInfo{$Cols[1]}{"cdsStart"} = $Cols[6];
		$GeneInfo{$Cols[1]}{"cdsEnd"} = $Cols[7];
		$GeneInfo{$Cols[1]}{"exonCount"} = $Cols[8];
		$Cols[9] =~ s/,$//;
		$GeneInfo{$Cols[1]}{"exonStarts"} = $Cols[9];
		$Cols[10] =~ s/,$//;
		$GeneInfo{$Cols[1]}{"exonEnds"} = $Cols[10];
		$GeneInfo{$Cols[1]}{"cdsStartStat"} = $Cols[13];
		$GeneInfo{$Cols[1]}{"cdsEndStat"} = $Cols[14];
		$GeneInfo{$Cols[1]}{"exonFrames"} = $Cols[15];
		
		if($i == 0)
		{
			$MinCoord = $Cols[4];
			$MaxCoord = $Cols[5];
		}
		else
		{
			$MinCoord = $Cols[4] if($Cols[4] < $MinCoord);
			$MaxCoord = $Cols[5] if($Cols[5] > $MaxCoord);
		}
	}
	$GeneInfo{"All"}{"MinCoord"} = $MinCoord;
	$GeneInfo{"All"}{"MaxCoord"} = $MaxCoord;
	
	return \%GeneInfo;
}

# 获取转录本的CDs区间及转录本碱基序列;
sub CdsOfTranscript
{
	my ($refGene,$TranscriptId,$RefGen,$BedtoolsBin) = @_;
	my $CdSeq = "";
	
	die "[ Error ] refGene not exist ($refGene).\n" unless($refGene && -s $refGene);
	die "[ Error ] No transcript id specified ($TranscriptId).\n" unless($TranscriptId);
	
	my $Return = ();
	$Return = `zcat $refGene | awk '{if(\$2 == "$TranscriptId"){print \$0}}' | head -n1` if($refGene =~ /\.gz$/);
	$Return = `cat $refGene | awk '{if(\$2 == "$TranscriptId"){print \$0}}' | head -n1` unless($refGene =~ /\.gz$/);
	chomp $Return;
	my @Cols = split /\t/, $Return;
	my ($Chr,$Ori,$CdStart,$CdEnd,$ExonStart,$ExonEnd) = (@Cols[2 .. 3],@Cols[6 .. 7],@Cols[9 .. 10]);
	
	# 获得编码区坐标;
	my (@CdStartArray,@CdEndArray) = ();
	if(1)
	{
		my @ExonStartArray = split /,/, $ExonStart;
		my @ExonEndArray = split /,/, $ExonEnd;
		for my $i (0 .. $#ExonStartArray)
		{
			if($ExonStartArray[$i] <= $CdEnd && $ExonEndArray[$i] >= $CdStart)
			{
				my $tStart = $ExonStartArray[$i];
				my $tEnd = $ExonEndArray[$i];
				$tStart = $CdStart if($CdStart > $tStart);
				$tEnd = $CdEnd if($CdEnd < $tEnd);
				
				# 将bed格式转换成1-based;
				$tStart ++;
				push @CdStartArray, $tStart;
				push @CdEndArray, $tEnd;
			}
		}
	}
	
	# 获取碱基序列;
	for my $i (0 .. $#CdStartArray)
	{
		# hg19_refGene.txt文件中默认是bed格式的区间，而RefGet接受的是1-based的坐标;
		$CdSeq .= &RefGet($RefGen,$Chr,$CdStartArray[$i],$CdEndArray[$i],$BedtoolsBin);
	}
	
	return ($CdSeq,$Chr,$Ori,\@CdStartArray,\@CdEndArray);
}

# 根据基因、断点坐标、融合侧方向，确定融合片段所属的外显子编号、当前断点所属的外显子或内含子编号;
sub SclipExonLocate
{
	my ($refGene,$Gene,$Pos,$Ori) = @_;
	my (@Transcript,@CurrentId,@ExonStartId,@ExonEndId) = ();
	
	# 默认截取5'端序列;
	$Ori = "+" unless($Ori);
	if($Ori !~ /\D/ && $Ori > 0)
	{
		$Ori = "+";
	}
	
	my %GeneInfo = %{&ExonCoordInfo($refGene,$Gene)};
	foreach my $NMId (keys %GeneInfo)
	{
		# 只处理“NM_”开头的转录本;
		next if($NMId eq "All" && $NMId !~ /^NM_/);
		
		# 获取Coding Sequence对应的外显子序号，注意不是整体的外显子编号;
		my (@CdStartArray,@CdEndArray,@CdId) = ();
		my @ExonStartArray = split /,/, $GeneInfo{$NMId}{"exonStarts"};
		my @ExonEndArray = split /,/, $GeneInfo{$NMId}{"exonEnds"};
		for my $i (0 .. $#ExonStartArray)
		{
			if($ExonStartArray[$i] <= $GeneInfo{$NMId}{"cdsEnd"} && $ExonEndArray[$i] >= $GeneInfo{$NMId}{"cdsStart"})
			{
				my $tStart = $ExonStartArray[$i];
				my $tEnd = $ExonEndArray[$i];
				$tStart = $GeneInfo{$NMId}{"cdsStart"} if($GeneInfo{$NMId}{"cdsStart"} > $tStart);
				$tEnd = $GeneInfo{$NMId}{"cdsEnd"} if($GeneInfo{$NMId}{"cdsEnd"} < $tEnd);
				
				# 注意这是bed格式的坐标，是0-based的;
				push @CdStartArray, $tStart;
				push @CdEndArray, $tEnd;
			}
		}
		# 假如都没有覆盖任何外显子，则忽略;
		next unless(@CdStartArray);
		
		# 给CD区域的外显子按编码方向编号;
		my $TotalCdNum = @CdStartArray;
		my $StrandOri = $GeneInfo{$NMId}{"strand"};
		if($StrandOri eq "+")
		{
			for my $i (0 .. $#CdStartArray)
			{
				$CdId[$i] = $i + 1;
			}
		}
		else
		{
			my $tId = 1;
			for(my $i = $#CdStartArray;$i >= 0;$i --)
			{
				$CdId[$i] = $tId;
				$tId ++;
			}
		}
		
		# 确定当前所处的Exon或Intron编号，以及所囊括的CD编号;
		my @Id = ();
		if($Ori eq "+")
		{
			# 保留的是5'端;
			for(my $i = $#CdStartArray;$i >= 0;$i --)
			{
				if($CdStartArray[$i] <= $Pos)
				{
					push @Id, $CdId[$i];
				}
			}
		}
		else
		{
			# 保留的是3'端;
			for my $i (0 .. $#CdStartArray)
			{
				if($CdEndArray[$i] >= $Pos)
				{
					push @Id, $CdId[$i];
				}
			}
		}
		# 假如一个外显子都没有覆盖，则记录0;
		push @Id, "0" unless(@Id);
		
		# 确定断点所处的位置;
		my $cId = "-";
		for my $i (0 .. $#CdStartArray)
		{
			last unless($cId eq "-");
			
			if($CdStartArray[$i] > $Pos)
			{
				# 处于上一个外显子和当前外显子坐标之间;
				if($i == 0)
				{
					$cId = "5'-UTR" if($StrandOri eq "+");
					$cId = "3'-UTR" unless($StrandOri eq "+");
				}
				else
				{
					$cId = "Intron" . $CdId[$i - 1] if($StrandOri eq "+");
					$cId = "Intron" . $CdId[$i] unless($StrandOri eq "+");
				}
			}
			elsif($CdEndArray[$i] > $Pos)
			{
				# 处于当前外显子内;
				$cId = "Exon" . $CdId[$i];
			}
		}
		# 都不是就是另一侧UTR;
		if($cId eq "-")
		{
			$cId = "3'-UTR" if($StrandOri eq "+");
			$cId = "5'-UTR" unless($StrandOri eq "+");
		}
		
		push @Transcript, $NMId;
		push @CurrentId, $cId;
		push @ExonStartId, $Id[0];
		push @ExonEndId, $Id[-1];
	}
	
	return (\@Transcript,\@CurrentId,\@ExonStartId,\@ExonEndId);
}

# 确定属于哪个基因;
sub WhichGeneInArea
{
	my ($refGene,$Chr,$From,$To) = @_;
	$To = $From unless($To);
	
	my $Return = "";
	$Return = `zcat $refGene | awk '{if(\$3 == \"$Chr\" && \$5 <= $To && \$6 >= $From){print \$0}}' | cut -f 13 | sort | uniq` if($refGene =~ /\.gz$/);
	$Return = `cat $refGene | awk '{if(\$3 == \"$Chr\" && \$5 <= $To && \$6 >= $From){print \$0}}' | cut -f 13 | sort | uniq` unless($refGene =~ /\.gz$/);
	chomp $Return;
	my @Gene = split /\n/, $Return;
	$Return = join(",",@Gene);
	
	return $Return;
}

# 将碱基序列翻译成氨基酸序列;
sub Nucleo2Amino
{
	my ($BaseSeq,$InitialId,$Ori) = @_;
	$InitialId = 0 unless($InitialId);
	if($InitialId)
	{
		die "[ Error ] InitialId should only be numbers.\n" if($InitialId =~ /\D/);
		die "[ Error ] InitialId not 1 or 2.\n" if($InitialId == 1 || $InitialId == 2);
	}
	
	# Codon Table;
	my(%AACode) = ('TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'TTC' => 'F', 'TTT' => 'F', 'TTA' => 'L', 'TTG' => 'L', 'TAC' => 'Y',
    'TAT' => 'Y', 'TAA' => '-', 'TAG' => '-', 'TGC' => 'C', 'TGT' => 'C', 'TGA' => '-', 'TGG' => 'W', 'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L',
    'CTT' => 'L', 'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P', 'CAC' => 'H', 'CAT' => 'H', 'CAA' => 'Q', 'CAG' => 'Q', 'CGA' => 'R',
    'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I', 'ATG' => 'M', 'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T',
    'ACT' => 'T', 'AAC' => 'N', 'AAT' => 'N', 'AAA' => 'K', 'AAG' => 'K', 'AGC' => 'S', 'AGT' => 'S', 'AGA' => 'R', 'AGG' => 'R', 'GTA' => 'V',
    'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V', 'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A', 'GAC' => 'D', 'GAT' => 'D', 'GAA' => 'E',
    'GAG' => 'E', 'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G');
	
	# Ori3Flag指转录本负链编码，默认是正链，但是没有指定的话就是负链;
	my $Ori3Flag = 0;
	if(length($Ori) == 0)
	{
		$Ori3Flag = 1;
	}
	elsif($Ori eq "-")
	{
		$Ori3Flag = 1;
	}
	elsif($Ori =~ /^\d+$/)
	{
		$Ori3Flag = 1 if($Ori <= 0);
	}
	if($Ori3Flag)
	{
		$BaseSeq = reverse $BaseSeq;
		$BaseSeq = uc($BaseSeq);
		$BaseSeq =~ tr/ATCG/TAGC/;
	}
	
	my $AASeq = "";
	$BaseSeq = uc($BaseSeq);
	my $MaxLen = length($BaseSeq);
	my $MaxId = $MaxLen - 1;
	my $CurrentId = $InitialId;
	while($CurrentId + 2 <= $MaxId)
	{
		my $Codon = substr($BaseSeq,$CurrentId,3);
		die "[ Error ] No amino acid for codon $Codon .\n" unless($AACode{$Codon});
		
		if($AASeq)
		{
			last if($AACode{$Codon} eq "-");
			$AASeq .= $AACode{$Codon};
		}
		elsif($AACode{$Codon} eq "M")
		{
			$AASeq = $AACode{$Codon};
		}
		
		$CurrentId += 3;
	}
	
	return $AASeq;
}

# 将变异整合进碱基序列中并做好标记;
sub NucleoBaseRevise
{
	die "[ Error ] Not 7 paramters.\n" unless($#_ == 6);
	my $NucString = $_[0];
	my $RefGen = $_[1];
	my $Chr = $_[2];
	my $Ori = $_[3];
	my @Start = @{$_[4]};
	my @End = @{$_[5]};
	# 变异的格式必须是“chr1,1234,A,-”;
	my @Var = @{$_[6]};
	
	# Ori3Flag指转录本负链编码，默认是正链，但是没有指定的话就是负链;
	my $Ori3Flag = 0;
	if(length($Ori) == 0)
	{
		$Ori3Flag = 1;
	}
	elsif($Ori eq "-")
	{
		$Ori3Flag = 1;
	}
	elsif($Ori =~ /^\d+$/)
	{
		$Ori3Flag = 1 if($Ori <= 0);
	}
	
	my @Base = split //, $NucString;
	my $BaseId = 0;
	my %BaseH = ();
	my %BaseIndex = ();
	for my $i (0 .. $#Start)
	{
		for my $j ($Start[$i] .. $End[$i])
		{
			die "[ Error ] The number of nucleotides not equal with bed sum.\n" if($BaseId > $#Base);
			$BaseH{$j} = $Base[$BaseId];
			$BaseId ++;
			# 坐标转换;
			$BaseIndex{$j} = $BaseId;
		}
	}
	
	my (@OriStart,@OriEnd,@Alt) = ();
	for my $i (0 .. $#Var)
	{
		my ($tChr,$tFrom,$tRef,$tAlt) = split /,/, $Var[$i];
		die "[ Error ] chr in var string not euqal with the nucleotides ($tChr vs. $Chr).\n" if($tChr ne $Chr);
		my $tTo = $tFrom;
		# snv或者insertion，from都等于to，只有deletion时需要处理;
		$tTo = $tFrom + length($tRef) - 1 if(length($tRef) > 1);
		# 确定ref和坐标是否能够对上;
		die "[ Error ] Ref not correct ($tChr,$tFrom,$tTo,$tRef)\n" unless(RefConfirm($RefGen,$tChr,$tFrom,$tTo,$tRef));
		unless($tRef eq $tAlt)
		{
			# 简化;
			($tChr,$tFrom,$tTo,$tRef,$tAlt) = VarSimplify($tChr,$tFrom,$tTo,$tRef,$tAlt);
			# 变异坐标标准化，5'或者3'转移;
			($tChr,$tFrom,$tTo,$tRef,$tAlt) = VarUniform($RefGen,$tChr,$tFrom,$tTo,$tRef,$tAlt,$Ori3Flag);
		}
		
		# 变异有可能与读码框只是有交集，但暂时不考虑，因为这需要额外考虑外显子边界的问题;
		die "[ Error ] $Var[$i] ($tFrom,$BaseIndex{$tFrom}   $tTo,$BaseIndex{$tTo}) not in $NucString\n" unless($BaseIndex{$tFrom} && $BaseIndex{$tTo});
		$tAlt = $Base[$BaseIndex{$tFrom} - 1] . $tAlt if($tRef eq "-" || ! $tRef);
		# 1-based的坐标转换成0-based;
		push @OriStart, $BaseIndex{$tFrom} - 1;
		push @OriEnd, $BaseIndex{$tTo} - 1;
		push @Alt, $tAlt;
	}
	
	
	my %RevisedH = ();
	for my $i (0 .. $#Alt)
	{
		# 将ref对应的碱基替换成ref+alt，并且排除变异坐标交叉的情况;
		for my $j ($OriStart[$i] .. $OriEnd[$i])
		{
			die "[ Error ] Var intersection.\n" if($RevisedH{$j});
			$Base[$j] = "";
			$RevisedH{$j} = 1;
		}
		
		if($Alt[$i] ne "-" && $Alt[$i])
		{
			$Base[$OriStart[$i]] = $Alt[$i];
		}
	}
	# 统计alt在更新后序列里的相对坐标，为了便于后续的人工核查，定为1-based;
	my (@RevisedStart,@RevisedEnd) = ();
	for my $i (0 .. $#Alt)
	{
		# RevisedStart和RevisedEnd分别记录左侧和右侧最靠近的没有改变的坐标;
		$RevisedStart[$i] = 0;
		$RevisedEnd[$i] = 0;
		
		if($Ori3Flag)
		{
			for(my $j = $#Base; $j > $OriEnd[$i]; $j --)
			{
				$RevisedStart[$i] += length($Base[$j]);
			}
		}
		else
		{
			for my $j (0 .. $OriStart[$i] - 1)
			{
				$RevisedStart[$i] += length($Base[$j]);
			}
		}
		
		my $DiffLen = 0;
		for my $j ($OriStart[$i] .. $OriEnd[$i])
		{
			$DiffLen += length($Base[$j]);
		}
		$RevisedEnd[$i] = $RevisedStart[$i] + $DiffLen + 1;
	}
	my $NewNucString = join("",@Base);
	
	return $NewNucString,\@RevisedStart,\@RevisedEnd;
}

# 计算refence指定窗口序列的GC含量;
sub GCPercentOfReference
{
	my ($RefGen,$Chr,$Pos,$WindowSize,$BedtoolsBin) = @_;
	$BedtoolsBin = "bedtools" unless($BedtoolsBin);
	my $GCPercent = "-";
	die "[ Error ] Inappropriate window size ($WindowSize)\n" unless($WindowSize > 0 && $WindowSize < 100000);
	
	my $From = int(($Pos - 1 ) / $WindowSize) * $WindowSize + 1;
	my $To = $From + $WindowSize - 1;
	my $Ref = &RefGet($RefGen,$Chr,$From,$To,$BedtoolsBin);
	my $GCNum = 0;
	my @Base = split //, $Ref;
	for my $i (0 .. $#Base)
	{
		$GCNum ++ if($Base[$i] eq "G" || $Base[$i] eq "C");
	}
	$GCPercent = $GCNum / $WindowSize;
	
	return $GCPercent;
}

# 根据Id信息过滤Fastq条目;
sub FqFilterById
{
	my %IdH = %{$_[0]};
	my $OriFq = $_[1];
	my $FltFq = $_[2];
	# 是否保留指定Id的reads信息，否则过滤掉;
	my $RevFlag = $_[3];
	
	Parameter::BinList::IfFileExist($OriFq);
	my $DirName = dirname $FltFq;
	$DirName = Parameter::BinList::IfDirExist($DirName);
	
	open(ORI,"cat $OriFq |") or die $! unless($OriFq =~ /\.gz$/);
	open(ORI,"zcat $OriFq |") or die $! if($OriFq =~ /\.gz$/);
	open(FLT,"> $FltFq") or die $! unless($FltFq =~ /\.gz$/);
	open(FLT,"| gzip > $FltFq |") or die $! if($FltFq =~ /\.gz$/);
	while(my $IdLine = <ORI>)
	{
		my $BaseLine = <ORI>;
		my $PlusLine = <ORI>;
		my $QualLine = <ORI>;
		
		chomp $IdLine;
		my @Cols = split / /, $IdLine;
		if($IdH{$Cols[0]})
		{
			if($RevFlag)
			{
				print FLT $IdLine,"\n";
				print FLT $BaseLine;
				print FLT $PlusLine;
				print FLT $QualLine;
			}
			else
			{
				next;
			}
		}
		elsif(!$RevFlag)
		{
			print FLT $IdLine,"\n";
			print FLT $BaseLine;
			print FLT $PlusLine;
			print FLT $QualLine;
		}
	}
	close ORI;
	close FLT;
	
	return 1;
}

1;