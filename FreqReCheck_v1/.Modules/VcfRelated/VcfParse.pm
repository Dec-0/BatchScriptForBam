# Package Name
package VcfRelated::VcfParse;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(VarSimplify VarUniform RefConfirm VarTrans ItemValueGet VarSplit SPColId VcfSingleVarInfoGet VcfSingleVarInfoGetFromLine VcfFormat);
use SeqRelated::Seq;

# 去除多余的碱基;
sub VarSimplify
{
	my ($Chr,$From,$To,$Ref,$Alt) = @_;
	
	if($Ref =~ /^$Alt/ && $Alt ne "*")
	{
		# deletion;
		$From += length($Alt);
		$Ref =~ s/^$Alt//;
		$Alt = "-";
	}
	elsif($Alt =~ /^$Ref/ && $Ref ne "*")
	{
		# insertion;
		$From = $To;
		$Alt =~ s/^$Ref//;
		$Ref = "-";
	}
	$Ref = uc($Ref);
	$Alt = uc($Alt);
	
	return $Chr,$From,$To,$Ref,$Alt;
}

# 将所有的变异都尽量往5'或者3'移动;
sub VarUniform
{
	my ($RefGen,$Chr,$From,$To,$Ref,$Alt,$Ori3Flag,$BedtoolsBin) = @_;
	die "[ Error ] Asterisk * found in Ref or Alt ($Chr,$From,$To,$Ref,$Alt) when var uniform.\n" if($Alt eq "*" || $Ref eq "*");
	
	($Chr,$From,$To,$Ref,$Alt) = &VarSimplify($Chr,$From,$To,$Ref,$Alt);
	if($Ref && $Alt eq "-")
	{
		# deletion;
		if($Ori3Flag)
		{
			# 3’端标准化;
			my $LeftBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$From,$From,$BedtoolsBin);
			my $RightBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$To + 1,$To + 1,$BedtoolsBin);
			while($LeftBase eq $RightBase)
			{
				$From ++;
				$To ++;
				$LeftBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$From,$From,$BedtoolsBin);
				$RightBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$To + 1,$To + 1,$BedtoolsBin);
			}
		}
		else
		{
			# 5’端标准化;
			my $LeftBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$From - 1,$From - 1,$BedtoolsBin);
			my $RightBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$To,$To,$BedtoolsBin);
			while($LeftBase eq $RightBase)
			{
				$From --;
				$To --;
				$LeftBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$From - 1,$From - 1,$BedtoolsBin);
				$RightBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$To,$To,$BedtoolsBin);
			}
		}
		$Ref = SeqRelated::Seq::RefGet($RefGen,$Chr,$From,$To,$BedtoolsBin);
	}
	elsif($Ref eq "-" && $Alt)
	{
		# insertion;
		my $AltLen = length($Alt);
		my @AltBase = split //, $Alt;
		if($Ori3Flag)
		{
			my $tPos = $From;
			my $tBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$tPos + 1,$tPos + 1,$BedtoolsBin);
			my $tNum = 0;
			my $tId = $tNum % $AltLen;
			while($tBase eq $AltBase[$tId])
			{
				$tPos ++;
				$tBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$tPos + 1,$tPos + 1,$BedtoolsBin);
				$tNum ++;
				$tId = $tNum % $AltLen;
			}
			
			if($tPos > $From)
			{
				$tBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$From + 1,$tPos,$BedtoolsBin);
				$tBase = $Alt . $tBase;
				$Alt = substr($tBase,$tPos - $From);
				$From = $tPos;
				$To = $From;
			}
		}
		else
		{
			my $tPos = $From;
			my $tBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$tPos,$tPos,$BedtoolsBin);
			my $tNum = 0;
			my $tId = $#AltBase - ($tNum % $AltLen);
			while($tBase eq $AltBase[$tId])
			{
				$tPos --;
				$tBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$tPos,$tPos,$BedtoolsBin);
				$tNum ++;
				$tId = $#AltBase - ($tNum % $AltLen);
			}
			
			if($tPos < $From)
			{
				$tBase = SeqRelated::Seq::RefGet($RefGen,$Chr,$tPos + 1,$From,$BedtoolsBin);
				$tBase .= $Alt;
				$Alt = substr($tBase,0,$AltLen);
				$From = $tPos;
				$To = $From;
			}
		}
	}
	$Ref = uc($Ref);
	$Alt = uc($Alt);
	
	return $Chr,$From,$To,$Ref,$Alt;
}

# 确认ref序列是否正确;
sub RefConfirm
{
	my ($RefGen,$Chr,$From,$To,$Ref,$BedtoolsBin) = @_;
	my $Flag = 1;
	
	die "[ Error ] The format of coordinate not correct ($Chr,$From,$To,$Ref).\n" if($From eq "-" || $To eq "-");
	die "[ Error in VcfRelated::VcfParse ] To smaller than From ($Chr,$From,$To,$Ref).\n" if($To =~ /\D/ || $From =~ /\D/ || $To < $From);
	if($Ref ne "-")
	{
		$Ref = uc($Ref);
		$Flag = 0 unless(length($Ref) == $To - $From + 1);
		
		if($Flag)
		{
			my $tSeq = SeqRelated::Seq::RefGet($RefGen,$Chr,$From,$To,$BedtoolsBin);
			if($tSeq ne $Ref)
			{
				print "[ Info ] Reference $tSeq not match with $Chr,$From,$To,$Ref\n";
				$Flag = 0;
			}
		}
	}
	
	return $Flag;
}

# 转换成顺反式判断统一的格式;
sub VarTrans
{
	my ($RefGen,$Chr,$From,$To,$Ref,$Alt) = @_;
	
	die "[ Error ] Ref seq not correct for $Chr,$From,$To,$Ref (Reference: $RefGen).\n" unless(&RefConfirm($RefGen,$Chr,$From,$To,$Ref));
	($Chr,$From,$To,$Ref,$Alt) = &VarUniform($RefGen,$Chr,$From,$To,$Ref,$Alt);
	
	$From -- unless($Ref eq "-");
	my $TransVar = join(",",$Chr,$From,$Ref,$Alt);
	
	return $TransVar;
}

# 活动比如AD值;
sub ItemValueGet
{
	my ($ItemLine,$ValueLine,$Item) = @_;
	my $Value = "None";
	
	my @Cols = split /:/, $ItemLine;
	my $Id = -1;
	for my $i (0 .. $#Cols)
	{
		if($Cols[$i] eq $Item)
		{
			$Id = $i;
			last;
		}
	}
	return $Value unless($Id >= 0);
	
	my @VCols = split /:/, $ValueLine;
	$Value = $VCols[$Id];
	
	return $Value;
}

# 将多个alt拆开，并且合并的snv可以被拆开;
sub VarSplit
{
	my $Line = $_[0];
	my @Var = ();
	
	chomp $Line;
	my @Cols = split /\t/, $Line;
	
	if($Cols[4] =~ /,/)
	{
		my @Alt = split /,/, $Cols[4];
		my $AltNum = @Alt;
		
		for my $i (0 .. $#Alt)
		{
			my $tLine = join("\t",@Cols[0 .. 3],$Alt[$i]);
			for my $j (5 .. $#Cols)
			{
				my $tmp = "";
				my $SplitChar = ";";
				$SplitChar = ":" if($Cols[$j] =~ /:/);
				my @Items = split /$SplitChar/, $Cols[$j];
				for my $k (0 .. $#Items)
				{
					my $Name = "";
					$Name = $1 if($Items[$k] =~ s/^([^=]+)=//);
					
					if($Items[$k] =~ /,/)
					{
						my @Num = split /,/, $Items[$k];
						if($#Num == $AltNum - 1)
						{
							$Items[$k] = $Num[$i];
						}
						elsif($#Num == $AltNum)
						{
							$Items[$k] = $Num[0] . "," . $Num[$i + 1];
						}
					}
					elsif($Items[$k] =~ /\//)
					{
						my @Num = split /\//, $Items[$k];
						if($#Num == $AltNum - 1)
						{
							$Items[$k] = $Num[0];
						}
						elsif($#Num == $AltNum)
						{
							$Items[$k] = $Num[0] . "/" . $Num[1];
						}
					}
					$Items[$k] = $Name . "=" . $Items[$k] if($Name);
					
					$tmp .= $SplitChar . $Items[$k];
				}
				$tmp =~ s/^$SplitChar//;
				$tLine .= "\t" . $tmp;
			}
			
			push @Var, $tLine;
		}
	}
	else
	{
		push @Var, $Line;
	}
	
	
	# split the combined snvs;
	my @FinalVar = ();
	for my $i (0 .. $#Var)
	{
		my @Cols = split /\t/, $Var[$i];
		if(length($Cols[3]) == length($Cols[4]) && length($Cols[3]) > 1 && $Cols[3] !~ /-/ && $Cols[4] !~ /-/)
		{
			my @RefItem = split //, $Cols[3];
			my @AltItem = split //, $Cols[4];
			
			for my $j (0 .. $#RefItem)
			{
				next if($RefItem[$j] eq $AltItem[$j]);
				
				my $tPos = $Cols[1] + $j;
				my $tLine = join("\t",$Cols[0],$tPos,$Cols[2],$RefItem[$j],$AltItem[$j],@Cols[5 .. $#Cols]);
				push @FinalVar, $tLine;
			}
		}
		else
		{
			push @FinalVar, $Var[$i];
		}
	}
	
	my %DupFlag = ();
	@Var = ();
	for my $i (0 .. $#FinalVar)
	{
		my @Cols = split /\t/, $FinalVar[$i];
		# 防止同一个变异记录两次;
		my $DupItem = join("\t",@Cols[0 .. 4]);
		next if($DupFlag{$DupItem});
		$DupFlag{$DupItem} = 1;
		
		if($Cols[3] =~ /^$Cols[4]/ && length($Cols[4]) > 1)
		{
			# deletion;
			my $tbase = substr($Cols[4],-1);
			$Cols[1] += length($Cols[4]) - 1;
			$Cols[3] =~ s/^$Cols[4]/$tbase/;
			$Cols[4] = $tbase;
			push @Var, join("\t",@Cols);
		}
		elsif($Cols[4] =~ /^$Cols[3]/ && length($Cols[3]) > 1)
		{
			# insertion;
			my $tbase = substr($Cols[3],-1);
			$Cols[1] += length($Cols[3]) - 1;
			$Cols[4] =~ s/^$Cols[3]/$tbase/;
			$Cols[3] = $tbase;
			push @Var, join("\t",@Cols);
		}
		else
		{
			push @Var, $FinalVar[$i];
		}
	}
	
	return \@Var;
}

# 对于有对照的vcf，定位FORMAT以及样本对应的列;
sub SPColId
{
	my $Vcf = $_[0];
	my $SP = $_[1];
	my $SPbak = $_[2];
	my $SPbak2 = $_[3];
	my ($ColId4Format,$ColId4SP) = (0,0);
	
	my $Return = "";
	$Return = `zcat $Vcf | grep ^# | tail -n1` if($Vcf =~ /\.gz$/);
	$Return = `cat $Vcf | grep ^# | tail -n1` unless($Vcf =~ /\.gz$/);
	chomp $Return;
	my @Cols = split /\t/, $Return;
	for my $i (0 .. $#Cols)
	{
		if($Cols[$i] eq "FORMAT")
		{
			$ColId4Format = $i + 1;
		}
		elsif($Cols[$i] eq $SP)
		{
			$ColId4SP = $i + 1;
			last;
		}
	}
	if(!$ColId4SP && $SPbak)
	{
		for my $i (0 .. $#Cols)
		{
			if($Cols[$i] eq "FORMAT")
			{
				$ColId4Format = $i + 1;
			}
			elsif($Cols[$i] eq $SPbak)
			{
				$ColId4SP = $i + 1;
				last;
			}
		}
	}
	if(!$ColId4SP && $SPbak2)
	{
		for my $i (0 .. $#Cols)
		{
			if($Cols[$i] eq "FORMAT")
			{
				$ColId4Format = $i + 1;
			}
			elsif($Cols[$i] eq $SPbak2)
			{
				$ColId4SP = $i + 1;
				last;
			}
		}
	}
	print "[ Warning ] Cannot find column FORMAT($ColId4Format) and $SP($ColId4SP) in $Vcf\n" unless($ColId4Format && $ColId4SP);
	
	return ($ColId4Format,$ColId4SP);
}

# 从vcf文件中获取变异相关的深度、频率信息;
sub VcfSingleVarInfoGet
{
	my ($File,$Chr,$From,$To,$Ref,$Alt,$RefGen,$BedtoolsBin,$SP,$SPBak) = @_;
	my ($Freq,$FullDp,$RefDp,$AltDp,$RefF,$RefR,$AltF,$AltR,$TLod) = ("-","-","-","-","-","-","-","-","-");
	
	die "[ Error ] Not vcf name ($File)\n" unless($File =~ /\.vcf$/ || $File =~ /\.vcf\.gz$/);
	die "[ Error ] Vcf not exist ($File).\n" unless($File && -e $File);
	
	my ($InfoId,$FormatId,$SPColId) = ("-3","-2","-1");
	if($SP)
	{
		my @Items = split /-/, $SP;
		($FormatId,$SPColId) = &SPColId($File,$SP,$Items[0]) unless($SPBak);
		($FormatId,$SPColId) = &SPColId($File,$SP,$Items[0],$SPBak) if($SPBak);
		$InfoId = $FormatId - 1;
		die "[ Error ] Can not locate columnes for FORMAT and $SP in $File\n" unless($InfoId && $FormatId && $SPColId);
		$InfoId --;
		$FormatId --;
		$SPColId --;
	}
	
	# ref核实;
	$To = $From + length($Ref) - 1 unless($To && $To !~ /\D/);
	$To = $From if($Ref eq "-");
	die "[ Error ] Reference $Ref not match with $Chr,$From,$To\n" unless(&RefConfirm($RefGen,$Chr,$From,$To,$Ref,$BedtoolsBin));
	($Chr,$From,$To,$Ref,$Alt) = &VarUniform($RefGen,$Chr,$From,$To,$Ref,$Alt);
	my $GrepItem = "^'$Chr'\$'\\t'";
	$GrepItem .= "'$From' | grep \$'\\t''$Ref'\$'\\t'" if($Ref ne "-" && length($Ref) == 1 && $Alt ne "-" && length($Alt) == 1);
	open(VCF,"cat $File | grep -v ^# | grep $GrepItem |") or die $! unless($File =~ /\.gz$/);
	open(VCF,"zcat $File | grep -v ^# | grep $GrepItem |") or die $! if($File =~ /\.gz$/);
	while(my $Line = <VCF>)
	{
		chomp $Line;
		my @VarLine = @{&VarSplit($Line)};
		
		my $MatchFlag = 0;
		for my $i (0 .. $#VarLine)
		{
			my @Cols = split /\t/, $VarLine[$i];
			my ($tChr,$tFrom,$tTo,$tRef,$tAlt) = @Cols[0 .. 4];
			next if($tAlt eq "*");
			
			$tTo = $tFrom + length($tRef) - 1;
			($tChr,$tFrom,$tTo,$tRef,$tAlt) = &VarUniform($RefGen,$tChr,$tFrom,$tTo,$tRef,$tAlt);
			
			if($tFrom == $From && $tTo == $To && $tRef eq $Ref && $tAlt eq $Alt)
			{
				if($Cols[$InfoId] =~ /DP4=(\d+),(\d+),(\d+),(\d+)/)
				{
					($RefF,$RefR,$AltF,$AltR) = ($1,$2,$3,$4);
					$RefDp = $RefF + $RefR;
					$AltDp = $AltF + $AltR;
					$FullDp = $RefDp + $AltDp;
					$Freq = sprintf("%.5f",$AltDp / $FullDp) if($FullDp);
				}
				
				unless($RefDp ne "-" && $AltDp ne "-")
				{
					my $Return = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"AD");
					if($Return ne "None")
					{
						($RefDp,$AltDp) = split /,/, $Return;
						$FullDp = $RefDp + $AltDp;
						$Freq = sprintf("%.5f",$AltDp / $FullDp) if($FullDp);
					}
				}
				
				unless($AltF ne "-" && $AltR ne "-")
				{
					my $tF = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"ALT_F1R2");
					my $tR = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"ALT_F2R1");
					($AltF,$AltR) = ($tF,$tR) if($tF ne "None" && $tR ne "None");
				}
				
				unless($RefF ne "-" && $RefR ne "-")
				{
					my $tF = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"REF_F1R2");
					my $tR = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"REF_F2R1");
					($RefF,$RefR) = ($tF,$tR) if($tF ne "None" && $tR ne "None");
				}
				
				unless($RefF ne "-" && $RefR ne "-" && $AltF ne "-" && $AltR ne "-")
				{
					$tF = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"F1R2");
					$tR = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"F2R1");
					if($tF ne "None" && $tR ne "None")
					{
						my @Items = split /,/, $tF;
						($RefF,$AltF) = @Items[0 .. 1];
						@Items = split /,/, $tR;
						($RefR,$AltR) = @Items[0 .. 1];
					}
				}
				
				$RefDp = $RefF + $RefR if($RefF ne "-" && $RefR ne "-" && $RefDp eq "-");
				$AltDp = $AltF + $AltR if($AltF ne "-" && $AltR ne "-" && $AltDp eq "-");
				$FullDp = $RefDp + $AltDp if($RefDp ne "-" && $AltDp ne "-" && $FullDp eq "-");
				if($AltDp ne "-" && $FullDp ne "-" && $Freq eq "-")
				{
					$Freq = sprintf("%.5f",$AltDp / $FullDp) if($FullDp);
				}
				
				# lod;
				if($Cols[$InfoId] =~ /TLOD=([\.\d]+)/ && $Freq ne "-" && $Freq > 0)
				{
					$TLod = $1;
					# this method came from "/annoroad/data1/bioinfo/PROJECT/Commercial/Medical/Leukemia/bin/Leu_IVD/subscript/lod_z.pl";
					$TLod = -(log($Freq) / log(2)) * $TLod + 100;
					$TLod = sprintf("%.2f",$TLod);
				}
				
				$MatchFlag = 1;
				last;
			}
		}
		
		last if($MatchFlag);
	}
	close VCF;
	
	return ($Freq,$FullDp,$RefDp,$AltDp,$RefF,$RefR,$AltF,$AltR,$TLod);
}

# 从vcf文件中获取变异相关的深度、频率信息;
sub VcfSingleVarInfoGetFromLine
{
	my ($File,$VarLine,$SP,$SPBak) = @_;
	my ($Freq,$FullDp,$RefDp,$AltDp,$RefF,$RefR,$AltF,$AltR,$TLod) = ("-","-","-","-","-","-","-","-","-");
	
	die "[ Error ] Not vcf name ($File)\n" unless($File =~ /\.vcf$/ || $File =~ /\.vcf\.gz$/);
	die "[ Error ] Vcf not exist ($File).\n" unless($File && -e $File);
	
	chomp $VarLine;
	my @Cols = split /\t/, $VarLine;
	my ($Chr,$From,$To,$Ref,$Alt) = @Cols[0 .. 4];
	
	my ($InfoId,$FormatId,$SPColId) = ("-3","-2","-1");
	if($SP)
	{
		my @Items = split /-/, $SP;
		($FormatId,$SPColId) = &SPColId($File,$SP,$Items[0]) unless($SPBak);
		($FormatId,$SPColId) = &SPColId($File,$SP,$Items[0],$SPBak) if($SPBak);
		$InfoId = $FormatId - 1;
		die "[ Error ] Can not locate columnes for FORMAT and $SP in $File\n" unless($InfoId && $FormatId && $SPColId);
		$InfoId --;
		$FormatId --;
		$SPColId --;
	}
	
	# 提取信息;
	if($Cols[$InfoId] =~ /DP4=(\d+),(\d+),(\d+),(\d+)/)
	{
		($RefF,$RefR,$AltF,$AltR) = ($1,$2,$3,$4);
		$RefDp = $RefF + $RefR;
		$AltDp = $AltF + $AltR;
		$FullDp = $RefDp + $AltDp;
		$Freq = sprintf("%.5f",$AltDp / $FullDp) if($FullDp);
	}
	
	unless($RefDp ne "-" && $AltDp ne "-")
	{
		my $Return = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"AD");
		if($Return ne "None")
		{
			($RefDp,$AltDp) = split /,/, $Return;
			$FullDp = $RefDp + $AltDp;
			$Freq = sprintf("%.5f",$AltDp / $FullDp) if($FullDp);
		}
	}
	
	unless($AltF ne "-" && $AltR ne "-")
	{
		my $tF = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"ALT_F1R2");
		my $tR = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"ALT_F2R1");
		($AltF,$AltR) = ($tF,$tR) if($tF ne "None" && $tR ne "None");
	}
	
	unless($RefF ne "-" && $RefR ne "-")
	{
		my $tF = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"REF_F1R2");
		my $tR = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"REF_F2R1");
		($RefF,$RefR) = ($tF,$tR) if($tF ne "None" && $tR ne "None");
	}
	
	unless($RefF ne "-" && $RefR ne "-" && $AltF ne "-" && $AltR ne "-")
	{
		$tF = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"F1R2");
		$tR = &ItemValueGet($Cols[$FormatId],$Cols[$SPColId],"F2R1");
		if($tF ne "None" && $tR ne "None")
		{
			my @Items = split /,/, $tF;
			($RefF,$AltF) = @Items[0 .. 1];
			@Items = split /,/, $tR;
			($RefR,$AltR) = @Items[0 .. 1];
		}
	}
	
	$RefDp = $RefF + $RefR if($RefF ne "-" && $RefR ne "-" && $RefDp eq "-");
	$AltDp = $AltF + $AltR if($AltF ne "-" && $AltR ne "-" && $AltDp eq "-");
	$FullDp = $RefDp + $AltDp if($RefDp ne "-" && $AltDp ne "-" && $FullDp eq "-");
	if($AltDp ne "-" && $FullDp ne "-" && $Freq eq "-")
	{
		$Freq = sprintf("%.5f",$AltDp / $FullDp) if($FullDp);
	}
	
	# lod;
	if($Cols[$InfoId] =~ /TLOD=([\.\d]+)/ && $Freq ne "-" && $Freq > 0)
	{
		$TLod = $1;
	}
	
	return ($Freq,$FullDp,$RefDp,$AltDp,$RefF,$RefR,$AltF,$AltR,$TLod);
}

# 获得vcf文件的版本号;
sub VcfFormat
{
	my $Vcf = $_[0];
	my $VersionNum = 0;
	
	my $HeadLine = "";
	$HeadLine = `zcat $Vcf | head -n 1` if($Vcf =~ /.gz$/);
	$HeadLine = `cat $Vcf | head -n 1` unless($Vcf =~ /.gz$/);
	
	if($HeadLine =~ /VCFv([\d\.]+)$/i)
	{
		$VersionNum = $1;
	}
	else
	{
		print "[ Warning ] Unknown version number ($VersionNum) for vcf: $Vcf.\n";
	}
	
	return $VersionNum;
}

1;