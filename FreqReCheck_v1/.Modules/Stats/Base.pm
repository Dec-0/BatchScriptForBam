# Package Name
package Stats::Base;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(MeanAndSd NormalizedMeanAndSd CI4Index GType2Freq PercentileValue PercentileValueOnHash StrandBiasOfGATKSB);
use Sort::PureNum;

sub MeanAndSd
{
	my @NumString = @{$_[0]};
	
	my $Sum = 0;
	for $i (0 .. $#NumString)
	{
		$Sum += $NumString[$i];
	}
	my $TotalNum = @NumString;
	die "[ Error ] The number string for mean and sd calculation was empty.\n" unless($TotalNum > 0);
	my $Mean = $Sum / $TotalNum;
	
	$Sum = 0;
	for $i (0 .. $#NumString)
	{
		$Sum += ($NumString[$i] - $Mean) ** 2;
	}
	$Sum = $Sum / $#NumString if($#NumString);
	my $Sd = sqrt($Sum);
	
	return $Mean,$Sd;
}

sub NormalizedMeanAndSd
{
	my @NumString = @{$_[0]};
	my $MultiNum = $_[1];
	$MultiNum = 0 unless($MultiNum);
	
	my $Sum = 0;
	for $i (0 .. $#NumString)
	{
		$Sum += $NumString[$i];
	}
	my $TotalNum = @NumString;
	die "[ Error ] The number string for mean and sd calculation was empty.\n" unless($TotalNum > 0);
	my $Mean = $Sum / $TotalNum;
	
	$Sum = 0;
	for $i (0 .. $#NumString)
	{
		$Sum += ($NumString[$i] - $Mean) ** 2;
	}
	$Sum = $Sum / $#NumString if($#NumString);
	my $Sd = sqrt($Sum);
	$Sd = $Sd / $Mean unless($MultiNum);
	$Sd = $Sd / $MultiNum if($MultiNum);
	$Mean = 1 unless($MultiNum);
	$Mean = $Mean / $MultiNum if($MultiNum);
	
	return $Mean,$Sd;
}

# 用于计算比如Sensitivity、PPV对应95%置信区间的分布范围;
sub CI4Index
{
	my ($Index,$PositiveCaseNum,$CI,$FormatFlag) = @_;
	my ($Min,$Max) = ($Index,$Index);
	$CI = 0.95 unless($CI);
	
	die "[ Error ] Index ($Index) exceed 0 or 1.\n" if($Index < 0 || $Index > 1);
	die "[ Error ] The number of positive cases should not be zero.\n" if($PositiveCaseNum <= 0);
	die "[ Error ] Only CI in (0.9,0.95,0.99) will be accepted.\n" unless($CI == 0.90 || $CI == 0.95 || $CI == 0.99);
	my $tValue = 1 - $Index;
	$tValue = $tValue * $Index;
	$tValue = $tValue / $PositiveCaseNum;
	my $StdError = sqrt($tValue);
	# 90%、95%、99%分别对应1.64、1.96、2.58;
	my $Diff = $StdError * 1.96;
	$Diff = $StdError * 1.64 if($CI == 0.90);
	$Diff = $StdError * 2.58 if($CI == 0.99);
	$Min = $Index - $Diff;
	$Max = $Index + $Diff;
	
	if($FormatFlag)
	{
		$Index = sprintf("%.1f",$Index * 100) . "%";
		$Min = sprintf("%.1f",$Min * 100) . "%";
		$Max = sprintf("%.1f",$Max * 100) . "%";
		$CI = $CI * 100;
		$Index = $Index . "($CI%CI," . $Min . "-" . $Max . ")";
		return $Index;
	}
	
	return $Min,$Max;
}

# 用于通过基因型（母本和子代用“-”分隔）及Fetal Fraction来计算理论频率;
# 计算的是B所代表的频率;
sub GType2Freq
{
	my ($GTypePair,$FF) = @_;
	
	die "[ Error ] Fetal fraction should be in [0,1].\n" unless($FF >= 0 && $FF <= 1);
	my @GType = split /-/, $GTypePair;
	my @Per = (1 - $FF,$FF);
	my $Freq = 0;
	for my $i (0 .. $#GType)
	{
		my @Base = split //, $GType[$i];
		my $BNum = 0;
		for my $j (0 .. $#Base)
		{
			$BNum ++ if($Base[$j] eq "B" || $Base[$j] eq "b");
		}
		$Freq += $Per[$i] * $BNum / length($GType[$i]);
	}
	
	return $Freq;
}

# 用于取得数组一系列的分位值，比如1/2、1/4、0.01等;
sub PercentileValue
{
	my $ParaNum = @_;
	die "[ Error ] There should be 2 parameters for PercentileValue\n" unless($ParaNum == 2);
	# 默认是排序过的;
	my @Items = @{$_[0]};
	my @Per = @{$_[1]};
	for my $i (0 .. $#Per)
	{
		die "[ Error ] Percentile should between 0~1\n" unless($Per[$i] >= 0 && $Per[$i] <= 1);
	}
	@Per = @{PureNumSort(\@Per)};
	my $MaxPerId = $#Per;
	
	my @Value = ();
	my $Total = @Items;
	my ($AccumNum,$PerId) = (0,0);
	for my $i (0 .. $#Items)
	{
		$AccumNum ++;
		my $tPer = $AccumNum / $Total;
		while($tPer >= $Per[$PerId] && $PerId <= $MaxPerId)
		{
			push @Value, $Items[$i];
			$PerId ++;
		}
	}
	for my $i ($PerId .. $MaxPerId)
	{
		push @Value, $Items[-1];
	}
	
	return \@Value;
}

# 用于取得数组一系列的分位值，比如1/2、1/4、0.01等;
sub PercentileValueOnHash
{
	my $ParaNum = @_;
	die "[ Error ] There should be 2 parameters for PercentileValue\n" unless($ParaNum == 2);
	# 默认是排序过的;
	my %Items = %{$_[0]};
	my @Per = @{$_[1]};
	for my $i (0 .. $#Per)
	{
		die "[ Error ] Percentile should between 0~1\n" unless($Per[$i] >= 0 && $Per[$i] <= 1);
	}
	@Per = @{PureNumSort(\@Per)};
	my $MaxPerId = $#Per;
	
	my @SizeD = keys %Items;
	@SizeD = @{PureNumSort(\@SizeD)};
	
	my @Value = ();
	my $Total = 0;
	for my $i (0 .. $#SizeD)
	{
		$Total += $Items{$SizeD[$i]};
	}
	my ($AccumNum,$PerId) = (0,0);
	for my $i (0 .. $#SizeD)
	{
		$AccumNum += $Items{$SizeD[$i]};
		my $tPer = $AccumNum / $Total;
		while($tPer >= $Per[$PerId] && $PerId <= $MaxPerId)
		{
			push @Value, $SizeD[$i];
			$PerId ++;
		}
	}
	for my $i ($PerId .. $MaxPerId)
	{
		push @Value, $SizeD[-1];
	}
	
	return \@Value;
}

# 计算链偏（GATK-SB方法）；
sub StrandBiasOfGATKSB
{
	my ($a,$b,$c,$d) = @_;
	$a = 0 unless($a);
	$b = 0 unless($b);
	$c = 0 unless($c);
	$d = 0 unless($d);
	if($a + $b == 0 || $c + $d == 0)
	{
		return "-";
	}
	
	my $SB1 = ($b / ($a + $b)) * ($c / ($c + $d));
	$SB1 = $SB1 / (($a + $c) / ($a + $b + $c + $d));
	
	my $SB2 = ($a / ($a + $b)) * ($d / ($c + $d));
	$SB2 = $SB2 / (($a + $c) / ($a + $b + $c + $d));
	
	my $SB = $SB1;
	$SB = $SB2 if($SB2 > $SB);
	
	return $SB;
}

1;